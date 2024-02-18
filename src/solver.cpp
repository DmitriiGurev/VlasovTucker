#include "solver.h"
#include "smoother.h"
#include "vtk.h"
#include "timer.h"
#include "full.h"
#include "tucker.h"

using namespace std;

namespace VlasovTucker
{
template class Solver<Full>;
template class Solver<Tucker>;

template <typename TensorType>
Solver<TensorType>::Solver(const Mesh* mesh,
                           const VelocityGrid* velocityGrid,
                           PlasmaParameters<TensorType>* plasmaParameters) :
    _mesh(mesh),
    _vGrid(velocityGrid),
    _plParams(plasmaParameters)
{
    // TODO: Set LogLevel in CMake
    _log = Log(LogLevel::AllText, "solver_");

    _PrecomputeNormalTensors();

    _poissonSolver = PoissonSolver(_mesh);

     // Label all faces as non-boundary
    _faceParticleBC = vector<ParticleBC<TensorType>>(_mesh->faces.size());

    // Set periodic BC for PDFs
    ParticleBC<TensorType> periodicBC;
    periodicBC.type = ParticleBCType::Periodic;
    for (auto planePair : _mesh->PeriodicBoundaries())
    {
        for (int boundaryMark : planePair)
            SetParticleBC(boundaryMark, periodicBC);
    }
}

template <typename TensorType>
void Solver<TensorType>::SetFieldBC(int boundaryInd, const FieldBC& bc)
{
    PoissonBC poissonBC;
    if (bc.type == FieldBCType::ConstantPotential)
    {
        poissonBC.type = PoissonBCType::Dirichlet;
        poissonBC.value = bc.potential;
    }
    else if (bc.type == FieldBCType::ChargedPlane)
    {
        poissonBC.type = PoissonBCType::Neumann;
        poissonBC.normalGrad = bc.chargeDenity / (2 * epsilon0); 
    }
    _poissonSolver.SetBC(boundaryInd, poissonBC);
}

template <typename TensorType>
void Solver<TensorType>::SetParticleBC(int boundaryInd, const ParticleBC<TensorType>& bc)
{
    for (int i = 0; i < _mesh->faces.size(); i++)
    {
        Face* face = _mesh->faces[i];
        if (face->entity == boundaryInd)
            _faceParticleBC[i] = bc;
    }
}

template <typename TensorType>
void Solver<TensorType>::Solve(double timeStep, int nIterations)
{
    _log << "Initialize the Poisson solver\n";
    Timer timer;
    _poissonSolver.Initialize();
    timer.PrintSectionTime("Poisson solver initialization");

    for (auto face : _mesh->faces)
    {
        if (_faceParticleBC[face->index].type == ParticleBCType::AbsorbingWall)
        {
            if (!_wallCharge.count(face->entity))
            {
                _wallCharge[face->entity] = 0;
                _wallArea[face->entity] = 0;
            }
            _wallArea[face->entity] += face->area;
        }
    }
    double wallCharge = 0;

    double comprErr = _plParams->CompressionError();
    double maxRank = _plParams->MaxRank();

    vector<double> rhoInitial = _plParams->Density();
    double averageDensity = 0;
    for (int i = 0; i < rhoInitial.size(); i++)
        averageDensity += rhoInitial[i];
    averageDensity /= rhoInitial.size();

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        Timer timer;
        _log << "\n" << "Iteration #" << it << "\n";
        
        _log << Indent(1) << "Compute the electric field\n";

        vector<double> rho = _plParams->Density();
        for (int i = 0; i < rho.size(); i++)
        {
            rho[i] -= averageDensity;
            rho[i] *= _plParams->charge;
        }

        _poissonSolver.Solve(rho);
        vector<double> phi = _poissonSolver.Potential(); // For debugging
        vector<array<double, 3>> field = move(_poissonSolver.ElectricField());
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Update the right-hand side\n";
        vector<TensorType> rhs;
        for (auto tet : _mesh->tets)
        {
            Tensor3d t3d(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
            t3d.setZero();

            TensorType tensor(t3d);
            // tensor.Compress(comprErr, maxRank);
            rhs.push_back(tensor);
        }

        _log << Indent(2) << "Boltzmann part\n";
        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];

            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;

                Face* face = tet->faces[i];

                if (_faceParticleBC[face->index].type == ParticleBCType::NonBoundary || 
                    _faceParticleBC[face->index].type == ParticleBCType::Periodic)
                {
                    rhs[tetInd] -= face->area / tet->volume *
                        0.5 * (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] +
                        _plParams->pdf[tetInd]) - _vNormalAbs[tetInd][i] *
                        (_plParams->pdf[adjTetInd] - _plParams->pdf[tetInd]));
                }
                else if (_faceParticleBC[face->index].type == ParticleBCType::AbsorbingWall)
                {
                    TensorType lost = face->area / tet->volume *
                        0.5 * (_vNormal[tetInd][i] * _plParams->pdf[tetInd] +
                        _vNormalAbs[tetInd][i] * _plParams->pdf[tetInd]);

                    rhs[tetInd] -= lost;

                    #pragma omp critical
                    _wallCharge[face->entity] += tet->volume * _plParams->charge * timeStep * lost.Sum() *
                        _vGrid->cellVolume;
                }
                else if (_faceParticleBC[face->index].type == ParticleBCType::ConstantSource)
                {
                    // Do nothing
                    continue;
                }
                else if (_faceParticleBC[face->index].type == ParticleBCType::Free)
                {
                    rhs[tetInd] -= face->area / tet->volume *
                        _vNormal[tetInd][i] * _plParams->pdf[tetInd];
                }

                rhs[tetInd].Compress(comprErr, maxRank);
            }
        }

        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(2) << "Vlasov part\n";
        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];

            // TODO: Add a constant electric field
            for (int k = 0; k < 3; k++)
            {
                double forceComponent = (_plParams->charge / _plParams->mass) * field[tetInd][k];
                rhs[tetInd] -= forceComponent * _PDFDerivative(tet, k);
            }
            rhs[tetInd].Compress(comprErr, maxRank);
        }
        
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time integration\n";

        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            _plParams->pdf[tetInd] += timeStep * rhs[tetInd];

            for (int i = 0; i < 4; i++)
            {
                Face* face = _mesh->tets[tetInd]->faces[i];
                if (_faceParticleBC[face->index].type == ParticleBCType::ConstantSource)
                {
                    _plParams->pdf[tetInd] = _faceParticleBC[face->index].sourcePDF;
                    break;
                }
            }
            _plParams->pdf[tetInd].Compress(comprErr, maxRank);
        }

        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time: " << it * timeStep << "\n";

        for (auto pair : _wallCharge)
        {
            int boundaryInd = pair.first;
            double charge = pair.second;
            double area = _wallArea[boundaryInd];

            FieldBC fieldBC;
            fieldBC.type = FieldBCType::ChargedPlane;
            fieldBC.chargeDenity = charge / area;
            SetFieldBC(boundaryInd, fieldBC);

            cout << charge << " " << area << "\n";
            cout << "Charge density: " << charge / area << "\n";
        }

        double averageSize = 0;
        for (int i = 0; i < _mesh->tets.size(); i++)
            averageSize += _plParams->pdf[i].Size();
        averageSize /= (double)_mesh->tets.size();
        _log << Indent(1) << "Average PDF size = " << averageSize << "\n";
        _log << Indent(1) << "(Uncompressed: " << _vGrid->nCellsTotal << ")\n";

        if (it % writeStep == 0)
        {
            vector<double> density = _plParams->Density();
            for (auto d : density)
                assert(d == d);

            string densityFile = "solution/density/density_" + to_string(it / writeStep);
            WriteCellScalarDataVTK(densityFile, *_mesh, density);

            string velocityFile = "solution/velocity/velocity_" + to_string(it / writeStep);
            WriteCellVectorDataVTK(velocityFile, *_mesh, _plParams->Velocity());

            string phiFile = "solution/phi/phi_" + to_string(it / writeStep);
            WriteCellScalarDataVTK(phiFile, *_mesh, phi);

            string fieldFile = "solution/field/e_" + to_string(it / writeStep);
            WriteCellVectorDataVTK(fieldFile, *_mesh, field);

            string distrFile = "solution/distribution/distribution_" + to_string(it / writeStep);
            WriteDistributionVTK(distrFile, *_vGrid,
                _plParams->pdf[_mesh->tets.size() / 2].Reconstructed());
        }
    }
}

template <typename TensorType>
void Solver<TensorType>::_PrecomputeNormalTensors()
{
    _vNormal = vector<array<TensorType, 4>>(_mesh->tets.size());
    _vNormalAbs = vector<array<TensorType, 4>>(_mesh->tets.size());

    double reduction = 0;
    double absReduction = 0;

    for (auto tet : _mesh->tets)
    {
        int tetInd = tet->index;
        for (int i = 0; i < 4; i++)
        {
            Point normal = tet->faces[i]->normal;
            Tensor3d vNormal = normal.coords[0] * _vGrid->v[0] +
                               normal.coords[1] * _vGrid->v[1] +
                               normal.coords[2] * _vGrid->v[2]; 

            Tensor3d vNormalAbs = vNormal.abs();

            _vNormal[tetInd][i] = TensorType(vNormal);
            _vNormalAbs[tetInd][i] = TensorType(vNormalAbs);

            _vNormal[tetInd][i].Compress(_plParams->CompressionError());
            _vNormalAbs[tetInd][i].Compress(_plParams->CompressionError(), 6);

            reduction += vNormal.size() / (double)_vNormal[tetInd][i].Size();
            absReduction += vNormalAbs.size() / (double)_vNormalAbs[tetInd][i].Size();
        }
    }
    reduction /= 4 * _mesh->tets.size();
    absReduction /= 4 * _mesh->tets.size();

    _log << "Normal speed size reduction: " << reduction << " times on average\n";
    _log << "Absolute normal speed size reduction: " << absReduction << " times on average\n";
}

template <>
Tucker Solver<Tucker>::_PDFDerivative(const Tet* tet, int ind) const
{
    int tetInd = tet->index;
 
    Tucker pdfCompr(_plParams->pdf[tetInd]);

    auto core = pdfCompr.Core();
    auto u = pdfCompr.U();

    u[ind] = (_vGrid->d[ind] * u[ind]).eval();

    return Tucker(core, u);
}

template <>
Full Solver<Full>::_PDFDerivative(const Tet* tet, int ind) const
{
    int tetInd = tet->index;

    Tensor3d pdfDer(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);

    array<int, 3> i;
    for (i[0] = 0; i[0] < _vGrid->nCells[0]; i[0]++)
    {
        for (i[1] = 0; i[1] < _vGrid->nCells[1]; i[1]++)
        {
            for (i[2] = 0; i[2] < _vGrid->nCells[2]; i[2]++)
            {
                array<int, 3> iPlus = {i[0], i[1], i[2]};
                array<int, 3> iMinus = {i[0], i[1], i[2]};

                if (i[ind] == 0)
                {
                    iPlus[ind] = 1;
                    iMinus[ind] = _vGrid->nCells[ind] - 1;
                }
                else if (i[ind] == _vGrid->nCells[ind] - 1)
                {
                    iPlus[ind] = 0;
                    iMinus[ind] = _vGrid->nCells[ind] - 2;
                }
                else
                {
                    iPlus[ind] += 1;
                    iMinus[ind] -= 1;
                }

                pdfDer(i[0], i[1], i[2]) = (
                    _plParams->pdf[tetInd](iPlus[0], iPlus[1], iPlus[2]) - 
                    _plParams->pdf[tetInd](iMinus[0], iMinus[1], iMinus[2])
                    ) / (2 * _vGrid->step[ind]);
            }
        }
    }
    return Full(pdfDer);
}
}