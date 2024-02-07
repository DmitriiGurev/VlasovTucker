#include "solver.h"
#include "poisson.h"
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

    // Set boundary conditions
    // _boundaryConditions = std::vector<ParticleBC>(_mesh->faces.size());
    
    _PrecomputeNormalTensors();
}

template <typename TensorType>
void Solver<TensorType>::Solve(int nIterations)
{
    // TODO: Calculate it using the Courant number
    double timeStep = 1.0 * 1.0e-4;

    _log << "Initialize the Poisson solver\n";

    Timer timer;
    PoissonSolver pSolver(_mesh);
    timer.PrintSectionTime("Poisson solver initialization");

    double comprErr = _plParams->CompressionError();
    double maxRank = _plParams->MaxRank();

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        Timer timer;
        _log << "\n" << "Iteration #" << it << "\n";
        
        _log << Indent(1) << "Compute the electric field\n";

        vector<double> rho = _plParams->Density();
        for (int i = 0; i < rho.size(); i++) 
            rho[i] *= _plParams->charge;

        pSolver.Solve(rho);
        // vector<double> phi = pSolver.Potential(); // For debugging
        vector<array<double, 3>> field = move(pSolver.ElectricField());
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

                rhs[tetInd] -= tet->faces[i]->area / tet->volume *
                    0.5 * (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] +
                    _plParams->pdf[tetInd]) - _vNormalAbs[tetInd][i] *
                    (_plParams->pdf[adjTetInd] - _plParams->pdf[tetInd]));

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
            _plParams->pdf[tetInd].Compress(comprErr, maxRank);
        }

        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time: " << it * timeStep << "\n";

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

            // string phiFile = "solution/phi/phi_" + to_string(it / writeStep);
            // WriteCellScalarDataVTK(phiFile, *_mesh, phi);

            // string fieldFile = "solution/field/e_" + to_string(it / writeStep);
            // WriteCellVectorDataVTK(fieldFile, *_mesh, field);

            string distrFile = "solution/distribution/distribution_" + to_string(it / writeStep);
            // WriteDistributionVTK(distrFile, *_vGrid, _plParams->pdf[_mesh->tets.size() / 2]);
            WriteDistributionVTK(distrFile, *_vGrid,
                _plParams->pdf[_mesh->tets.size() / 2].Reconstructed());

            // string analytFile = "solution/analytical/analytical_" + to_string(it / writeStep);
            // auto rhoFunc = [it, timeStep](const Point& p)
            // { 
            //     return 10 + 0.2 * sin(1 * (p.coords[0] - it * timeStep) * (2 * pi));
            //     // return 10 + 0.2 * (((p.coords[0] - it * timeStep < 0.2) &&
            //     //     (p.coords[0] - it * timeStep > 0.0)) ? 1 : 0);
            // };
            // WriteCellScalarDataVTK(analytFile, *_mesh, ScalarField(_mesh, rhoFunc));
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
    reduction /= _mesh->tets.size();
    absReduction /= _mesh->tets.size();

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