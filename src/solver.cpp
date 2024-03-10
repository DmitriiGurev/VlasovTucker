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
                           ParticleData<TensorType>* particleData) :
    _mesh(mesh),
    _vGrid(velocityGrid),
    _pData(particleData)
{
    _log = Log(LogLevel::Console);

    // Compute the tensors vNormal and vNormalAbs 
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

// Field boundary conditions
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
        poissonBC.normalGrad = bc.chargeDensity / (2 * epsilon0); 
    }
    _poissonSolver.SetBC(boundaryInd, poissonBC);
}

// Particle boundary conditions
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
void Solver<TensorType>::SetSparseSolverType(SparseSolverType type)
{
    _poissonSolver.SetSparseSolverType(type);
}

template <typename TensorType>
void Solver<TensorType>::Solve()
{
    // Calculate the area of absorbing walls, if any, and set the collected charge to zero
    _InitializeWallCharge();

    _log << "Initialize the Poisson solver\n";
    Timer timer;
    _poissonSolver.Initialize();
    timer.PrintSectionTime("Poisson solver initialization");

    _log << "Start the main loop\n";
    for (int iteration = 0; iteration < nIterations; iteration++)
    {
        _log << "\n" << "Iteration #" << iteration << "\n";
        _log << "Time: " << iteration * timeStep << "\n";
        
        _log << Indent(1) << "Compute the electric field\n";
        // Calculate the charge density
        _rho = _pData->Density();
        for (int i = 0; i < _rho.size(); i++)
        {
            _rho[i] *= _pData->charge;
            // Add the background charge
            if (!backgroundChargeDensity.empty())
                _rho[i] += backgroundChargeDensity[i];
        }

        // Solve the Poisson equation
        _poissonSolver.Solve(_rho);
        _phi = _poissonSolver.Potential();
        _field = _poissonSolver.ElectricField();

        timer.PrintSectionTime(Indent(1) + "Done");

        _log << Indent(1) << "Update the PDF\n";
        _UpdatePDF();
        timer.PrintSectionTime(Indent(1) + "Done");

        _log << Indent(1) << "Update the boundary conditions\n";
        // Set the updated Neumann condition on the absorbing boundaries 
        for (auto pair : _wallCharge)
        {
            int boundaryInd = pair.first;
            double charge = pair.second;
            
            double area = _wallArea[boundaryInd];

            FieldBC fieldBC;
            fieldBC.type = FieldBCType::ChargedPlane;
            fieldBC.chargeDensity = charge / area;

            SetFieldBC(boundaryInd, fieldBC);
        }

        // Write the results
        if (iteration % writeStep == 0)
            _WriteResults(iteration);
    }
}

template <typename TensorType>
void Solver<TensorType>::_UpdatePDF()
{
    // Get the tensor compression parameters
    double comprErr = _pData->CompressionError();
    double maxRank = _pData->MaxRank();

    _log << Indent(2) << "Compute the right-hand side\n";
    // Initialize a zero RHS
    vector<TensorType> rhs;
    for (auto tet : _mesh->tets)
    {
        Tensor3d zero(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
        zero.setZero();
        TensorType tensor(zero);
        rhs.push_back(tensor);
    }

    _log << Indent(3) << "Boltzmann part\n";
    #pragma omp parallel for 
    for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
    {
        Tet* tet = _mesh->tets[tetInd];
        for (int f = 0; f < 4; f++)
        {
            Face* face = tet->faces[f];
            const auto& bc = _faceParticleBC[face->index];
            TensorType flux = _Flux(tet, f, bc.type);
            rhs[tetInd] -= face->area / tet->volume * flux;

            // Update the charge of absorbing planes
            if (bc.type == ParticleBCType::Absorbing && bc.collectCharge)
            {
                double particlesAbsorbed = timeStep * face->area * flux.Sum() *
                    _vGrid->cellVolume;

                #pragma omp critical
                _wallCharge[face->entity] += _pData->charge * particlesAbsorbed;
                    
            }

            // Recompress the RHS tensor
            rhs[tetInd].Compress(comprErr, maxRank);
        }
    }

    _log << Indent(3) << "Vlasov part\n";
    #pragma omp parallel for 
    for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
    {
        Tet* tet = _mesh->tets[tetInd];
        for (int k = 0; k < 3; k++)
        {
            double electricField = _field[tetInd][k] + externalField[k];
            double forceComponent = (_pData->charge / _pData->mass) * electricField;
            rhs[tetInd] -= forceComponent * _PDFDerivative(tet, k);
        }

        // Recompress the RHS tensor
        rhs[tetInd].Compress(comprErr, maxRank);
    }

    _log << Indent(2) << "Time integration\n";

    #pragma omp parallel for 
    for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
    {
        _pData->pdf[tetInd] += timeStep * rhs[tetInd];

        // Recompress the PDF
        _pData->pdf[tetInd].Compress(comprErr, maxRank);
    }
}

template <typename TensorType>
void Solver<TensorType>::_WriteResults(int iteration)
{
    _log << "Particle species: " << _pData->species << "\n";

    // Print information about tensors
    double averageTensorSize = 0;
    for (int i = 0; i < _mesh->tets.size(); i++)
        averageTensorSize += _pData->pdf[i].Size();
    averageTensorSize /= (double)_mesh->tets.size();
    _log << Indent(1) << "Average PDF size = " << averageTensorSize << "\n";
    _log << Indent(1) << "(Uncompressed: " << _vGrid->nCellsTotal << ")\n";

    string prefix = "solution/";
    int fileNumber = iteration / writeStep;
    string postfix = "_" + to_string(fileNumber);

    // Write the density
    string densityFile = prefix + "density/density_" + _pData->species + postfix;
    WriteCellScalarDataVTK(densityFile, *_mesh, _pData->Density());

    // Write the average velocity
    string velocityFile = prefix + "velocity/velocity_" + _pData->species + postfix;
    WriteCellVectorDataVTK(velocityFile, *_mesh, _pData->Velocity());

    // Write the PDF
    int tetInd = _mesh->tets.size() / 2;
    string distrFile = prefix + "distribution/distribution_" + _pData->species + postfix;
    WriteDistributionVTK(distrFile, *_vGrid, _pData->pdf[tetInd].Reconstructed());

    // Write the charge density
    string cDensityFile = prefix + "charge_density/charge_density_" + postfix;
    WriteCellScalarDataVTK(cDensityFile, *_mesh, _rho);

    // Write the electric potential
    string phiFile = prefix + "phi/phi_" + postfix;
    WriteCellScalarDataVTK(phiFile, *_mesh, _phi);

    // Write the electric field
    string fieldFile = prefix + "field/e_" + postfix;
    WriteCellVectorDataVTK(fieldFile, *_mesh, _field);
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
        for (int f = 0; f < 4; f++)
        {
            Point normal = tet->faces[f]->normal;
            Tensor3d vNormal = normal[0] * _vGrid->v[0] +
                               normal[1] * _vGrid->v[1] +
                               normal[2] * _vGrid->v[2]; 

            Tensor3d vNormalAbs = vNormal.abs();

            _vNormal[tetInd][f] = TensorType(vNormal);
            _vNormalAbs[tetInd][f] = TensorType(vNormalAbs);

            _vNormal[tetInd][f].Compress(_pData->CompressionError());
            _vNormalAbs[tetInd][f].Compress(_pData->CompressionError(), 6);

            reduction += vNormal.size() / (double)_vNormal[tetInd][f].Size();
            absReduction += vNormalAbs.size() / (double)_vNormalAbs[tetInd][f].Size();
        }
    }
    reduction /= 4 * _mesh->tets.size();
    absReduction /= 4 * _mesh->tets.size();

    _log << "Normal speed size reduction: " << reduction << " times on average\n";
    _log << "Absolute normal speed size reduction: " << absReduction << " times on average\n";
}

template <typename TensorType>
void Solver<TensorType>::_InitializeWallCharge()
{
    for (auto face : _mesh->faces)
    {
        const auto& bc = _faceParticleBC[face->index];
        if (bc.type == ParticleBCType::Absorbing && bc.collectCharge)
        {
            if (!_wallCharge.count(face->entity))
            {
                _wallCharge[face->entity] = 0;
                _wallArea[face->entity] = 0;
            }
            _wallArea[face->entity] += face->area;
        }
    }
}

template <typename TensorType>
TensorType Solver<TensorType>::_Flux(const Tet* tet, int f, ParticleBCType bcType) const
{
    TensorType flux;

    int tetInd = tet->index;
    int adjTetInd = tet->adjTets[f]->index;

    Face* face = tet->faces[f];

    if (bcType == ParticleBCType::NonBoundary || bcType == ParticleBCType::Periodic)
    {
        flux = 0.5 * (_vNormal[tetInd][f] * (_pData->pdf[adjTetInd] +
            _pData->pdf[tetInd]) - _vNormalAbs[tetInd][f] *
            (_pData->pdf[adjTetInd] - _pData->pdf[tetInd]));
    }
    else if (bcType == ParticleBCType::Absorbing)
    {
        flux = 0.5 * (_vNormal[tetInd][f] * _pData->pdf[tetInd] +
            _vNormalAbs[tetInd][f] * _pData->pdf[tetInd]);
    }
    else if (bcType == ParticleBCType::Source)
    {
        flux = 0.5 * (_vNormal[tetInd][f] * (_faceParticleBC[face->index].sourcePDF +
            _pData->pdf[tetInd]) - _vNormalAbs[tetInd][f] *
            (_faceParticleBC[face->index].sourcePDF - _pData->pdf[tetInd]));
    }
    else if (bcType == ParticleBCType::Free)
    {
        flux = _vNormal[tetInd][f] * _pData->pdf[tetInd];
    }

    return flux;
}

template <>
Tucker Solver<Tucker>::_PDFDerivative(const Tet* tet, int ind) const
{
    int tetInd = tet->index;
 
    Tucker pdfCompr(_pData->pdf[tetInd]);

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
                    _pData->pdf[tetInd](iPlus[0], iPlus[1], iPlus[2]) - 
                    _pData->pdf[tetInd](iMinus[0], iMinus[1], iMinus[2])
                    ) / (2 * _vGrid->step[ind]);
            }
        }
    }
    return Full(pdfDer);
}
}