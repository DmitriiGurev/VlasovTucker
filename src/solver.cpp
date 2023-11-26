#include "solver.h"
#include "vtk.h"
#include "poisson.h"
#include "smoother.h"

using namespace std;

Solver::Solver(const Mesh* mesh,
               const VelocityGrid<Tensor>* velocityGrid,
               PlasmaParameters* plasmaParameters) :
    _mesh(mesh),
    _vGrid(velocityGrid),
    _plParams(plasmaParameters)
{
    // TODO: Set LogLevel in CMake
    _log = Log(LogLevel::AllText, "solver_");
    
    _ComputeNormalTensors();

    int memorySpent = _vGrid->nCellsTotal * _mesh->tets.size() * 4 * 8 / pow(10, 6);
    _log << "The normal velocity tensors take up " << memorySpent << " MB of RAM\n"; 
}

void Solver::Solve(int nIterations)
{
    // TODO: Calculate it using the Courant number
    double timeStep = 1.0 * 1.0e-4;

    _log << "Initialize the Poisson solver\n";
    PoissonSolver pSolver(_mesh);

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        _log << "\n" << "Iteration #" << it << "\n";
        
        _log << Indent(1) << "Compute the electric field\n";
        vector<double> rho = _plParams->Density();
        for (int i = 0; i < rho.size(); i++) 
            rho[i] *= _plParams->charge;

        // Smooth the charge density field (?)
        Smoother smoother(_mesh);
        smoother.factor = 0.7;
        smoother.nRounds = 5;
        // smoother.SmoothField(rho);

        vector<double> phi = pSolver.Solve(rho);
        vector<array<double, 3>> field = pSolver.ElectricField(phi);

        _log << Indent(1) << "Compute the right-hand side\n";
        vector<Tensor> rhs(_mesh->tets.size());
        for (auto tet : _mesh->tets)
        {
            rhs[tet->index] = Tensor(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
            rhs[tet->index].setZero();
        }

        _log << Indent(2) << "Boltzmann part\n";
        for (auto tet : _mesh->tets)
        {
            int tetInd = tet->index;
            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;
                
                rhs[tet->index] -= 0.5 * tet->faces[i]->area / tet->volume * 
                    (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] + _plParams->pdf[tetInd]) -
                    _vNormalAbs[tetInd][i] * (_plParams->pdf[adjTetInd] - _plParams->pdf[tetInd]));
            }
        }

        _log << Indent(2) << "Vlasov part\n";
        for (auto tet : _mesh->tets)
        {
            // TODO: Add a constant electric field
            for (int i = 0; i < 3; i++)
            {
                double forceComponent = _plParams->charge * field[tet->index][i];
                rhs[tet->index] -= forceComponent * _PDFDerivative(tet, i);
            }
        }

        _log << Indent(1) << "Time integration\n";
        if (timeIntegrationScheme == TimeIntegrationScheme::Explicit)
        {
            _log << Indent(2) << "Scheme: Explicit\n";
            for (auto tet : _mesh->tets)
                _plParams->pdf[tet->index] += timeStep * rhs[tet->index];
        }
        
        if (timeIntegrationScheme == TimeIntegrationScheme::Implicit)
        {
            _log << Indent(2) << "Scheme: Implicit (LU-SGS)\n";
            vector<Tensor> delta = move(rhs);

            Tensor unity(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
            unity.setConstant(1);

            _log << Indent(3) << "Backward sweep\n";
            for (int tetInd = _mesh->tets.size() - 1; tetInd >= 0; tetInd--)
            {
                Tet* tet = _mesh->tets[tetInd];
                for (int i = 0; i < 4; i++)
                {
                    int adjTetInd = tet->adjTets[i]->index;
                    if (adjTetInd > tetInd)
                    {   
                        delta[tetInd] -= 0.5 * tet->faces[i]->area / tet->volume * 
                            (_vNormal[tetInd][i] - _vNormalAbs[tetInd][i]) * delta[adjTetInd];
                    }
                }

                Tensor diagonalTensor = unity / timeStep;
                for (int i = 0; i < 4; i++)
                {
                    diagonalTensor += 0.5 * tet->faces[i]->area / tet->volume *
                        (_vNormal[tetInd][i] + _vNormalAbs[tetInd][i]);
                }

                delta[tetInd] /= diagonalTensor;
            }

            _log << Indent(3) << "Forward sweep\n";
            for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
            {
                Tensor increment(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
                increment.setZero();

                Tet* tet = _mesh->tets[tetInd];
                for (int i = 0; i < 4; i++)
                {
                    int adjTetInd = tet->adjTets[i]->index;
                    if (adjTetInd < tetInd)
                    {   
                        increment -= 0.5 * tet->faces[i]->area / tet->volume * 
                            (_vNormal[tetInd][i] - _vNormalAbs[tetInd][i]) * delta[adjTetInd];
                    }

                    Tensor diagonalTensor = unity / timeStep;
                    for (int i = 0; i < 4; i++)
                    {
                        diagonalTensor += 0.5 * tet->faces[i]->area / tet->volume *
                            (_vNormal[tetInd][i] + _vNormalAbs[tetInd][i]);
                    }

                    delta[tetInd] += increment / diagonalTensor;
                }
            }

            for (auto tet : _mesh->tets)
                _plParams->pdf[tet->index] += delta[tet->index];
        }

        double nSumm = 0.0;
        vector<double> density = _plParams->Density();
        for (auto tet : _mesh->tets)
            nSumm += density[tet->index];

        _log << Indent(1) << "Time: " << it * timeStep << "\n";
        _log << Indent(1) << "Total density: " << nSumm << "\n";

        if (abs(nSumm) > 1e8)
            throw runtime_error("Solution diverged");

        // Write VTK files
        VTK::WriteCellScalarData("potential", *_mesh, phi);
        VTK::WriteCellVectorData("field", *_mesh, field);
        VTK::WriteCellScalarData("density", *_mesh, _plParams->Density());

        if (it % writeStep == 0)
        {
            string densityFile = "solution/density/density_" + to_string(it / writeStep);
            VTK::WriteCellScalarData(densityFile, *_mesh, density);

            string phiFile = "solution/phi/phi_" + to_string(it / writeStep);
            VTK::WriteCellScalarData(phiFile, *_mesh, phi);

            string fieldFile = "solution/field/e_" + to_string(it / writeStep);
            VTK::WriteCellVectorData(fieldFile, *_mesh, field);

            string distrFile = "solution/distribution/distribution_" + to_string(it / writeStep);
            VTK::WriteDistribution(distrFile, *_vGrid, _plParams->pdf[100]);

            // Test
            vector<double> analytical(_mesh->tets.size());
            for (auto tet : _mesh->tets)
            {
                analytical[tet->index] =
                    10 + sin((-2 * it * timeStep + 2 * tet->centroid.coords[0]) * (2 * pi));
            }

            string analyticalFile = "solution/analytical/analytical_" + to_string(it / writeStep);
            VTK::WriteCellScalarData(analyticalFile, *_mesh, analytical);
        }
    }
}

void Solver::_ComputeNormalTensors()
{
    _vNormal = vector<array<Tensor, 4>>(_mesh->tets.size());
    _vNormalAbs = vector<array<Tensor, 4>>(_mesh->tets.size());

    for (auto tet : _mesh->tets)
    {
        int tetInd = tet->index;
        for (int i = 0; i < 4; i++)
        {
            Point normal = tet->faces[i]->normal;
            _vNormal[tetInd][i] = normal.coords[0] * _vGrid->v[0] +
                                  normal.coords[1] * _vGrid->v[1] +
                                  normal.coords[2] * _vGrid->v[2]; 

            _vNormalAbs[tetInd][i] = _vNormal[tetInd][i].abs();
        }
    }
}

Tensor Solver::_PDFDerivative(const Tet* tet, int ind) const
{
    int tetInd = tet->index;

    Tensor pdfDer(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);

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

    return pdfDer;
}