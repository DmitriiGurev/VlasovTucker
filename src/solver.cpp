#include "solver.h"
#include "poisson.h"
#include "smoother.h"
#include "vtk.h"
#include "timer.h"

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

    Timer timer;
    PoissonSolver pSolver(_mesh);
    timer.PrintSectionTime("Poisson solver initialization");

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        Timer timer;
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
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Update the right-hand side\n";
        vector<Tensor> rhs(_mesh->tets.size());
        for (auto tet : _mesh->tets)
        {
            rhs[tet->index] = Tensor(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
            rhs[tet->index].setZero();
        }

        _log << Indent(2) << "Boltzmann part\n";
        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];

            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;
                
                rhs[tetInd] -= 0.5 * tet->faces[i]->area / tet->volume * 
                    (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] + _plParams->pdf[tetInd]) -
                    _vNormalAbs[tetInd][i] * (_plParams->pdf[adjTetInd] - _plParams->pdf[tetInd]));
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
                double forceComponent = _plParams->charge * field[tetInd][k];
                rhs[tetInd] -= forceComponent * _PDFDerivative(tet, k);
            }
        }
        
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time integration\n";

        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];
            _plParams->pdf[tetInd] += timeStep * rhs[tetInd];
        }

        timer.PrintSectionTime(Indent(2) + "Done");

        double nSumm = 0.0;
        vector<double> density = _plParams->Density();
        for (auto tet : _mesh->tets)
            nSumm += density[tet->index];

        _log << Indent(1) << "Time: " << it * timeStep << "\n";
        _log << Indent(1) << "Total density: " << nSumm << "\n";

        if (abs(nSumm) > 1e8)
            throw runtime_error("Solution diverged");

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