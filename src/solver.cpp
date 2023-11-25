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
    /// TODO: Set LogLevel in CMake
    _log = Log(LogLevel::AllText, "solver_");
    
    _ComputeNormalTensors();

    _log << "The normal velocity tensors take up " <<
        _vGrid->nCellsTotal * _mesh->tets.size() * 4 * 8 / pow(10, 6) <<
        " MB of RAM\n"; 
}

void Solver::Solve(int nIterations)
{
    /// TODO: Calculate it using the Courant number
    double timeStep = 0.5 * 1.0e-4;

    _log << "Initialize the Poisson solver\n";
    PoissonSolver pSolver(_mesh);

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        _log << "Iteration #" << it << "\n";
        
        _log << "Compute the electric field\n";
        double sum = 0;
        vector<double> rho = _plParams->Density();
        for (int i = 0; i < rho.size(); i++) 
        {
            rho[i] *= _plParams->charge;
            sum += rho[i] * _mesh->tets.at(i)->volume;
        }
        cout << "Total charge = " << sum << "\n";

        // Smoothing (?)
        // Smoother smoother(_mesh);
        // smoother.factor = 0.7;
        // smoother.nRounds = 5;
        // smoother.SmoothField(rho);

        vector<double> phi = pSolver.Solve(rho);
        vector<array<double, 3>> field = pSolver.ElectricField(phi);

        VTK::WriteCellScalarData("potential", *_mesh, phi);
        VTK::WriteCellVectorData("field", *_mesh, field);
        // return;

        _log << "Boltzmann part\n";
        for (auto tet : _mesh->tets)
        {
            int tetInd = tet->index;
            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;
                
                _plParams->pdf[tet->index] += 
                    -(timeStep / tet->volume) * 0.5 * tet->faces[i]->area * 
                    (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] + _plParams->pdf[tetInd]) -
                    _vNormalAbs[tetInd][i] * (_plParams->pdf[adjTetInd] - _plParams->pdf[tetInd]));
            }
        }

        _log << "Vlasov part\n";
        for (auto tet : _mesh->tets)
        {
            array<double, 3> force = field[tet->index];
            for (int i = 0; i < 3; i++)
                force[i] *= _plParams->charge; // ?

            array<Tensor, 3> pdfDers;
            pdfDers[0] = move(_PDFDerivative(tet, 0));
            pdfDers[1] = move(_PDFDerivative(tet, 1));
            pdfDers[2] = move(_PDFDerivative(tet, 2));

            for (int i = 0; i < 3; i++)
                _plParams->pdf[tet->index] += -timeStep * force[i] * pdfDers[i];
        }

        double nSumm = 0.0;
        vector<double> density = _plParams->Density();
        for (auto tet : _mesh->tets)
            nSumm += density[tet->index];

        _log << "Time: " << it * timeStep << "\n";
        _log << "Total density: " << nSumm << "\n";

        if (abs(nSumm) > 1e8)
            throw runtime_error("Solution diverged");

        int writeStep = 5;
        if (it % writeStep == 0)
        {
            VTK::WriteCellScalarData("solution/density/density_" + to_string(it / writeStep), 
                *_mesh, density);
            VTK::WriteCellScalarData("solution/phi/phi_" + to_string(it / writeStep),
                *_mesh, phi);
            VTK::WriteCellVectorData("solution/field/e_" + to_string(it / writeStep),
                *_mesh, field);
            VTK::WriteDistribution("solution/distribution/distribution_" + to_string(it / writeStep),
                *_vGrid, _plParams->pdf[100]);
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

// for (int i0 = 0; i0 < n0; i0++)
// {
//     for (int i1 = 0; i1 < n1; i1++)
//     {
//         for (int i2 = 0; i2 < n2; i2++)
//         {
//             // 0
//             if (i0 == 0)
//             {
//                 pdfDer[0](0, i1, i2) = (_plParams->pdf[tetInd](1, i1, i2) -
//                     _plParams->pdf[tetInd](n0 - 1, i1, i2));
//             }
//             else if (i0 == n0 - 1)
//             {
//                 pdfDer[0](n0 - 1, i1, i2) = (_plParams->pdf[tetInd](0, i1, i2) -
//                     _plParams->pdf[tetInd](n0 - 2, i1, i2));
//             }
//             else
//             {
//                 pdfDer[0](i0, i1, i2) = (_plParams->pdf[tetInd](i0 + 1, i1, i2) -
//                     _plParams->pdf[tetInd](i0 - 1, i1, i2));
//             }

//             // 1
//             if (i1 == 0)
//             {
//                 pdfDer[1](i0, i1, i2) = (_plParams->pdf[tetInd](i0, 1, i2) -
//                     _plParams->pdf[tetInd](i0, n1 - 1, i2));
//             }
//             else if (i1 == n1 - 1)
//             {
//                 pdfDer[1](i0, i1, i2) = (_plParams->pdf[tetInd](i0, 0, i2) -
//                     _plParams->pdf[tetInd](i0, n1 - 2, i2));
//             }
//             else
//             {
//                 pdfDer[1](i0, i1, i2) = (_plParams->pdf[tetInd](i0, i1 + 1, i2) -
//                     _plParams->pdf[tetInd](i0, i1 - 1, i2));
//             }

//             // 2
//             if (i2 == 0)
//             {
//                 pdfDer[2](i0, i1, i2) = (_plParams->pdf[tetInd](i0, i1, 1) -
//                     _plParams->pdf[tetInd](i0, i1, n2 - 1));
//             }
//             else if (i2 == n2 - 1)
//             {
//                 pdfDer[2](i0, i1, i2) = (_plParams->pdf[tetInd](i0, i1, 0) -
//                     _plParams->pdf[tetInd](i0, i1, n2 - 2));
//             }
//             else
//             {
//                 pdfDer[2](i0, i1, i2) = (_plParams->pdf[tetInd](i0, i1, i2 + 1) -
//                     _plParams->pdf[tetInd](i0, i1, i2 - 1));
//             }
//         }
//     }
// }

// // X
// Eigen::MatrixXd der1D = Eigen::MatrixXd::Zero(n0, n0);
// der1D(0, 0) = -2;
// der1D(0, 1) = 1;
// der1D(0, n0 - 1) = 1;

// for (int i = 1; i < n0 - 1; i++)
// {
//     der1D(i, i - 1) = 1;
//     der1D(i, i) = -2;
//     der1D(i, i + 1) = 1;
// }

// der1D(n0 - 1, 0) = 1;
// der1D(n0 - 1, n0 - 2) = 1;
// der1D(n0 - 1, n0 - 1) = -2;

// cout << der1D;