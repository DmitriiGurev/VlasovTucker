#include "solver.h"
#include "vtk.h"

using namespace std;

Solver::Solver(const Mesh* mesh,
               const VelocityGrid<Tensor>* velocityGrid,
               PlasmaParameters* plasmaParameters) :
    _mesh(mesh),
    _vGrid(velocityGrid),
    _plasmaParams(plasmaParameters)
{
    _vNormal = vector<array<Tensor, 4>>(_mesh->tets.size());
    _vNormalAbs = vector<array<Tensor, 4>>(_mesh->tets.size());

    for (auto tet : _mesh->tets)
    {
        int tetInd = tet->index;
        for (int i = 0; i < 4; i++)
        {
            Point normal = tet->faces[i]->normal;
            _vNormal[tetInd][i] = normal.x * _vGrid->v[0] +
                                  normal.y * _vGrid->v[1] +
                                  normal.z * _vGrid->v[2]; 

            _vNormalAbs[tetInd][i] = _vNormal[tetInd][i].abs();
        }
    }

    cout << "The normal velocity tensors take up " <<
        _vGrid->nCellsTotal * _mesh->tets.size() * 4 * 8 / pow(10, 6) << " MB of RAM\n"; 
}

void Solver::Solve(int nIterations)
{
    double timeStep = 1.0e-4;

    for (int it = 0; it < nIterations; it++)
    {
        for (auto tet : _mesh->tets)
        {
            int tetInd = tet->index;
            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;
                _plasmaParams->pdf[tet->index] += -(timeStep / tet->volume) *
                    0.5 * tet->faces[i]->area *
                    (_vNormal[tetInd][i] *
                    (_plasmaParams->pdf[adjTetInd] +
                    _plasmaParams->pdf[tetInd]) -
                    _vNormalAbs[tetInd][i] *
                    (_plasmaParams->pdf[adjTetInd] -
                    _plasmaParams->pdf[tetInd]));
            }

            array<double, 3> force = {1.0, 0.0, 0.0};
            
            int n0 = _vGrid->nCells[0];
            int n1 = _vGrid->nCells[1];
            int n2 = _vGrid->nCells[2];
            
            Tensor pdfDer(n0, n1, n2);
            pdfDer.setZero();
            for (int i0 = 0; i0 < n0; i0++)
            {
                for (int i1 = 0; i1 < n1; i1++)
                {
                    for (int i2 = 0; i2 < n2; i2++)
                    {
                        if (i0 == 0)
                        {
                            pdfDer(0, i1, i2) = (_plasmaParams->pdf[tetInd](1, i1, i2) -
                                _plasmaParams->pdf[tetInd](n0 - 1, i1, i2));
                        }
                        else if (i0 == n0 - 1)
                        {
                            pdfDer(n0 - 1, i1, i2) = (_plasmaParams->pdf[tetInd](0, i1, i2) -
                                _plasmaParams->pdf[tetInd](n0 - 2, i1, i2));
                        }
                        else
                        {
                            pdfDer(i0, i1, i2) = (_plasmaParams->pdf[tetInd](i0 + 1, i1, i2) -
                                _plasmaParams->pdf[tetInd](i0 - 1, i1, i2));
                        }
                    }
                }
            }
            _plasmaParams->pdf[tet->index] += -timeStep *
                    force[0] * pdfDer / (2 * _vGrid->step[0]);
        }

        int writeStep = 50;
        if (it % writeStep == 0)
        {
            double nSumm = 0.0;
            vector<double> density = _plasmaParams->Density();
            for (auto tet : _mesh->tets)
                nSumm += density[tet->index];

            cout << it * timeStep << " ";
            cout << nSumm / _mesh->tets.size() << " ";
            cout << it << " ";
            cout << "\n";

            VTK::WriteCellData("solution/density/density_" + to_string(it / writeStep),
                *_mesh, density);
            VTK::WriteDistribution("solution/distribution/distribution_" + to_string(it / writeStep),
                *_vGrid, _plasmaParams->pdf[100]);
        }
    }
}