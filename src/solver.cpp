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
        }

        int writeStep = 100;
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

            VTK::WriteCellData("solution/density/density_" + to_string(it / writeStep), *_mesh, density);
        }
    }
}