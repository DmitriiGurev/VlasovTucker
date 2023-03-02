#include <iostream>
#include <string>

#include "mesh.h"
#include "tensor_type.h"
#include "plasma_parameters.h"
#include "vtk.h"
#include "velocity_grid.h"

using namespace std;

int main(/*int argc, char *argv[]*/)
{
    string meshFileName = "../data/meshes/test_mesh_periodic.msh";
    Mesh mesh(meshFileName);
    cout << mesh.tets.size() << "\n";

    int nCellsVX = 30;
    double minVX = -1.0e6;
    double maxVX = 1.0e6;
    VelocityGrid<Tensor> vGrid({ nCellsVX, nCellsVX, nCellsVX },
                               { maxVX,    maxVX,    maxVX    },
                               { minVX,    minVX,    minVX    });

    PlasmaParameters plasmaParams;
    plasmaParams.pSpecies = ParticleType::Electron;
    plasmaParams.pMass    = elMass;
    plasmaParams.pCharge  = -elCharge;

    PlasmaParameters::UnifomMaxwell paramsPDF;
    paramsPDF.physDensity = 1.0e0;
    paramsPDF.temperature = 1.0 * electronvolt;
    paramsPDF.averageV    = { 0.0, 0.0, 0.0 };

    plasmaParams.SetPDF<PlasmaParameters::UnifomMaxwell>(mesh, vGrid, paramsPDF);

    VTK::WriteDistribution("distribution", vGrid, plasmaParams.pPDF[0]);

    // Solver solver(mesh, vGrid, plasmaParams);
    // Solver::Solution solution = solver.Solve(100, 1.0e-5);
}