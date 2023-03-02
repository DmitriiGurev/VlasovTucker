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

    int nCellsVX = 35;
    double minVX = -1.0e6;
    double maxVX = 1.0e6;
    VelocityGrid<Tensor> vGrid({ nCellsVX, nCellsVX, nCellsVX },
                               { maxVX,    maxVX,    maxVX    },
                               { minVX,    minVX,    minVX    });

    PlasmaParameters plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Electron;
    plasmaParams.mass    = elMass;
    plasmaParams.charge  = -elCharge;

    PlasmaParameters::UnifomMaxwell paramsPDF;
    paramsPDF.physDensity = 1.0e0;
    paramsPDF.temperature = 1.0 * electronvolt;
    paramsPDF.averageV    = { 0.0, 0.0, 0.0 };

    cout << "The PDF takes up " <<
        vGrid.nCellsTotal * mesh.tets.size() * 8 / pow(10, 9) << " GB of RAM\n"; 
        
    plasmaParams.SetPDF<PlasmaParameters::UnifomMaxwell>(paramsPDF);

    VTK::WriteDistribution("distribution", vGrid, plasmaParams.pdf[0]);

    // Solver solver(mesh, vGrid, plasmaParams);
    // Solver::Solution solution = solver.Solve(100, 1.0e-5);
}