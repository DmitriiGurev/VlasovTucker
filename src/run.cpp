#include <iostream>
#include <string>

#include "mesh.h"
#include "typedefs.h"
#include "plasma_parameters.h"
#include "vtk.h"
#include "velocity_grid.h"
#include "solver.h"
#include "log.h"
#include "timer.h"

using namespace std;

int main(int argc, char *argv[])
{
    string meshFileName = "../data/meshes/fully_periodic_mesh.msh";

    Timer timer;
    timer.StartSection();
    Mesh mesh(meshFileName);
    timer.PrintSectionTime("Mesh initialization");

    cout << mesh.faces.size() << " faces, ";
    cout << mesh.tets.size() << " tets\n";

    int nCells = 21;
    double minVX = -5;
    double maxVX = 5;
    VelocityGrid<Tensor> vGrid({ nCells, 3,    3   },
                               { maxVX,  0.1,  0.1 },
                               { minVX, -0.1, -0.1 });

    PlasmaParameters plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Custom;
    plasmaParams.mass    = 1;
    plasmaParams.charge  = 10;

    PlasmaParameters::MaxwellPDF paramsPDF;
    for (auto tet : mesh.tets)
    {
        Point c = tet->centroid;
        double rho = 10 + 0.2 * sin(1 * c.coords[0] * (2 * pi));
        paramsPDF.physDensity.push_back(rho);
    }
    paramsPDF.temperature = 0;
    paramsPDF.mostProbableV = { 0, 0, 0 };

    cout << "The PDF takes up " << vGrid.nCellsTotal * mesh.tets.size() * 8 / pow(10, 6) <<
        " MB of RAM\n"; 
        
    plasmaParams.SetPDF<PlasmaParameters::MaxwellPDF>(paramsPDF);

    VTK::WriteDistribution("initial_distribution", vGrid, plasmaParams.pdf[500]);
    VTK::WriteCellScalarData("initial_density", mesh, plasmaParams.Density());

    Solver solver(&mesh, &vGrid, &plasmaParams);
    solver.writeStep = 50;
    solver.Solve(10000);
}