#include <iostream>
#include <string>

#include "mesh.h"
#include "typedefs.h"
#include "plasma_parameters.h"
#include "vtk.h"
#include "velocity_grid.h"
#include "solver.h"
#include "log.h"

using namespace std;

int main(int argc, char *argv[])
{
    string meshFileName = "../data/meshes/fully_periodic_fine.msh";
    Mesh mesh(meshFileName);
    cout << mesh.faces.size() << " faces, ";
    cout << mesh.tets.size() << " tets\n";

    int nCellsVX = 9;
    double minVX = -1;
    double maxVX = 1;
    VelocityGrid<Tensor> vGrid({ nCellsVX, nCellsVX, nCellsVX },
                               { maxVX,    maxVX,    maxVX    },
                               { minVX,    minVX,    minVX    });

    PlasmaParameters plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Custom;
    plasmaParams.mass    = 1;
    plasmaParams.charge  = 10;

    // PlasmaParameters::MaxwellPDF paramsPDF;
    // paramsPDF.physDensity = ConstDensity(&mesh, 1.0e0);
    // paramsPDF.temperature = 0;
    // paramsPDF.mostProbableV = { 0.5, 0.0, 0.0 };

    PlasmaParameters::MaxwellPDF paramsPDF;
    double sum = 0;
    double vol = 0;
    for (auto tet : mesh.tets)
    {
        Point c = tet->centroid;
        // double rho = exp(-pow((c - Point({0.4, 0.5, 0.5})).Abs(), 2) * 40) +
                    //  exp(-pow((c - Point({0.6, 0.5, 0.5})).Abs(), 2) * 40);

        double rho = ((c - Point({0.4, 0.5, 0.5})).Abs() < 0.2) ? 1 : 0;
        rho += ((c - Point({0.6, 0.5, 0.5})).Abs() < 0.2) ? 1 : 0;

        // double rho =  sin(c.coords[0] * (2 * pi));
        
        // double rho = sin(c.coords[0] * 2 * pi) *
        //              sin(c.coords[1] * 2 * pi) *
        //              sin(c.coords[2] * 2 * pi);

        paramsPDF.physDensity.push_back(rho);
    }

    paramsPDF.temperature = 0;
    paramsPDF.mostProbableV = { 0, 0, 0 };

    cout << "The PDF takes up " << 
        vGrid.nCellsTotal * mesh.tets.size() * 8 / pow(10, 6) << " MB of RAM\n"; 
        
    plasmaParams.SetPDF<PlasmaParameters::MaxwellPDF>(paramsPDF);

    VTK::WriteDistribution("initial_distribution", vGrid, plasmaParams.pdf[500]);
    VTK::WriteCellScalarData("initial_density", mesh, plasmaParams.Density());

    Solver solver(&mesh, &vGrid, &plasmaParams);
    solver.writeStep = 50;
    solver.Solve(10000);
}