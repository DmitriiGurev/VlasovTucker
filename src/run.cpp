#include <iostream>
#include <string>

#include "mesh.h"
#include "typedefs.h"
#include "plasma_parameters.h"
#include "vtk.h"
#include "velocity_grid.h"
#include "solver.h"

using namespace std;

int main(/*int argc, char *argv[]*/)
{
    string meshFileName = "../data/meshes/fully_periodic_coarse.msh";
    Mesh mesh(meshFileName);
    cout << mesh.faces.size() << " faces, ";
    cout << mesh.tets.size() << " tets\n";

    int nCellsVX = 9;
    double minVX = -1.0;
    double maxVX = 1.0;
    VelocityGrid<Tensor> vGrid({ nCellsVX, nCellsVX, nCellsVX },
                               { maxVX,    maxVX,    maxVX    },
                               { minVX,    minVX,    minVX    });

    PlasmaParameters plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Electron;
    plasmaParams.mass    = elMass;
    plasmaParams.charge  = -elCharge;

    PlasmaParameters::Maxwell paramsPDF;
    // paramsPDF.physDensity = ConstDensity(&mesh, 1.0e0);
    for (auto tet : mesh.tets)
    {
        Point c = tet->centroid;
        paramsPDF.physDensity.push_back(exp(-pow((c - Point(0.5, 0.5, 0.5)).Abs(), 2) * 10));
    }
    paramsPDF.temperature = 1.0e-9;
    // paramsPDF.temperature = 0.0;
    paramsPDF.averageV    = { 0.0, 0.0, 0.0 };

    cout << "The PDF takes up " <<
        vGrid.nCellsTotal * mesh.tets.size() * 8 / pow(10, 6) << " MB of RAM\n"; 
        
    plasmaParams.SetPDF<PlasmaParameters::Maxwell>(paramsPDF);

    VTK::WriteDistribution("initial_distribution", vGrid, plasmaParams.pdf[500]);
    VTK::WriteCellData("initial_density", mesh, plasmaParams.Density());

    Solver solver(&mesh, &vGrid, &plasmaParams);
    solver.Solve(10000);
}