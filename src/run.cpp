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

using namespace VlasovTucker;
using namespace std;

int main(int argc, char *argv[])
{
    string meshFileName = "../data/meshes/rectangle_very_fine.msh";

    Timer timer;
    timer.StartSection();
    Mesh mesh(meshFileName);
    timer.PrintSectionTime("Mesh initialization");

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";

    VelocityGrid<Tensor3d> vGrid({11, 3, 3}, {-3, -0.1, -0.1}, {3, 0.1, 0.1});

    PlasmaParameters plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Custom;
    plasmaParams.mass = 1;
    plasmaParams.charge = 10;

    PlasmaParameters::MaxwellPDF paramsPDF;

    auto rhoFunc = [](const Point& p) { return 10 + 0.2 * sin(1 * p.coords[0] * (2 * pi)); };

    paramsPDF.physDensity = move(ScalarField(&mesh, rhoFunc));
    paramsPDF.temperature = 0;
    paramsPDF.mostProbableV = {0, 0, 0};

    plasmaParams.SetPDF<PlasmaParameters::MaxwellPDF>(paramsPDF);

    WriteDistributionVTK("initial_distribution", vGrid, plasmaParams.pdf[mesh.tets.size() / 2]);
    WriteCellScalarDataVTK("initial_density", mesh, plasmaParams.Density());

    Solver solver(&mesh, &vGrid, &plasmaParams);
    solver.comprPrecision = 1e-6;

    solver.writeStep = 10;
    solver.Solve(100000);
}