#include <iostream>
#include <string>

#include "header.h"

using namespace VlasovTucker;
using namespace std;

// Tensor format
using TensorType = Full;

// TODO: Update

int main(int argc, char *argv[])
{
    Timer timer;
    timer.StartSection();
    
    string meshFileName = "../data/meshes/rectangle_very_fine.msh";
    Mesh mesh(meshFileName);
    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{1, 2}, {3, 4}, {5, 6}});
    mesh.Reconstruct();

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    WriteMeshVTK("mesh", mesh);

    timer.PrintSectionTime("Mesh initialization");

    VelocityGrid vGrid({11, 11, 11}, {-3, -0.1, -0.1}, {3, 0.1, 0.1});

    PlasmaParameters<TensorType> plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Custom;
    plasmaParams.mass = 1;
    plasmaParams.charge = 10;

    MaxwellPDF paramsPDF;
    auto rhoFunc = [](const Point& p) { return 10 + 0.2 * sin(1 * p.coords[0] * (2 * pi)); };
    paramsPDF.physDensity = move(ScalarField(&mesh, rhoFunc));
    paramsPDF.temperature = 0;
    paramsPDF.mostProbableV = {0, 0, 0};

    plasmaParams.SetCompressionError(1e-6);
    plasmaParams.SetMaxwellPDF(paramsPDF);

    WriteDistributionVTK("initial_distribution", vGrid,
        plasmaParams.pdf[mesh.tets.size() / 2].Reconstructed());
    WriteCellScalarDataVTK("initial_density", mesh, plasmaParams.Density());

    Solver<TensorType> solver(&mesh, &vGrid, &plasmaParams);
    solver.writeStep = 100;
    solver.Solve(1e-4, 100000);
}