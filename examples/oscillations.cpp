#include <iostream>
#include <string>

#include "header.h"

using namespace VlasovTucker;
using namespace std;

// Tensor format
using TensorType = Full;

int main()
{
    Timer timer;
    timer.StartSection();
    
    string meshFileName = "../data/meshes/rectangle_fine.msh";
    Mesh mesh(meshFileName);
    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{1, 2}, {3, 4}, {5, 6}});
    mesh.Reconstruct();

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    WriteMeshVTK("mesh", mesh);

    timer.PrintSectionTime("Mesh initialization");

    VelocityGrid vGrid({11, 11, 11}, {-3, -0.1, -0.1}, {3, 0.1, 0.1});

    ParticleData<TensorType> particleData(&mesh, &vGrid);
    particleData.species = "custom";
    particleData.mass = 1;
    particleData.charge = 10;

    MaxwellPDF paramsPDF;
    auto rhoFunc = [](const Point& p) { return 10 + 0.2 * sin(1 * p[0] * (2 * pi)); };
    paramsPDF.physDensity = move(ScalarField(&mesh, rhoFunc));
    paramsPDF.temperature = 0;
    paramsPDF.mostProbableV = {0, 0, 0};

    particleData.SetCompressionError(1e-6);
    particleData.SetMaxwellPDF(paramsPDF);

    WriteDistributionVTK("initial_distribution", vGrid,
        particleData.pdf[mesh.tets.size() / 2].Reconstructed());
    WriteCellScalarDataVTK("initial_density", mesh, particleData.Density());

    Solver<TensorType> solver(&mesh, &vGrid, &particleData);
    solver.writeStep = 100;
    solver.timeStep = 1e-4;
    solver.nIterations = 1e5;  
    solver.Solve();
}