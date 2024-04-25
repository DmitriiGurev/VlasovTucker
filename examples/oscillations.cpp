#include <iostream>
#include <string>

#include "header.h"

using namespace VlasovTucker;
using namespace std;

// Tensor format
using TensorType = Tucker;

int main()
{
    Timer timer;
    timer.StartSection();

    string meshFileName = "../data/meshes/fully_periodic_coarse.msh";
    Mesh mesh(meshFileName);
    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{1, 2}, {3, 4}, {5, 6}});
    mesh.Reconstruct(4 * pi);

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    cout << "Mesh size: " << mesh.AverageCellSize() << "\n";
    WriteMeshVTK("mesh", mesh);

    timer.PrintSectionTime("Mesh initialization");

    VelocityGrid vGrid({21, 21, 21}, {-1, -1, -1}, {1, 1, 1});

    ParticleData<TensorType> particleData(&mesh, &vGrid);
    particleData.species = "custom";
    particleData.mass = 1;
    particleData.charge = 1;

    MaxwellPDF paramsPDF;
    // auto rhoFunc = [](const Point& p) { return 1 + 0.01 * sin(0.5 * 1 / sqrt(3) * (p[0] + p[1] + p[2])); };
    auto rhoFunc = [](const Point& p) { return 1 /*+ 0.5 * sin(0.5 * p[0]) * sin(0.5 * p[1]) * sin(0.5 * p[2])*/; };
    paramsPDF.physDensity = move(ScalarField(&mesh, rhoFunc));
    paramsPDF.temperature = 0.0;
    paramsPDF.mostProbableV = {0, 0, 0};

    particleData.SetCompressionError(1e-6);
    particleData.SetMaxwellPDF(paramsPDF);

    WriteDistributionVTK("initial_distribution", vGrid,
        particleData.pdf[mesh.tets.size() / 2].Reconstructed());
    WriteCellScalarDataVTK("initial_density", mesh, particleData.Density());

    Solver<TensorType> solver(&mesh, &vGrid, &particleData);
    // solver.SetSparseSolverType(SparseSolverType::ConjugateGradient);
    solver.writeStep = 10;
    solver.timeStep = 1e-2;
    solver.nIterations = 50 / solver.timeStep;  
    // solver.nIterations = 1000; 
    solver.Solve();
}