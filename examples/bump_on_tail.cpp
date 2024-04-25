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
    
    string meshFileName = "../data/meshes/landau_damping_tests/rectangle_0.0088.msh";
    Mesh mesh(meshFileName);
    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{1, 2}, {3, 4}, {5, 6}});
    mesh.Reconstruct(2 * pi / 0.3);

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    cout << "Mesh size: " << mesh.AverageCellSize() << "\n";
    WriteMeshVTK("mesh", mesh);

    timer.PrintSectionTime("Mesh initialization");

    VelocityGrid vGrid({401, 3, 3}, {-10, -10, -10}, {10, 10, 10});

    ParticleData<TensorType> particleData(&mesh, &vGrid);
    particleData.species = "custom";
    particleData.mass = 1;
    particleData.charge = 1;

    // MaxwellPDF paramsPDF;
    auto rhoFunc = [](const Point& p) { return 1 + 0.04 * sin(0.3 * p[0]); };
    // paramsPDF.physDensity = move(ScalarField(&mesh, rhoFunc));
    // paramsPDF.temperature = 1.0;
    // paramsPDF.mostProbableV = {0, 0, 0};

    particleData.SetCompressionError(1e-6);
    // particleData.SetMaxwellPDF(paramsPDF);

    particleData.pdf.resize(mesh.tets.size());
    for (int tetInd = 0; tetInd < mesh.tets.size(); tetInd++)
    {
        Point c = mesh.tets[tetInd]->centroid;
        Tensor3d pdf(vGrid.nCells[0], vGrid.nCells[1], vGrid.nCells[2]);
        for (int i = 0; i < vGrid.nCells[0]; i++)
        {
            Vector3d v = vGrid.At(i, vGrid.nCells[1] / 2, vGrid.nCells[2] / 2);
            double vSq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            pdf(i, vGrid.nCells[1] / 2, vGrid.nCells[2] / 2) = rhoFunc(c) * 1 / sqrt(2 * pi) *
                (0.9 * exp(-vSq / 2) + 0.2 * exp(-pow(v[0] - 4.5, 2) / (2 * pow(0.5, 2)))) /
                (vGrid.step[1] * vGrid.step[2]);

        }
        particleData.pdf[tetInd] = TensorType(pdf);
        particleData.pdf[tetInd].Compress(1e-6, vGrid.nCells[0]);
    }


    WriteDistributionVTK("initial_distribution", vGrid,
        particleData.pdf[mesh.tets.size() / 2].Reconstructed());
    WriteCellScalarDataVTK("initial_density", mesh, particleData.Density());

    Solver<TensorType> solver(&mesh, &vGrid, &particleData);
    // solver.SetSparseSolverType(SparseSolverType::ConjugateGradient);
    solver.writeStep = 1000;
    solver.timeStep = 2e-3;
    solver.nIterations = 100 / solver.timeStep + 1;
    // solver.nIterations = 1000; 
    solver.Solve();
}