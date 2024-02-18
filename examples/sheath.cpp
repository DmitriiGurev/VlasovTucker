#include <iostream>
#include <string>

#include "header.h"

using namespace VlasovTucker;
using namespace std;

// Tensor format
using TensorType = Full;

int main()
{
    // 0. Set plasma parameters
    double elTemperature = 1 * electronvolt; // [K]
    double elDensity = 1e17; // [m^-3]

    double debyeLength = DebyeLength(elTemperature, elDensity, elCharge);
    cout << "Debye length: " << debyeLength << " [m]\n";

    double plasmaFrequency = PlasmaFrequency(elDensity, elCharge, elMass);
    cout << "Plasma frequency: " << plasmaFrequency << " [rad/s]\n";

    double plasmaT = 2 * pi / plasmaFrequency;
    cout << "Characteristic time: " << plasmaT << " [s]\n";

    // 1. Load the mesh
    string meshFileName = "../data/meshes/rectangle_very_fine.msh";
    Mesh mesh(meshFileName);

    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{3, 4}, {5, 6}});

    // Scale down to 10 Debye lengths
    double scaleFactor = 10 * debyeLength;
    mesh.Reconstruct(scaleFactor);

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    WriteMeshVTK("mesh", mesh);

    // 2. Create a velocity grid
    // Let PDF(maxV) = 1e-6 PDF(0)
    double maxV = sqrt(-log(1e-6) * 2 * boltzConst * elTemperature / elMass);
    cout << "Characteristic velocity range: (" << -maxV << ", " << maxV << ") [m/s]\n";

    VelocityGrid vGrid({25, 7, 7}, {-5 * maxV, -maxV, -maxV}, {5 * maxV, maxV, maxV});

    // 3. Set the initial PDF
    PlasmaParameters<TensorType> plasmaParams(&mesh, &vGrid);
    plasmaParams.species = ParticleType::Electron;
    plasmaParams.mass = elMass;
    plasmaParams.charge = -elCharge;

    // Maxwell distribution
    MaxwellPDF paramsPDF;
    auto rhoFunc = [elDensity](const Point& p) { return elDensity; };
    paramsPDF.physDensity = move(ScalarField(&mesh, rhoFunc));
    paramsPDF.temperature = elTemperature;
    paramsPDF.mostProbableV = {0, 0, 0};

    // Tensor compression error
    plasmaParams.SetCompressionError(1e-6);
    
    plasmaParams.SetMaxwellPDF(paramsPDF);

    // 4. Initialize the solver
    Solver<TensorType> solver(&mesh, &vGrid, &plasmaParams);

    // Set boundary conditions
    // Charged plane (Neumann BC)
    FieldBC fieldBC1;
    fieldBC1.type = FieldBCType::ChargedPlane;
    fieldBC1.chargeDenity = 0;
    solver.SetFieldBC(1, fieldBC1);

    // Zero potential (Dirichlet BC)
    FieldBC fieldBC2;
    fieldBC2.type = FieldBCType::ConstantPotential;
    fieldBC2.potential = 0;
    solver.SetFieldBC(2, fieldBC2);

    // Charging wall (Absorbing BC) 
    ParticleBC<TensorType> particleBC1;
    particleBC1.type = ParticleBCType::AbsorbingWall;
    solver.SetParticleBC(1, particleBC1);

    // Free wall (Free BC)
    ParticleBC<TensorType> particleBC2;
    particleBC2.type = ParticleBCType::Free;
    solver.SetParticleBC(2, particleBC2);

    // 5. Solve
    double timeStep = 1e-4 * plasmaT;
    int nSteps = 100000;

    solver.writeStep = 1000;
    solver.Solve(timeStep, nSteps);
}