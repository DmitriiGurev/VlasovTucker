#include <iostream>
#include <string>

#include "header.h"

using namespace VlasovTucker;
using namespace std;

// Tensor format
using TensorType = Full;

int main()
{
    // 1. Set plasma parameters
    double elTemperature = 1 * electronvolt; // [K]
    double ionTemperature = 400;             // [K]
    double density = 1e17;                   // [m^-3]
    double ionMass = atomicMass;             // [kg]

    double debyeLength = DebyeLength(elTemperature, density, elCharge);
    cout << "Debye length: " << debyeLength << " [m]\n";

    double plasmaFrequency = PlasmaFrequency(density, elCharge, elMass);
    cout << "Plasma frequency: " << plasmaFrequency << " [rad/s]\n";

    double plasmaT = 2 * pi / plasmaFrequency;
    cout << "Characteristic time: " << plasmaT << " [s]\n";

    // 2. Load the mesh
    string meshFileName = "../data/meshes/rectangle_fine.msh";
    Mesh mesh(meshFileName);

    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{3, 4}, {5, 6}});

    // Scale down to 22 Debye lengths
    double scaleFactor = 22 * debyeLength;
    mesh.Reconstruct(scaleFactor);

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    WriteMeshVTK("mesh", mesh);

    // 3. Create the velocity grids
    // Let PDF(maxV) = 1e-6 PDF(0)
    double maxVE = sqrt(-log(1e-6) * 2 * boltzConst * elTemperature / elMass);
    cout << "Characteristic electron velocity range:" <<
        " (" << -maxVE << ", " << maxVE << ") [m/s]\n";

    VelocityGrid vGridE({50, 5, 5},
                        {-4 * maxVE, -maxVE, -maxVE},
                        {4 * maxVE, maxVE, maxVE});

    double maxVI = sqrt(-log(1e-6) * 2 * boltzConst * ionTemperature / ionMass);
    cout << "Characteristic ion velocity range:" <<
        " (" << -maxVI << ", " << maxVI << ") [m/s]\n";

    VelocityGrid vGridI({50, 5, 5},
                        {-4 * maxVI, -maxVI, -maxVI},
                        {4 * maxVI, maxVI, maxVI});


    // 4. Set the initial PDFs
    auto rhoFunc = [density](const Point& p) { return density; };

    ParticleData<TensorType> particleDataE(&mesh, &vGridE);
    particleDataE.species = "electron";
    particleDataE.mass = elMass;
    particleDataE.charge = -elCharge;

    // Maxwell distribution
    MaxwellPDF maxwellE;
    maxwellE.physDensity = move(ScalarField(&mesh, rhoFunc));
    maxwellE.temperature = elTemperature;
    maxwellE.mostProbableV = {0, 0, 0};

    // Tensor compression error
    particleDataE.SetCompressionError(1e-6);
    particleDataE.SetMaxwellPDF(maxwellE);

    // Same for the ions
    ParticleData<TensorType> particleDataI(&mesh, &vGridI);
    particleDataI.species = "ion";
    particleDataI.mass = ionMass;
    particleDataI.charge = elCharge;

    MaxwellPDF maxwellI;
    maxwellI.physDensity = move(ScalarField(&mesh, rhoFunc));
    maxwellI.temperature = ionTemperature;
    maxwellI.mostProbableV = {0, 0, 0};

    particleDataI.SetCompressionError(1e-6);
    particleDataI.SetMaxwellPDF(maxwellI);

    // 5. Initialize the solver
    Solver<TensorType> solverE(&mesh, &vGridE, &particleDataE);
    Solver<TensorType> solverI(&mesh, &vGridI, &particleDataI);

    // Set boundary conditions
    // Charging wall (Absorbing BC) 
    ParticleBC<TensorType> particleBC1;
    particleBC1.type = ParticleBCType::Absorbing;
    particleBC1.collectCharge = true;
    solverE.SetParticleBC(1, particleBC1);
    solverI.SetParticleBC(1, particleBC1);

    // Free wall (Free BC)
    ParticleBC<TensorType> particleBC2;
    particleBC2.type = ParticleBCType::Free;
    solverE.SetParticleBC(2, particleBC2);
    solverI.SetParticleBC(2, particleBC2);

    // Field BCs can be set just for one solver
    // Charged plane (Neumann BC)
    FieldBC fieldBC1;
    fieldBC1.type = FieldBCType::ChargedPlane;
    // Initial boundary charge density
    fieldBC1.chargeDensity = 0;
    solverE.SetFieldBC(1, fieldBC1);

    // Zero potential (Dirichlet BC)
    FieldBC fieldBC2;
    fieldBC2.type = FieldBCType::ConstantPotential;
    fieldBC2.potential = 0;
    solverE.SetFieldBC(2, fieldBC2);

    // 6. Solve
    MulticomponentSolver<TensorType> multiSolver(&solverE);
    multiSolver.AddSolver(&solverI);

    multiSolver.timeStep = 1e-4 * plasmaT;
    multiSolver.stepMultipliers[&solverI] = 10;
    multiSolver.stepMultipliers[&solverE] = 1;

    multiSolver.nIterations = 1e6;   
    multiSolver.writeStep = 1e3;

    multiSolver.Solve();
}