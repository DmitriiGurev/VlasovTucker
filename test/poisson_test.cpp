#include "poisson.h"
#include "particle_data.h"
#include "vtk.h"
#include "timer.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace VlasovTucker;
using namespace std;

vector<string> cubicMeshes = {"../data/meshes/poisson_tests/box_4955_tets.msh",
                              "../data/meshes/poisson_tests/box_10021_tets.msh",
                              "../data/meshes/poisson_tests/box_14324_tets.msh",
                              "../data/meshes/poisson_tests/box_21756_tets.msh",
                              "../data/meshes/poisson_tests/box_37523_tets.msh",
                              "../data/meshes/poisson_tests/box_171081_tets.msh",
                              "../data/meshes/poisson_tests/box_562821_tets.msh"};

vector<string> sphericalMeshes = {"../data/meshes/poisson_tests/sphere_2697_tets.msh",
                                  "../data/meshes/poisson_tests/sphere_5955_tets.msh",
                                  "../data/meshes/poisson_tests/sphere_19938_tets.msh",
                                  "../data/meshes/poisson_tests/sphere_38931_tets.msh",
                                  "../data/meshes/poisson_tests/sphere_90835_tets.msh",
                                  "../data/meshes/poisson_tests/sphere_154054_tets.msh",
                                  "../data/meshes/poisson_tests/sphere_300334_tets.msh"};

void RunTest1()
{
    cout << "Start test #1\n";

    Timer timer;

    ofstream results("poisson_test_1_results.txt");

    for (string meshFileName : cubicMeshes)
    {
        cout << "Mesh: " << meshFileName << "\n";
        
        Mesh mesh(meshFileName);
        mesh.SetPeriodicBounaries({{5, 6}});
        mesh.Reconstruct();
        WriteMeshVTK("mesh", mesh);
        cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
        cout << "Average cell size: " << mesh.AverageCellSize() << "\n";
        timer.PrintSectionTime("Mesh initialization");

        PoissonSolver solver(&mesh);
        solver.SetSparseSolverType(SparseSolverType::ConjugateGradient);

        solver.SetBC(1, PoissonBC({PoissonBCType::Dirichlet, 0, 0}));
        solver.SetBC(2, PoissonBC({PoissonBCType::Dirichlet, 0, 0}));
        solver.SetBC(3, PoissonBC({PoissonBCType::Neumann, 0, 0}));
        solver.SetBC(4, PoissonBC({PoissonBCType::Neumann, 0, 0}));

        solver.Initialize();

        timer.PrintSectionTime("Initializing the solver");

        auto rhoFunc = [](const Point& p) 
        {
            return -epsilon0 * cos(2 * pi * p[1]) *
                (-pow(2 * pi * p[0], 2) + pow(2 * pi, 2) * p[0] + 2);
        };

        vector<double> rho = ScalarField(&mesh, rhoFunc);
        WriteCellScalarDataVTK("rho", mesh, rho);

        int nRounds = 5;
        for (int r = 0; r < nRounds; r++)
            solver.Solve(rho);

        timer.PrintSectionTime("Solving the system");

        vector<double> phi = solver.Potential();
        vector<Vector3d> field = solver.ElectricField();

        WriteCellScalarDataVTK("phi", mesh, phi);
        WriteCellVectorDataVTK("electric_field", mesh, field);

        // Analytical solution
        vector<double> phiA(mesh.tets.size());
        for (auto tet : mesh.tets)
        {
            Point p = tet->centroid;
            phiA[tet->index] = p[0] * (p[0] - 1) * cos(2 * pi * p[1]);
        }
        WriteCellScalarDataVTK("phi_analytical", mesh, phiA);

        vector<Vector3d> fieldA(mesh.tets.size());
        for (auto tet : mesh.tets)
        {
            Point p = tet->centroid;
            fieldA[tet->index][0] = -(2 * p[0] - 1) * cos(2 * pi * p[1]);
            fieldA[tet->index][1] = 2 * pi * (p[0] * p[0] - p[0]) * sin(2 * pi * p[1]);
            fieldA[tet->index][2] = 0;
        }
        WriteCellVectorDataVTK("electric_field_analytical", mesh, fieldA);

        // Calculate the difference
        double phiMSE = 0;
        for (int i = 0; i < mesh.tets.size(); i++)
            phiMSE += pow(phi[i] - phiA[i], 2);
        phiMSE = sqrt(phiMSE / mesh.tets.size());
        cout << "Phi MSE = " << phiMSE << "\n";

        double fieldMSE = 0;
        for (int i = 0; i < mesh.tets.size(); i++)
        {
            fieldMSE += pow(field[i][0] - fieldA[i][0], 2) +
                        pow(field[i][1] - fieldA[i][1], 2) +
                        pow(field[i][2] - fieldA[i][2], 2);
        }
        fieldMSE = sqrt(fieldMSE / mesh.tets.size());
        cout << "Field MSE = " << fieldMSE << "\n";

        results << mesh.AverageCellSize() << " " << phiMSE << " " << fieldMSE << "\n" << flush;
        cout << "\n\n";
    }
}

void RunTest2() 
{
    cout << "Start test #2\n";

    Timer timer;
    
    ofstream results("poisson_test_2_results.txt");

    for (string meshFileName : cubicMeshes)
    {
        cout << "Mesh: " << meshFileName << "\n";

        Mesh mesh(meshFileName);
        mesh.SetPeriodicBounaries({{1, 2}, {3, 4}, {5, 6}});
        mesh.Reconstruct();
        WriteMeshVTK("mesh", mesh);
        cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
        cout << "Average cell size: " << mesh.AverageCellSize() << "\n";
        timer.PrintSectionTime("Mesh initialization");

        PoissonSolver solver(&mesh);
        solver.SetSparseSolverType(SparseSolverType::ConjugateGradient);

        solver.Initialize();

        timer.PrintSectionTime("Initializing the solver");

        auto rhoFunc = [](const Point& p) 
        {
            return epsilon0 * 4 * pi * pi * 6 * sin(2 * pi * (p[0] + p[1] + 2 * p[2]));
        };

        vector<double> rho = ScalarField(&mesh, rhoFunc);
        WriteCellScalarDataVTK("rho", mesh, rho);

        int nRounds = 5;
        for (int r = 0; r < nRounds; r++)
            solver.Solve(rho);

        timer.PrintSectionTime("Solving the system");

        vector<double> phi = solver.Potential();
        vector<Vector3d> field = solver.ElectricField();

        WriteCellScalarDataVTK("phi", mesh, phi);
        WriteCellVectorDataVTK("electric_field", mesh, field);

        // Analytical solution
        vector<double> phiA(mesh.tets.size());
        for (auto tet : mesh.tets)
        {
            Point p = tet->centroid;
            phiA[tet->index] = sin(2 * pi * (p[0] + p[1] + 2 * p[2]));
        }
        WriteCellScalarDataVTK("phi_analytical", mesh, phiA);

        vector<Vector3d> fieldA(mesh.tets.size());
        for (auto tet : mesh.tets)
        {
            Point p = tet->centroid;
            fieldA[tet->index][0] = -2 * pi * 1 * cos(2 * pi * (p[0] + p[1] + 2 * p[2]));
            fieldA[tet->index][1] = -2 * pi * 1 * cos(2 * pi * (p[0] + p[1] + 2 * p[2]));
            fieldA[tet->index][2] = -2 * pi * 2 * cos(2 * pi * (p[0] + p[1] + 2 * p[2]));
        }
        WriteCellVectorDataVTK("electric_field_analytical", mesh, fieldA);

        // Calculate the difference
        double phiMSE = 0;
        for (int i = 0; i < mesh.tets.size(); i++)
            phiMSE += pow(phi[i] - phiA[i], 2);
        phiMSE = sqrt(phiMSE / mesh.tets.size());
        cout << "Phi MSE = " << phiMSE << "\n";

        double fieldMSE = 0;
        for (int i = 0; i < mesh.tets.size(); i++)
        {
            fieldMSE += pow(field[i][0] - fieldA[i][0], 2) +
                        pow(field[i][1] - fieldA[i][1], 2) +
                        pow(field[i][2] - fieldA[i][2], 2);
        }
        fieldMSE = sqrt(fieldMSE / mesh.tets.size());
        cout << "Field MSE = " << fieldMSE << "\n";

        results << mesh.AverageCellSize() << " " << phiMSE << " " << fieldMSE << "\n" << flush;
        cout << "\n\n";
    }
}

void RunTest3() 
{
    cout << "Start test #3\n";

    Timer timer;
    
    ofstream results("poisson_test_3_results.txt");

    for (string meshFileName : sphericalMeshes)
    {
        cout << "Mesh: " << meshFileName << "\n";

        Mesh mesh(meshFileName);
        mesh.Reconstruct();
        WriteMeshVTK("mesh", mesh);
        cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
        cout << "Average cell size: " << mesh.AverageCellSize() << "\n";
        timer.PrintSectionTime("Mesh initialization");

        PoissonSolver solver(&mesh);
        solver.SetSparseSolverType(SparseSolverType::ConjugateGradient);

        solver.SetBC(1, PoissonBC({PoissonBCType::Dirichlet, 0, 0}));

        solver.Initialize();

        timer.PrintSectionTime("Initializing the solver");

        auto rhoFunc = [](const Point& p) 
        {
            return -epsilon0 * pow(p.Abs(), 2);
        };

        vector<double> rho = ScalarField(&mesh, rhoFunc);
        WriteCellScalarDataVTK("rho", mesh, rho);

        int nRounds = 5;
        for (int r = 0; r < nRounds; r++)
            solver.Solve(rho);

        timer.PrintSectionTime("Solving the system");

        vector<double> phi = solver.Potential();
        vector<Vector3d> field = solver.ElectricField();

        WriteCellScalarDataVTK("phi", mesh, phi);
        WriteCellVectorDataVTK("electric_field", mesh, field);

        // Analytical solution
        vector<double> phiA(mesh.tets.size());
        for (auto tet : mesh.tets)
        {
            Point p = tet->centroid;
            phiA[tet->index] = (pow(p.Abs(), 4) - 1) / 20.;
        }
        WriteCellScalarDataVTK("phi_analytical", mesh, phiA);

        vector<Vector3d> fieldA(mesh.tets.size());
        for (auto tet : mesh.tets)
        {
            Point p = tet->centroid;
            fieldA[tet->index][0] = -pow(p.Abs(), 3) / 5. * p[0] / p.Abs();
            fieldA[tet->index][1] = -pow(p.Abs(), 3) / 5. * p[1] / p.Abs();
            fieldA[tet->index][2] = -pow(p.Abs(), 3) / 5. * p[2] / p.Abs();
        }
        WriteCellVectorDataVTK("electric_field_analytical", mesh, fieldA);

        // Calculate the difference
        double phiMSE = 0;
        for (int i = 0; i < mesh.tets.size(); i++)
            phiMSE += pow(phi[i] - phiA[i], 2);
        phiMSE = sqrt(phiMSE / mesh.tets.size());
        cout << "Phi MSE = " << phiMSE << "\n";

        double fieldMSE = 0;
        for (int i = 0; i < mesh.tets.size(); i++)
        {
            fieldMSE += pow(field[i][0] - fieldA[i][0], 2) +
                        pow(field[i][1] - fieldA[i][1], 2) +
                        pow(field[i][2] - fieldA[i][2], 2);
        }
        fieldMSE = sqrt(fieldMSE / mesh.tets.size());
        cout << "Field MSE = " << fieldMSE << "\n";

        results << mesh.AverageCellSize() << " " << phiMSE << " " << fieldMSE << "\n" << flush;
        cout << "\n\n";
    }
}

int main() 
{
    RunTest1();

    RunTest2();
    
    RunTest3();
}