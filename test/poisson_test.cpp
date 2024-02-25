#include "poisson.h"
#include "plasma_parameters.h"
#include "vtk.h"
#include "timer.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace VlasovTucker;
using namespace std;

int main() 
{
    Timer timer;

    string meshFileName = "../data/meshes/fully_periodic_coarse.msh";
    Mesh mesh(meshFileName);
    mesh.PrintBoundaryLabels();
    mesh.SetPeriodicBounaries({{3, 4}, {5, 6}});
    mesh.Reconstruct();

    cout << mesh.faces.size() << " faces, " << mesh.tets.size() << " tets\n";
    WriteMeshVTK("mesh", mesh);

    timer.PrintSectionTime("Mesh initialization");

    double averageSize = 0;
    for (auto tet : mesh.tets)
        averageSize += pow(tet->volume * 6 * sqrt(2), 1 / 3.);
    averageSize /= mesh.tets.size();
    cout << "Average cell size: " << averageSize << "\n";

    cout << "Initialize the solver...\n";

    PoissonSolver solver(&mesh);
 
    // Set boundary conditions (periodic BCs are set automatically)
    PoissonBC dirichletBC;
    dirichletBC.type = PoissonBCType::Dirichlet;
    dirichletBC.value = 1;

    PoissonBC neumannBC;
    neumannBC.type = PoissonBCType::Neumann;
    // neumannBC.gradient = Point({0, 1, 0});
    neumannBC.normalGrad = 1;

    // solver.SetBC(5, dirichletBC);
    // solver.SetBC(6, neumannBC);

    solver.SetBC(1, PoissonBC({PoissonBCType::Neumann, 0, 2}));
    solver.SetBC(2, PoissonBC({PoissonBCType::Dirichlet, 0, 0}));

    solver.Initialize();

    timer.PrintSectionTime("Initializing the solver");

    auto rhoFunc = [](const Point& p) 
    {
        return -1;
        // return 2 * 4 * pi * pi * sin((p[0] + p[1]) * 2 * pi);
    };

    vector<double> rho = ScalarField(&mesh, rhoFunc);
    WriteCellScalarDataVTK("rho", mesh, rho);

    for (int r = 0; r < 3; r++)
        solver.Solve(rho);

    vector<double> phi = solver.Potential();
    timer.PrintSectionTime("Solving the system");

    // double phiAvg = 0;
    // for (int i = 0; i < mesh.tets.size(); i++)
    //     phiAvg += phi[i];
    // phiAvg /= mesh.tets.size();
    // for (int i = 0; i < mesh.tets.size(); i++)
    //     phi[i] -= phiAvg;

    vector<double> phiAnalytical;
    double err = 0.0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        Point c = mesh.tets[i]->centroid;
        double analytical = sin((c[0] + c[1]) * 2 * M_PI);
        double actual = phi[i];
        err += abs(analytical - actual);

        phiAnalytical.push_back(analytical);
    }
    cout << "Phi error: " << err / mesh.tets.size() << "\n";
    WriteCellScalarDataVTK("phi", mesh, phi);
    WriteCellScalarDataVTK("phi_analytical", mesh, phiAnalytical);

    vector<Vector3d> field = solver.ElectricField();
    vector<Vector3d> fieldAnalytical;
    double errField = 0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        Point c = mesh.tets[i]->centroid;
        double cosine = cos((c[0] + c[1]) * 2 * M_PI);
        Point analytical({1, 1, 0});
        analytical = analytical * cosine * (2 * M_PI);

        Vector3d actual = field[i];
        errField += (analytical - actual).Abs();

        fieldAnalytical.push_back(analytical.coords);
    }
    cout << "Field error: " << errField / mesh.tets.size() << "\n";
    WriteCellVectorDataVTK("electric_field", mesh, field);
    WriteCellVectorDataVTK("electric_field_analytical", mesh, fieldAnalytical);
    timer.PrintSectionTime("Writing the results");
}