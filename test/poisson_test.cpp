#include "poisson.h"
#include "vtk.h"
#include "timer.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

int main()
{
    Timer timer;

    Mesh mesh("../data/meshes/fully_periodic_fine.msh");
    timer.PrintSectionTime("Reading the mesh");
    cout << mesh.points.size() << endl;
    cout << mesh.tets.size() << endl;
    
    PoissonSolver solver(&mesh);
    timer.PrintSectionTime("Initializing the solver");
    auto rhoFunc = [&](const Point& p)
    {
        return eps0 * 12 * M_PI * M_PI * 
               sin(p.coords[0] * 2 * M_PI) *
               sin(p.coords[1] * 2 * M_PI) *
               sin(p.coords[2] * 2 * M_PI);
    };

    vector<double> rho(mesh.tets.size());

    for (int i = 0; i < mesh.tets.size(); i++)
        rho[i] = rhoFunc(mesh.tets[i]->centroid);

    vector<double> phi = solver.Solve(rho);
    timer.PrintSectionTime("Solving the system");

    double err = 0.0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        double analytical = sin(mesh.tets[i]->centroid.coords[0] * 2 * M_PI) * 
                            sin(mesh.tets[i]->centroid.coords[1] * 2 * M_PI) * 
                            sin(mesh.tets[i]->centroid.coords[2] * 2 * M_PI);

        double actual = phi[i];
        err += abs(analytical - actual);
    }
    cout << err / mesh.tets.size() << "\n";

    vector<array<double, 3>> field = solver.ElectricField(phi);

    double errDX = 0.0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        double analyticalDX = 2 * M_PI * cos(mesh.tets[i]->centroid.coords[0] * 2 * M_PI) * 
                              sin(mesh.tets[i]->centroid.coords[1] * 2 * M_PI) * 
                              sin(mesh.tets[i]->centroid.coords[2] * 2 * M_PI);

        double actualDX = -field[i][0];
        errDX += abs(analyticalDX - actualDX);

        // cout << analyticalDX << " " << actualDX << "\n";
    }
    cout << errDX / mesh.tets.size() << "\n";

    VTK::WriteCellScalarData("phi", mesh, phi);
    VTK::WriteCellVectorData("electric_field", mesh, field);
    timer.PrintSectionTime("Writing the results");
}