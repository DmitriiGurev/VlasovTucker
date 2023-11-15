#include "poisson.h"
#include "vtk.h"
#include "timer.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

int main() 
{
    Timer timer;

    Mesh mesh("../data/meshes/rectangle.msh");
    timer.PrintSectionTime("Reading the mesh");
    cout << mesh.points.size() << endl;
    cout << mesh.tets.size() << endl;

    VTK::WriteMesh("mesh", mesh);
    
    cout << "Initialize the solver...\n";
    PoissonSolver solver(&mesh);
    timer.PrintSectionTime("Initializing the solver");
    auto rhoFunc = [&](const Point& p) 
    {
        // return exp(-pow((p - Point({0.5, 0.5, 0.5})).Abs(), 2) * 10);
        
        // double x = p.coords[0];
        // return 10 + (x < 0.5) ? (-0.25 + x) : (0.75 - x);

        double pi = 3.1415926;
        double rho = sin(p.coords[0] * 2 * pi); /* *
                     sin(p.coords[1] * 2 * pi) *
                     sin(p.coords[2] * 2 * pi); */

        return rho;
    };

    vector<double> rho(mesh.tets.size());

    double sum = 0;
    double vol = 0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        rho[i] = rhoFunc(mesh.tets[i]->centroid);
        sum += rho[i] * mesh.tets[i]->volume;
        vol += mesh.tets[i]->volume;
    }
    cout << sum << "\n";
    cout << vol << "\n";

    vector<double> phi = solver.Solve(rho);
    timer.PrintSectionTime("Solving the system");

    double err = 0.0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        double analytical = sin(mesh.tets[i]->centroid.coords[0] * 2 * M_PI); //* 
                            // sin(mesh.tets[i]->centroid.coords[1] * 2 * M_PI) * 
                            // sin(mesh.tets[i]->centroid.coords[2] * 2 * M_PI);

        double actual = phi[i];
        err += abs(analytical - actual);
    }
    cout << err / mesh.tets.size() << "\n";

    vector<array<double, 3>> field = solver.ElectricField(phi);

    double errDX = 0.0;
    VTK::WriteCellScalarData("rho", mesh, rho);
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        double analyticalDX = 2 * M_PI * cos(mesh.tets[i]->centroid.coords[0] * 2 * M_PI); // * 
                            //   sin(mesh.tets[i]->centroid.coords[1] * 2 * M_PI) * 
                            //   sin(mesh.tets[i]->centroid.coords[2] * 2 * M_PI);
        double actualDX = -field[i][0];
        errDX += abs(analyticalDX - actualDX);

        // cout << analyticalDX << " " << actualDX << "\n";
    }
    cout << errDX / mesh.tets.size() << "\n";

    VTK::WriteCellScalarData("phi", mesh, phi);
    VTK::WriteCellVectorData("electric_field", mesh, field);
    timer.PrintSectionTime("Writing the results");
}