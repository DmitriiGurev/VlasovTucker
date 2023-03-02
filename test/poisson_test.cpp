#include "poisson.h"
#include "vtk.h"
#include "timer.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

int main()
{
    Timer timer;

    Mesh mesh("../data/meshes/test_mesh_mixed.msh");
    timer.PrintSectionTime("Reading the mesh");
    cout << mesh.points.size() << endl;
    cout << mesh.tets.size() << endl;
    
    PoissonSolver solver(mesh);
    timer.PrintSectionTime("Initializing the solver");
    auto rhoFunc = [&](const Point& p)
    {
        return eps0 * (4 + 0.25) * M_PI * M_PI * 
               sin(p.x * 2 * M_PI) *
               sin(p.y * 2 * M_PI) *
               cos(p.z * M_PI * 0.5);
    };

    vector<double> rho(mesh.tets.size());

    for (int i = 0; i < mesh.tets.size(); i++)
        rho[i] = rhoFunc(mesh.tets[i]->centroid);

    map<string, vector<double>> data;
    data["phi"] = solver.Solve(rho);
    timer.PrintSectionTime("Solving the system");

    double err = 0.0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        double analytical = sin(mesh.tets[i]->centroid.x * 2 * M_PI) * 
                            sin(mesh.tets[i]->centroid.y * 2 * M_PI) * 
                            cos(mesh.tets[i]->centroid.z * 0.5 * M_PI);

        double actual = data["phi"][i];
        err += abs(analytical - actual);
    }
    cout << err / mesh.tets.size() << "\n";

    WriteToVTK(ModeVTK::CellData, "out", mesh, data);
    timer.PrintSectionTime("Writing the results");
}