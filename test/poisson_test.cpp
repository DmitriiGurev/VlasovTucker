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

    Mesh mesh("../data/meshes/rect_box_11102.msh");
    timer.PrintSectionTime("Reading the mesh");
    cout << mesh.points.size() << endl;
    cout << mesh.tets.size() << endl;

    double averageSize = 0;
    for (auto tet : mesh.tets)
    {
        averageSize += pow(tet->volume * 6 * sqrt(2), 1 / 3.);
    }
    averageSize /= mesh.tets.size();
    cout << "Average size: " << averageSize << "\n";

    WriteMeshVTK("mesh", mesh);
    
    cout << "Initialize the solver...\n";
    PoissonSolver solver(&mesh);
    timer.PrintSectionTime("Initializing the solver");

    auto rhoFunc = [](const Point& p) 
    {
        return 2 * 4 * pi * pi * sin((p.coords[0] + p.coords[1]) * 2 * pi);
    };

    vector<double> rho = ScalarField(&mesh, rhoFunc);
    WriteCellScalarDataVTK("rho", mesh, rho);

    for (int r = 0; r < 3; r++)
        solver.Solve(rho);

    vector<double> phi = solver.Potential();
    timer.PrintSectionTime("Solving the system");

    double phiAvg = 0;
    for (int i = 0; i < mesh.tets.size(); i++)
        phiAvg += phi[i];
    phiAvg /= mesh.tets.size();
    for (int i = 0; i < mesh.tets.size(); i++)
        phi[i] -= phiAvg;

    vector<double> phiAnalytical;
    double err = 0.0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        Point c = mesh.tets[i]->centroid;
        double analytical = sin((c.coords[0] + c.coords[1]) * 2 * M_PI);
        double actual = phi[i];
        err += abs(analytical - actual);

        phiAnalytical.push_back(analytical);
    }
    cout << "Phi error: " << err / mesh.tets.size() << "\n";
    WriteCellScalarDataVTK("phi", mesh, phi);
    WriteCellScalarDataVTK("phi_analytical", mesh, phiAnalytical);

    vector<array<double, 3>> field = solver.ElectricField();
    vector<array<double, 3>> fieldAnalytical;
    double errField = 0;
    for (int i = 0; i < mesh.tets.size(); i++)
    {
        Point c = mesh.tets[i]->centroid;
        double cosine = cos((c.coords[0] + c.coords[1]) * 2 * M_PI);
        Point analytical({1, 1, 0});
        analytical = analytical * cosine * (2 * M_PI);

        Point actual(field[i]);
        errField += (analytical - actual).Abs();

        fieldAnalytical.push_back(analytical.coords);
    }
    cout << "Field error: " << errField / mesh.tets.size() << "\n";
    WriteCellVectorDataVTK("electric_field", mesh, field);
    WriteCellVectorDataVTK("electric_field_analytical", mesh, fieldAnalytical);
    timer.PrintSectionTime("Writing the results");
}