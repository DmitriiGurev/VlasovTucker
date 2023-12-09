#include "poisson.h"
#include "plasma_parameters.h"
#include "vtk.h"
#include "timer.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

int main() 
{
    Timer timer;

    Mesh mesh("../data/rect_box_1989.msh");
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

    VTK::WriteMesh("mesh", mesh);
    
    cout << "Initialize the solver...\n";
    PoissonSolver solver(&mesh);
    timer.PrintSectionTime("Initializing the solver");

    auto rhoFunc = [](const Point& p) 
    {
        return 2 * 4 * M_PI * M_PI * sin((p.coords[0] + p.coords[1]) * 2 * M_PI);
    };

    vector<double> rho = move(ScalarField(&mesh, rhoFunc));
    VTK::WriteCellScalarData("rho", mesh, rho);

    vector<double> phi = solver.Solve(rho);
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
    VTK::WriteCellScalarData("phi", mesh, phi);
    VTK::WriteCellScalarData("phi_analytical", mesh, phiAnalytical);

    vector<array<double, 3>> field = solver.ElectricField(phi);
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
    VTK::WriteCellVectorData("electric_field", mesh, field);
    VTK::WriteCellVectorData("electric_field_analytical", mesh, fieldAnalytical);
    timer.PrintSectionTime("Writing the results");
}