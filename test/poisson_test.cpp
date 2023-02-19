#include "poisson.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>

using namespace std;

int main()
{
    Mesh mesh("../data/meshes/test_mesh_fine.msh");
    cout << mesh.points.size() << endl;
    cout << mesh.tets.size() << endl;
    
    PoissonSolver solver(mesh);

    auto rhoFunc = [&](const Point& p)
    {
        return eps0 * M_PI * M_PI * 3.0 * 
               sin(p.x * M_PI) *
               sin(p.y * M_PI) *
               sin(p.z * M_PI);
    };

    vector<double> rho(mesh.tets.size());

    for (int i = 0; i < mesh.tets.size(); i++)
        rho[i] = rhoFunc(mesh.tets[i]->centroid);

    vector<double> phi = solver.Solve(rho);

    double err = 0.0;
    for (int i = 0; i < 500; i++)
    {
        double analytical = sin(mesh.tets[i]->centroid.x * M_PI) * 
                            sin(mesh.tets[i]->centroid.y * M_PI) * 
                            sin(mesh.tets[i]->centroid.z * M_PI);

        double actual = phi[i];

        err += abs(analytical - actual);
    }
    cout << err / 500 << "\n";

//     cout << "start writing\n";
//     ofstream out;
//     out.open("out.vtk");
//     out << "# vtk DataFile Version 2.0\n";
//     out << "Poisson_test\n";
//     out << "ASCII\n";
//     out << "DATASET UNSTRUCTURED_GRID\n";
//     out << "POINTS " << mesh.points.size() <<  " float\n";
//     for (auto p : mesh.points)
//         out << p->x << " " << p->y << " " << p->z << "\n";
//     out << "CELLS " << mesh.tets.size() << " " << mesh.tets.size() * 5 << "\n";
//     for (auto t : mesh.tets)
//         out << 4 << " "
//             << t->points[0]->index << " " 
//             << t->points[1]->index << " "
//             << t->points[2]->index << " "
//             << t->points[3]->index << "\n";
//     out << "\n";
//     out << "CELL_TYPES " << mesh.tets.size() << "\n";
//     for (auto t : mesh.tets)
//         out << 4 << "\n";
//     out.close();
//     cout << "done writing\n";
}