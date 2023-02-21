#include "vtk.h"

using namespace std;

void WriteToVTK(string fileName,
                const Mesh& mesh,
                const map<string, vector<double>>& data)
{
    ofstream out;
    out.open(fileName);

    out << "# vtk DataFile Version 2.0\n";
    out << "Poisson_test\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << mesh.points.size() <<  " float\n";
    for (auto p : mesh.points)
        out << p->x << " " << p->y << " " << p->z << "\n";

    out << "CELLS " << mesh.tets.size() << " " << mesh.tets.size() * 5 << "\n";
    for (auto t : mesh.tets)
        out << 4 << " "
            << t->points[0]->index << " " 
            << t->points[1]->index << " "
            << t->points[2]->index << " "
            << t->points[3]->index << "\n";
    out << "\n";

    out << "CELL_TYPES " << mesh.tets.size() << "\n";
    for (auto t : mesh.tets)
        out << 10 << "\n";
        
    for (const auto& pair : data)
    {
        out << "CELL_DATA " << mesh.tets.size() << "\n";
        out << "SCALARS " << pair.first << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (int i = 0; i < mesh.tets.size(); i++)
            out << pair.second[i] << "\n";
    }

    out.close();
}