#include "vtk.h"

using namespace std;

namespace VlasovTucker
{
void WriteCellScalarDataVTK(string fileName,
                            const Mesh& mesh,
                            const vector<double>& data)
{
    ofstream out;
    out.open(fileName + ".vtk");

    out << "# vtk DataFile Version 2.0\n";
    out << "Vlasov-T\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    
    out << "POINTS " << mesh.points.size() <<  " float\n";
    for (auto p : mesh.points)
        out << p->coords[0] << " " << p->coords[1] << " " << p->coords[2] << "\n";

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

    if (data.size() > 0)
    {
        out << "CELL_DATA " << mesh.tets.size() << "\n";
        out << "SCALARS " << "data" << " double 1\n";
        out << "LOOKUP_TABLE default\n";
        for (int i = 0; i < mesh.tets.size(); i++)
            out << data[i] << "\n";
    }
    
    out.close();
}

void WriteCellVectorDataVTK(string fileName,
                            const Mesh& mesh,
                            const vector<array<double, 3>>& data)
{
    ofstream out;
    out.open(fileName + ".vtk");

    out << "# vtk DataFile Version 2.0\n";
    out << "Vlasov-T\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << mesh.points.size() <<  " float\n";
    for (auto p : mesh.points)
        out << p->coords[0] << " " << p->coords[1] << " " << p->coords[2] << "\n";

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
        
    if (data.size() > 0)
    {
        out << "CELL_DATA " << mesh.tets.size() << "\n";
        out << "SCALARS " << "data" << " double 3\n";
        out << "LOOKUP_TABLE default\n";
        for (int i = 0; i < mesh.tets.size(); i++)
            out << data[i][0] << " " << data[i][1] << " " << data[i][2] << "\n";
    }
    
    out.close();
}

void WriteMeshVTK(string fileName,
                  const Mesh& mesh)
{
    // Write faces
    ofstream out;
    out.open(fileName + "_faces.vtk");

    out << "# vtk DataFile Version 2.0\n";
    out << "Vlasov-T\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << mesh.points.size() <<  " float\n";
    for (auto p : mesh.points)
        out << p->coords[0] << " " << p->coords[1] << " " << p->coords[2] << "\n";

    out << "CELLS " << mesh.faces.size() << " " << mesh.faces.size() * 4 << "\n";
    for (auto f : mesh.faces)
        out << 3 << " "
            << f->points[0]->index << " " 
            << f->points[1]->index << " "
            << f->points[2]->index << "\n";
    out << "\n";

    out << "CELL_TYPES " << mesh.faces.size() << "\n";
    for (auto f : mesh.faces)
        out << 5 << "\n";

    out.close();

    // Write tetrahedra
    out.open(fileName + "_tets.vtk");

    out << "# vtk DataFile Version 2.0\n";
    out << "Vlasov-T\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << mesh.points.size() <<  " float\n";
    for (auto p : mesh.points)
        out << p->coords[0] << " " << p->coords[1] << " " << p->coords[2] << "\n";

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
}

void WriteDistributionVTK(string fileName,
                          const VelocityGrid<Tensor3d>& velocityGrid,
                          const Tensor3d& distribution)
{
    ofstream out;
    out.open(fileName + ".vtk");

    out << "# vtk DataFile Version 2.0\n";
    out << "Vlasov-T\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    int n0 = velocityGrid.nCells[0];
    int n1 = velocityGrid.nCells[1];
    int n2 = velocityGrid.nCells[2];

    int nTotal = velocityGrid.nCellsTotal;

    out << "POINTS " << nTotal <<  " float\n";
    for (int i0 = 0; i0 < n0; i0++)
    {
        for (int i1 = 0; i1 < n1; i1++)
        {
            for (int i2 = 0; i2 < n2; i2++)
            {
                out << velocityGrid.At(i0, i1, i2)[0] << " " <<
                       velocityGrid.At(i0, i1, i2)[1] << " " << 
                       velocityGrid.At(i0, i1, i2)[2] << "\n";
            }
        }
    }
 
    out << "CELLS " << nTotal << " " << nTotal * 2 << "\n";
    for (int i0 = 0; i0 < n0; i0++)
    {
        for (int i1 = 0; i1 < n1; i1++)
        {
            for (int i2 = 0; i2 < n2; i2++)
            {
                out << 1 << " " << i2 + i1 * n1 + i0 * n1 * n0 << "\n";
            }
        }
    } 

    out << "CELL_TYPES " << nTotal << "\n";
    for (int i = 0; i < nTotal; i++)
        out << 1 << "\n";
    
    out << "CELL_DATA " << nTotal << "\n";
    out << "SCALARS " << "data" << " double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int i0 = 0; i0 < n0; i0++)
    {
        for (int i1 = 0; i1 < n1; i1++)
        {
            for (int i2 = 0; i2 < n2; i2++)
            {
                out << distribution(i0, i1, i2) << "\n";
            }
        }
    }
    out.close();
}
}