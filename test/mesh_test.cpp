#include <iostream>

#include "mesh.h"
#include "vtk.h"

using namespace std;

int main()
{
    Mesh mesh("../data/meshes/test_mesh_periodic.msh");
    VTK::WriteMesh("mesh", mesh);
}