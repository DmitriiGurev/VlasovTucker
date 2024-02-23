#include <iostream>

#include "mesh.h"
#include "vtk.h"

using namespace VlasovTucker;
using namespace std;

int main(int argc, char *argv[])
{
    Mesh mesh(argv[1]);
    mesh.Reconstruct();
    WriteMeshVTK("mesh", mesh);
}