// Gmsh project created on Tue Feb 21 21:34:36 2023
SetFactory("OpenCASCADE");
//+
Sphere(1) = {-0, -0, -0, 1.0, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("Boundary: Dirichlet", 4) = {1};
//+
Physical Volume("Volume", 5) = {1};
