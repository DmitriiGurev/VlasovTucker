// Gmsh project created on Fri Feb 23 12:58:32 2024
SetFactory("OpenCASCADE");
//+
Torus(1) = {0, 0, 0, 1, 0.7, 2*Pi};
//+
Physical Volume("Volume", 3) = {1};
//+
Physical Surface("Poisson: Neumann", 4) = {1};
