// Gmsh project created on Mon Feb 26 15:11:59 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Volume("Volume", 13) = {1};
//+
Physical Surface("Periodic", 14) = {6, 5};
//+
Physical Surface("Dirichlet", 15) = {2, 1};
//+
Physical Surface("Neumann", 16) = {3, 4};
