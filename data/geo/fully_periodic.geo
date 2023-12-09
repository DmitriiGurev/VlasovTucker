// Gmsh project created on Tue Feb 21 21:52:33 2023
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Volume("Volume", 13) = {1};
//+
Physical Surface("Poisson: Periodic", 15) = {5, 3, 1, 4, 6, 2};
//+
Physical Surface("Boundary: Periodic A", 16) = {5, 6};
//+
Physical Surface("Boundary: Periodic B", 17) = {2, 1};
//+
Physical Surface("Boundary: Periodic C", 18) = {3, 4};
