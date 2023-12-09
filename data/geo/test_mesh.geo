// Gmsh project created on Thu Feb 16 18:31:53 2023
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Volume("Volume", 13) = {1};
//+
Physical Surface("Boundary: Periodic A", 14) = {1, 2};
//+
Physical Surface("Boundary: Periodic B", 16) = {4, 3};
//+
Physical Surface("Boundary: Dirichlet", 17) = {6, 5};
//+
Physical Surface("Poisson: Dirichlet", 18) = {6};
//+
Physical Surface("Poisson: Neumann", 19) = {5};
//+
Physical Surface("Poisson: Periodic", 20) = {1, 4, 2, 3};//+
Physical Surface("Boundary: Periodic C", 18) += {3, 2};
