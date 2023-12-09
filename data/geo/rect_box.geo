// Gmsh project created on Wed Dec 06 11:55:07 2023
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 0.2};
//+
Physical Volume("Volume", 13) = {1};
//+
Physical Surface("Poisson: Dirichlet", 14) = {3, 4, 1, 6, 5, 2};
//+
Physical Surface(" Poisson", 14) -= {6, 5};
//+
Physical Surface("Poisson: Periodic", 15) = {6, 5};
//+
Physical Surface(" Poisson", 14) -= {2, 3, 1, 4};
//+
Physical Surface("Poisson: Periodic", 15) += {2, 3, 1, 4};
//+
Physical Surface("Boundary: Periodic A", 16) = {6, 5};
//+
Physical Surface("Boundary: Periodic B", 17) = {3, 4};
//+
Physical Surface("Boundary: Periodic C", 18) = {1, 2};
