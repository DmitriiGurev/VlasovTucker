// Gmsh project created on Wed Nov 15 21:17:55 2023
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 0.05, 0.05};
//+
Physical Volume("Volume", 13) = {1};
//+
Physical Surface("Boundary: Periodic A", 14) = {6, 5};
//+
Physical Surface("Boundary: Periodic B", 15) = {4, 3};
//+
Physical Surface("Boundary: Periodic C", 16) = {2, 1};
//+
Dilate {{0, 0, 0}, {1, 0.1, 0.1}} {
  Volume{1}; 
}
//+
Physical Surface(17) = {6, 5, 2, 3, 4, 1};
//+
Physical Surface("Poisson: Periodic", 18) = {5, 1, 4, 6, 3, 2};
//+
Dilate {{0, 0, 0}, {1, 2, 2}} {
  Volume{1}; 
}
//+
Dilate {{0, 0, 0}, {1, 2, 2}} {
  Volume{1}; 
}
//+
Physical Surface(" Boundary", 14) -= {6};
//+
Physical Surface(" Boundary", 14) -= {4};
//+
Physical Surface(" Boundary", 14) -= {5};
//+
Physical Surface(" Boundary", 15) -= {3};
//+
Physical Surface(" Boundary", 15) -= {6};
//+
Physical Surface(" Boundary", 15) -= {6};
//+
Physical Surface(" Boundary", 15) -= {4};
//+
Physical Surface(" Boundary", 16) -= {3};
//+
Physical Surface(" Boundary", 16) -= {1};
//+
Physical Surface(" Boundary", 16) -= {2};
//+
Physical Surface(17) -= {6};
//+
Physical Surface(17) -= {5};
//+
Physical Surface(17) -= {2};
//+
Physical Surface(17) -= {1};
//+
Physical Surface(17) -= {4, 3};
//+
Physical Surface(" Poisson", 18) -= {6};
//+
Physical Surface(" Poisson", 18) -= {1, 4};
//+
Physical Surface(" Poisson", 18) -= {2};
//+
Physical Surface(" Poisson", 18) -= {3, 4};
//+
Physical Surface(" Poisson", 18) -= {5, 6};
//+
Physical Surface("Boundary: Periodic A", 19) = {3, 4};
//+
Physical Surface("Boundary: Periodic B", 20) = {5, 6};
//+
Physical Surface("Poisson: Periodic", 21) = {5, 4, 6, 3};
//+
Physical Surface("Poisson: Neumann", 22) = {2};
//+
Physical Surface("Poisson: Dirichlet", 23) = {1};
//+
Physical Surface("Particles: Absorbing", 24) = {2};
//+
Physical Surface("Particles: Dirichlet", 25) = {1};
//+
Physical Surface("Particles: Periodic", 26) = {6, 5, 4, 3};
