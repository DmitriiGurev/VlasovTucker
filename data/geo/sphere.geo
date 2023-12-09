// Gmsh project created on Tue Feb 21 21:34:36 2023
SetFactory("OpenCASCADE");
//+
Sphere(1) = { ,  ,  , 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Sphere(1) = {-0, -0, -0, 0.5, -Pi/2, Pi/2, 2*Pi};
//+
Physical Surface("Boundary: Dirichlet", 4) = {1};
//+
Physical Surface("Poisson: Dirichlet", 5) = {1};
//+
Physical Volume("Volume", 6) = {1};
