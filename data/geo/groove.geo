// Gmsh project created on Sun Mar 03 11:14:19 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 0.95, 0.1, 0.05};
//+
Box(2) = {0.95, 0.025, 0, 0.05, 0.05, 0.05};
//+
BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; }
//+
Physical Volume("Volume", 25) = {1};
//+
Physical Surface("Poisson: Dirichlet", 26) = {1};
//+
Physical Surface("Poisson: Neumann", 27) = {6, 8, 10, 9, 7};
//+
Physical Surface("Periodic 1", 28) = {3, 5};
//+
Physical Surface("Periodic 2", 29) = {2, 4};
