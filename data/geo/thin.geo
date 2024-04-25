// Gmsh project created on Sat Mar 30 17:59:27 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 0.01, 0.01};
//+
Physical Volume("Vol", 13) = {1};
//+
Physical Surface("Surf", 14) = {6, 4, 5, 3, 2, 1};
