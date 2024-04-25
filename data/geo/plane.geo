// Gmsh project created on Sat Mar 30 18:10:21 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 0.1};
//+
Physical Volume("Vol", 13) = {1};
//+
Physical Surface("Surf", 14) = {6, 5, 2, 3, 1, 4};
