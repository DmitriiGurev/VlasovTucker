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
Dilate {{0, 0, 0}, {1, 0.5, 0.5}} {
  Volume{1}; 
}
//+
Dilate {{0, 0, 0}, {1, 0.5, 0.5}} {
  Volume{1}; 
}
