#pragma once

#include <fstream>
#include <vector>
#include <map>
#include <string>

#include <mesh.h>

void WriteToVTK(std::string fileName,
                const Mesh& mesh,
                const std::map<std::string, std::vector<double>>& data);