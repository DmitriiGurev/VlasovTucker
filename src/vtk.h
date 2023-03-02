#pragma once

#include <fstream>
#include <vector>
#include <map>
#include <string>

#include "mesh.h"

enum class ModeVTK
{
    CellData,
    Mesh
};

void WriteToVTK(ModeVTK mode,
                std::string fileName,
                const Mesh& mesh,
                const std::map<std::string, std::vector<double>>& data = {});