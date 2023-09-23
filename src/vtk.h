#pragma once

#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <string>

#include "mesh.h"
#include "typedefs.h"
#include "velocity_grid.h"

namespace VTK
{
void WriteCellScalarData(std::string fileName,
                         const Mesh& mesh,
                         const std::vector<double>& data = {});

void WriteCellVectorData(std::string fileName,
                         const Mesh& mesh,
                         const std::vector<std::array<double, 3>>& data = {});

void WriteMesh(std::string fileName,
               const Mesh& mesh);

void WriteDistribution(std::string fileName,
                       const VelocityGrid<Tensor>& velocityGrid,
                       const Tensor& distribution);
}