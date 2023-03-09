#pragma once

#include <fstream>
#include <vector>
#include <map>
#include <string>

#include "mesh.h"
#include "typedefs.h"
#include "velocity_grid.h"

namespace VTK
{
void WriteCellData(std::string fileName,
                   const Mesh& mesh,
                   const std::vector<double>& data = {});

void WriteMesh(std::string fileName,
               const Mesh& mesh);

void WriteDistribution(std::string fileName,
                       const VelocityGrid<Tensor>& velocityGrid,
                       const Tensor& distribution);
}