#pragma once

#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <string>

#include "mesh.h"
#include "typedefs.h"
#include "velocity_grid.h"

namespace VlasovTucker
{
void WriteCellScalarDataVTK(std::string fileName,
                            const Mesh& mesh,
                            const std::vector<double>& data = {});

void WriteCellVectorDataVTK(std::string fileName,
                            const Mesh& mesh,
                            const std::vector<std::array<double, 3>>& data = {});

void WriteMeshVTK(std::string fileName,
                  const Mesh& mesh);

void WriteDistributionVTK(std::string fileName,
                          const VelocityGrid<Tensor3d>& velocityGrid,
                          const Tensor3d& distribution);
}