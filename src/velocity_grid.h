#pragma once

#include "typedefs.h"

#include <Eigen/Dense>

#include <array>

namespace VlasovTucker
{
struct VelocityGrid
{
    VelocityGrid(std::array<int, 3> nCells,
                 std::array<double, 3> minV,
                 std::array<double, 3> maxV);

    std::array<double, 3> At(int i0, int i1, int i2) const;

    std::array<int, 3> nCells;
    int nCellsTotal;

    std::array<double, 3> step;
    double cellVolume;

    std::array<double, 3> maxV;
    std::array<double, 3> minV;

    std::array<Tensor3d, 3> v;

    // Central difference matrices
    std::array<Eigen::MatrixXd, 3> d;
};
}