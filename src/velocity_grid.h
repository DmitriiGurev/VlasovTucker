#pragma once

#include <array>

#include "tensor_type.h"

template <class TensorType>
struct VelocityGrid
{
    VelocityGrid(std::array<int, 3> nCells,
                 std::array<double, 3> maxV,
                 std::array<double, 3> minV);

    std::array<double, 3> At(int i0, int i1, int i2) const;

    std::array<int, 3> nCells;
    int nCellsTotal;

    std::array<double, 3> step;
    double cellVolume;

    std::array<double, 3> maxV;
    std::array<double, 3> minV;

    std::array<TensorType, 3> v;
};