#include "velocity_grid.h"

#include <iostream>

using namespace std;

template <>
VelocityGrid<Tensor>::VelocityGrid(array<int, 3> nCells,
                                   array<double, 3> minV,
                                   array<double, 3> maxV) :
    nCells(nCells), minV(minV), maxV(maxV)
{
    nCellsTotal = nCells[0] * nCells[1] * nCells[2];

    for (int j = 0; j < 3; j++)
        step[j] = (maxV[j] - minV[j]) / (nCells[j] - 1);

    cellVolume = step[0] * step[1] * step[2];

    for (int j = 0; j < 3; j++)
    {
        v[j] = Tensor(nCells[0], nCells[1], nCells[2]);

        std::array<int, 3> ind;
        for (ind[0] = 0; ind[0] < nCells[0]; ind[0]++)
        {
            for (ind[1] = 0; ind[1] < nCells[1]; ind[1]++)
            {
                for (ind[2] = 0; ind[2] < nCells[2]; ind[2]++)
                {
                    v[j](ind[0], ind[1], ind[2]) = minV[j] + ind[j] * step[j];
                } 
            } 
        } 
    }

    // Compute the central difference matrices
    for (int j = 0; j < 3; j++)
    {
        d[j] = Eigen::MatrixXd::Zero(nCells[j], nCells[j]);
        d[j](0, 1) = 1;
        d[j](0, nCells[j] - 1) = -1;
        for (int i = 1; i < nCells[j] - 1; i++)
        {
            d[j](i, i + 1) = 1;
            d[j](i, i - 1) = -1;
        }
        d[j](nCells[j] - 1, 0) = 1;
        d[j](nCells[j] - 1, nCells[j] - 2) = -1;

        d[j] /= 2 * step[j];
    }
}

template <>
array<double, 3> VelocityGrid<Tensor>::At(int i0, int i1, int i2) const
{
    return {minV[0] + i0 * step[0],
            minV[1] + i1 * step[1],
            minV[2] + i2 * step[2]};
}