#pragma once

#include "mesh.h"
#include "plasma_parameters.h"
#include "typedefs.h"

#include <vector>

class Solver
{
public:
    Solver(const Mesh* mesh,
           const VelocityGrid<Tensor>* velocityGrid,
           PlasmaParameters* plasmaParameters);

    void Solve(int nIterations);

private:
    const Mesh* _mesh;
    const VelocityGrid<Tensor>* _vGrid;

    PlasmaParameters* _plasmaParams;

    std::vector<std::array<Tensor, 4>> _vNormal;
    std::vector<std::array<Tensor, 4>> _vNormalAbs;
};