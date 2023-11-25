#pragma once

#include "mesh.h"
#include "plasma_parameters.h"
#include "typedefs.h"
#include "log.h"

#include <vector>
#include <limits.h>

class Solver
{
public:
    Solver(const Mesh* mesh,
           const VelocityGrid<Tensor>* velocityGrid,
           PlasmaParameters* plasmaParameters);

    void Solve(int nIterations);

private:
    void _ComputeNormalTensors();
    Tensor _PDFDerivative(const Tet* tet, int ind) const;

public:
    int writeStep = INT_MAX;

private:
    const Mesh* _mesh;
    const VelocityGrid<Tensor>* _vGrid;

    Log _log;

    PlasmaParameters* _plParams;

    Log _log;

    // Normal velocity tensors
    std::vector<std::array<Tensor, 4>> _vNormal;
    std::vector<std::array<Tensor, 4>> _vNormalAbs;
};