#pragma once

#include "mesh.h"
#include "plasma_parameters.h"
#include "typedefs.h"
#include "log.h"

#include <vector>
#include <map>
#include <limits.h>

// TODO: Add BC

namespace VlasovTucker
{
template <typename TensorType>
class Solver
{
public:
    Solver(const Mesh* mesh,
           const VelocityGrid* velocityGrid,
           PlasmaParameters<TensorType>* plasmaParameters);

    void Solve(int nIterations);

private:
    void _PrecomputeNormalTensors();
    TensorType _PDFDerivative(const Tet* tet, int ind) const;

public:
    int writeStep = INT_MAX;

private:
    const Mesh* _mesh;
    const VelocityGrid* _vGrid;

    PlasmaParameters<TensorType>* _plParams;

    // std::vector<ParticleBC> _boundaryConditions;

    Log _log;

    // Normal velocity tensors
    std::vector<std::array<TensorType, 4>> _vNormal;
    std::vector<std::array<TensorType, 4>> _vNormalAbs;
};
}