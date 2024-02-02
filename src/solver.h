#pragma once

#include "mesh.h"
#include "plasma_parameters.h"
#include "typedefs.h"
#include "log.h"
#include "tucker.h"

#include <vector>
#include <map>
#include <limits.h>

// TODO: Add BC

class Solver
{
public:
    Solver(const Mesh* mesh,
           const VelocityGrid<Tensor>* velocityGrid,
           PlasmaParameters* plasmaParameters);

    void Solve(int nIterations);

private:
    void _PrecomputeNormalTensors();
    Tensor _PDFDerivative(const Tet* tet, int ind) const;

    void _PrecomputeGradCoeffs();
    std::array<Tensor, 3> _Gradient(Tet* tet) const;

public:
    int writeStep = INT_MAX;

private:
    const Mesh* _mesh;
    const VelocityGrid<Tensor>* _vGrid;

    PlasmaParameters* _plParams;

    // std::vector<ParticleBC> _boundaryConditions;

    Log _log;

    // Normal velocity tensors
    std::vector<std::array<Tensor, 4>> _vNormal;
    std::vector<std::array<Tensor, 4>> _vNormalAbs;

    // Compressed normal velocity tensors
    std::vector<std::array<Tucker, 4>> _comprVNormal;
    std::vector<std::array<Tucker, 4>> _comprVNormalAbs;

    // Coefficients for LSM gradient
    std::vector<Eigen::Matrix3d> _gradMatrices;
    // std::vector<Eigen::Vector3d> _gradRHSCoeffs;
    
    // ...
    // TODO: Change to unordered map
    std::map<std::tuple<Tet*, Tet*>, Point> _distances;
};