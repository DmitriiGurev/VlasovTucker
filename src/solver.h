#pragma once

#include "mesh.h"
#include "plasma_parameters.h"
#include "typedefs.h"
#include "log.h"

#include <vector>
#include <map>
#include <limits.h>

// TODO: Add BC

// Flux reconstruction scheme
enum class FluxScheme {
    Upwind,      // 1st order
    UpwindFROMM, // 2nd order
    UpwindQUICK, // 2nd order (?)
};

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

    FluxScheme fluxScheme = FluxScheme::Upwind;

private:
    const Mesh* _mesh;
    const VelocityGrid<Tensor>* _vGrid;

    PlasmaParameters* _plParams;

    Log _log;

    // Normal velocity tensors
    std::vector<std::array<Tensor, 4>> _vNormal;
    std::vector<std::array<Tensor, 4>> _vNormalAbs;

    // Coefficients for LSM gradient
    std::vector<Eigen::Matrix3d> _gradMatrices;
    // std::vector<Eigen::Vector3d> _gradRHSCoeffs;
    
    // ...
    // TODO: Change to unordered map
    std::map<std::tuple<Tet*, Tet*>, Point> _distances;
};