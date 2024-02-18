#pragma once

#include "mesh.h"
#include "plasma_parameters.h"
#include "poisson.h"
#include "typedefs.h"
#include "log.h"

#include <vector>
#include <map>
#include <limits.h>

namespace VlasovTucker
{
// Boundary conditions
enum class FieldBCType { ConstantPotential, ChargedPlane };

struct FieldBC
{
    FieldBCType type;
    double potential;
    double chargeDenity;
};

enum class ParticleBCType { NonBoundary, Periodic, ConstantSource, AbsorbingWall, Free };

template <typename TensorType>
struct ParticleBC
{
    ParticleBCType type = ParticleBCType::NonBoundary;
    TensorType sourcePDF;
};

template <typename TensorType>
class Solver
{
public:
    Solver(const Mesh* mesh,
           const VelocityGrid* velocityGrid,
           PlasmaParameters<TensorType>* plasmaParameters);

    void SetFieldBC(int boundaryInd, const FieldBC& bc);
    void SetParticleBC(int boundaryInd, const ParticleBC<TensorType>& bc);

    void Solve(double timeStep, int nIterations);

private:
    void _PrecomputeNormalTensors();
    TensorType _PDFDerivative(const Tet* tet, int ind) const;

public:
    int writeStep = INT_MAX;

private:
    const Mesh* _mesh;
    const VelocityGrid* _vGrid;

    PoissonSolver _poissonSolver;

    PlasmaParameters<TensorType>* _plParams;

    std::vector<ParticleBC<TensorType>> _faceParticleBC;
    std::unordered_map<int, double> _wallCharge;
    std::unordered_map<int, double> _wallArea;

    Log _log;

    // Normal velocity tensors
    std::vector<std::array<TensorType, 4>> _vNormal;
    std::vector<std::array<TensorType, 4>> _vNormalAbs;
};
}