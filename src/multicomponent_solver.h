#pragma once

#include "solver.h"

namespace VlasovTucker
{
template <typename TensorType>
class MulticomponentSolver
{
public:
    MulticomponentSolver(Solver<TensorType>* baseSolver);

    void AddSolver(Solver<TensorType>* solver);

    void Solve();

public:
    double timeStep = 0;

    // Update the particular PDF once in several steps
    std::map<Solver<TensorType>*, int> stepMultipliers;

    int nIterations = 0;
    int writeStep = INT_MAX;

private:
    std::vector<Solver<TensorType>*> _solvers;

    Log _log;
};
}