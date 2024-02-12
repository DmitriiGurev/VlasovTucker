#pragma once

#include <vector>
#include <array>

#include "mesh.h"
#include "constants.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;

namespace VlasovTucker
{
enum class PoissonBCType { NonBoundary, Neumann, Dirichlet, Periodic };

struct PoissonBC
{
    PoissonBCType type = PoissonBCType::NonBoundary;
    double value = 0;
    double normalGrad = 0;
};

class PoissonSolver
{
public:
    PoissonSolver(const Mesh* mesh);

    void SetBC(int boundaryInd, const PoissonBC& bc);

    void Initialize();

    void Solve(std::vector<double> rho);

    const std::vector<double>& Potential() const;
    // TODO: Change to Point
    std::vector<std::array<double, 3>> ElectricField() const;

private:
    std::vector<Point> Gradient();

private:
    const Mesh* _mesh;
    std::vector<PoissonBC> _faceBC;

    bool _solutionIsUnique;

	Eigen::SparseMatrix<double> _system;

	Eigen::SparseLU<
        Eigen::SparseMatrix<double>,
        Eigen::COLAMDOrdering<int>
        > _solver;

    std::vector<double> _solution;
    std::vector<Point> _gradient;
};
}