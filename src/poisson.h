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
    Point gradient = Point({0, 0, 0});
};

class PoissonSolver
{
public:
    PoissonSolver(const Mesh* mesh);

    void Solve(std::vector<double> rho);

    const std::vector<double>& Potential() const;
    // TODO: Change to Point
    std::vector<std::array<double, 3>> ElectricField() const;

    // void SetBC(plane, function)

private:
    std::vector<Point> Gradient();

private:
    const Mesh* _mesh;

    bool _solutionIsUnique;

    // map<plane, function> _faceBCFunctions
    std::vector<PoissonBC> _faceBC = std::vector<PoissonBC>(_mesh->faces.size());

	Eigen::SparseMatrix<double> _system;

	Eigen::SparseLU<
        Eigen::SparseMatrix<double>,
        Eigen::COLAMDOrdering<int>
        > _solver;

    std::vector<double> _solution /*= std::vector<double>(_mesh->tets.size())*/;
    std::vector<Point> _gradient /*= std::vector<Point>(_mesh->tets.size())*/;
};
}