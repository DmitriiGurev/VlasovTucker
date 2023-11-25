#pragma once

#include <vector>
#include <array>

#include "mesh.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;

// TODO: Take it from constants.h
const int eps0 = 1.0;

class PoissonSolver
{
public:
    PoissonSolver(const Mesh* mesh);

    std::vector<double> Solve(std::vector<double> rho) const;

    std::vector<array<double, 3>> ElectricField(const std::vector<double>& potential) const;

private:
    const Mesh* _mesh;

    enum class BCType
    {
        NonBoundary,
        Neumann,
        Dirichlet,
        Periodic
    };

    bool _solutionIsUnique;

    std::vector<BCType> _faceTypes = std::vector<BCType>(_mesh->faces.size());

    Eigen::VectorXd				_rhs;
	Eigen::SparseMatrix<double> _system;

	Eigen::SparseLU<Eigen::SparseMatrix<double>,
                    Eigen::COLAMDOrdering<int>> _solver;

    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
    //                          Eigen::Lower|Eigen::Upper> _solver;
};

std::array<double, 3> Gradient(Point p, double val,
                               const std::vector<Point>& neighPoints,
                               const std::vector<double>& neighVals);