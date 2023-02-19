#pragma once

#include <vector>

#include "mesh.h"

#include <Eigen/Sparse>

using namespace std;

// TODO: Take it from constants.h
const int eps0 = 1.0;

class PoissonSolver
{
public:
    PoissonSolver(const Mesh& mesh);

    std::vector<double> Solve(std::vector<double> rho) const;

private:
    Mesh _mesh;

    enum class BCType
    {
        NonBoundary,
        Neumann,
        Dirichlet,
        // Periodic
    };

    std::vector<BCType> _faceTypes = std::vector<BCType>(_mesh.faces.size());

    Eigen::VectorXd				_rhs;
	Eigen::SparseMatrix<double> _system;

	Eigen::SparseLU<Eigen::SparseMatrix<double>,
                    Eigen::COLAMDOrdering<int>> _solver;
};