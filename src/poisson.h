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

// (line, row, value)
typedef Eigen::Triplet<double> Triplet;

class PoissonSolver
{
public:
    PoissonSolver();
    PoissonSolver(const Mesh* mesh);

    PoissonSolver& operator=(PoissonSolver&& other);

    void SetBC(int boundaryInd, const PoissonBC& bc);

    void Initialize();

    void Solve(std::vector<double> rho);

    const std::vector<double>& Potential() const;
    std::vector<Vector3d> ElectricField() const;

private:
    std::vector<Vector3d> _Gradient();

    void _FillLine(vector<Triplet>& coeffs, int i);

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
    std::vector<Vector3d> _gradient;
};
}