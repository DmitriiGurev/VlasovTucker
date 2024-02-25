#pragma once

#include <vector>
#include <array>

#include "mesh.h"
#include "constants.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

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
    // Least-squares gradient
    Vector3d _TetLSG(Tet* tet) const;
    std::vector<Vector3d> _Gradient() const;
    Vector3d _WeightedGradient(Tet* tet, int f) const;

    void _FillLineCoeffs(std::vector<Triplet>& coeffs, int i) const;
    void _FillLineRHS(Eigen::VectorXd& rhs, int i) const;
    void _CorrectRHS(Eigen::VectorXd& rhs, int i) const;

    void _MakeNeutral(std::vector<double>& rho) const;

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

std::vector<double> EigenVectorToStdVector(const Eigen::VectorXd& eigenVector);

}