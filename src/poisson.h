#pragma once

#include <vector>
#include <array>

#include "mesh.h"
#include "constants.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace VlasovTucker
{

enum class SparseSolverType { SparseLU, ConjugateGradient };

class SparseSolver
{
public:
    SparseSolver(SparseSolverType type = SparseSolverType::SparseLU);

    void SetType(SparseSolverType type);

    void Compute(const Eigen::SparseMatrix<double>& system);

    Eigen::VectorXd Solve(Eigen::VectorXd rhs) const;

    SparseSolver& operator=(SparseSolver&& other);

public:
    Eigen::VectorXd guess;

private:
    SparseSolverType _type;

    Eigen::SparseLU<
        Eigen::SparseMatrix<double>,
        Eigen::COLAMDOrdering<int>
        > _solverLU;

    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<double>,
        Eigen::Upper
        > _solverCG;
};

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

    void SetBC(int boundaryInd, const PoissonBC& bc);

    void SetSparseSolverType(SparseSolverType type);

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

    // Initial guess for iterative methods
    void _SetGuess(const std::vector<double>& guess);
    
private:
    const Mesh* _mesh;
    std::vector<PoissonBC> _faceBC;

    bool _solutionIsUnique;

	Eigen::SparseMatrix<double> _system;

    SparseSolver _solver; 

    std::vector<double> _solution;
    std::vector<Vector3d> _gradient;
};

std::vector<double> EigenVectorToStdVector(const Eigen::VectorXd& eigenVector);

}