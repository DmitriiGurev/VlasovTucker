#include "poisson.h"

#include "log.h"
#include "smoother.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

namespace VlasovTucker
{
PoissonSolver::PoissonSolver() :
    _mesh(nullptr)
{}

PoissonSolver::PoissonSolver(const Mesh* mesh) :
    _mesh(mesh)
{
    _solutionIsUnique = false;

    // Label all faces as non-boundary
    _faceBC = vector<PoissonBC>(_mesh->faces.size());

    // Set periodic BC
    PoissonBC periodicBC;
    periodicBC.type = PoissonBCType::Periodic;
    for (auto planePair : _mesh->PeriodicBoundaries())
    {
        for (int boundaryMark : planePair)
            SetBC(boundaryMark, periodicBC);
    }
}

// Move assignment
PoissonSolver& PoissonSolver::operator=(PoissonSolver&& other)
{
    _mesh             = move(other._mesh);
    _faceBC           = move(other._faceBC);
    _solutionIsUnique = move(other._solutionIsUnique);
    _system           = move(other._system);
    _solution         = move(other._solution);
    _gradient         = move(other._gradient);

    // Cannot be moved
    if (_system.size() != 0)
        _solver.compute(_system);

    return *this;
}

void PoissonSolver::SetBC(int boundaryInd, const PoissonBC& bc)
{
    if (bc.type == PoissonBCType::Dirichlet)
        _solutionIsUnique = true;

    for (auto face : _mesh->EntityToFaces().at(boundaryInd))
        _faceBC[face->index] = bc;
}

void PoissonSolver::Initialize()
{
    int nEquations = _mesh->tets.size();

    // Assemble the system
    vector<Triplet> coeffs;
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        // Over-relaxed correction
        _FillLineCoeffs(coeffs, i);
    }

    _system = Eigen::SparseMatrix<double>(nEquations, nEquations);
	_system.setFromTriplets(coeffs.begin(), coeffs.end());
	_system.makeCompressed();

	_solver.analyzePattern(_system);
	_solver.factorize(_system);

    _solver.compute(_system);
}

void PoissonSolver::_FillLineCoeffs(vector<Triplet>& coeffs, int i) const
{
    if (!_solutionIsUnique && i == 0)
    {
        // If there are infinitely many solutions, choose the one with zero in the first
        // tetrahedron
        coeffs.push_back(Triplet(0, 0, 1));
        return;
    }

    for (int j = 0; j < 4; j++)
    {
        Tet* tet = _mesh->tets[i];
        Face* face = tet->faces[j];
        PoissonBCType bcType = _faceBC[face->index].type;

        if (bcType == PoissonBCType::NonBoundary)
        {
            Tet* adjTet = tet->adjTets[j];
            Point d = adjTet->centroid - tet->centroid;

            coeffs.push_back(Triplet(i, adjTet->index, face->area / (d.DotProduct(face->normal))));
            coeffs.push_back(Triplet(i, i, -face->area / (d.DotProduct(face->normal))));
            continue;
        }

        if (bcType == PoissonBCType::Periodic)
        {
            Tet* adjTet = tet->adjTets[j];
            Point d = (adjTet->centroid - tet->centroid);

            // Find the periodic counterpart of the common face
            int k = 0;
            while (adjTet->adjTets[k] != tet)
                k++;

            d = d + face->centroid - adjTet->faces[k]->centroid;

            coeffs.push_back(Triplet(i, adjTet->index, face->area / (d.DotProduct(face->normal))));
            coeffs.push_back(Triplet(i, i, -face->area / (d.DotProduct(face->normal))));
            continue;
        }

        if (bcType == PoissonBCType::Dirichlet)
        {
            Point d = face->centroid - tet->centroid;
            
            coeffs.push_back(Triplet(i, i, -face->area / (d.DotProduct(face->normal))));
            continue;
        }
    }
}

void PoissonSolver::Solve(vector<double> rho)
{
    // Make the total charge equal to zero
    if (!_solutionIsUnique)
        _MakeNeutral(rho);

    Eigen::VectorXd rhs(_mesh->tets.size());
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        rhs(i) = (-rho[i] / epsilon0) * _mesh->tets[i]->volume;

        // Over-relaxed correction
        _FillLineRHS(rhs, i);
    }

    if (_gradient.empty())
    {
        cout << Indent(2) << "Calculate the initial approximation of the gradient\n";

        // Solve without correction
        _solution = EigenVectorToStdVector(_solver.solve(rhs));
        _gradient = move(_Gradient());
    }

    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        // Add the cross-diffusion
        _CorrectRHS(rhs, i);
    }

    // Solve with correction
    _solution = EigenVectorToStdVector(_solver.solve(rhs));
    _gradient = move(_Gradient());
}

const vector<double>& PoissonSolver::Potential() const
{
    return _solution;
}

vector<Vector3d> PoissonSolver::ElectricField() const
{
    vector<Vector3d> field(_mesh->tets.size());
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        for (int k = 0; k < 3; k++)   
            field[i][k] = -_gradient[i][k];
    }
    return field;
}

void PoissonSolver::_MakeNeutral(std::vector<double>& rho) const
{
    double sum = 0;
    double vol = 0;

    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        sum += rho[i] * _mesh->tets[i]->volume;
        vol += _mesh->tets[i]->volume;
    }

    for (int i = 0; i < _mesh->tets.size(); i++)
        rho[i] -= sum / vol;
}

void PoissonSolver::_FillLineRHS(Eigen::VectorXd& rhs, int i) const
{
    if (!_solutionIsUnique && i == 0)
    {
        // If there are infinitely many solutions, choose the one with zero in the first
        // tetrahedron
        rhs(0) = 0;
        return;
    }

    Tet* tet = _mesh->tets[i];
    for (int j = 0; j < 4; j++)
    {
        Face* face = tet->faces[j];

        if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
        {
            Point d = face->centroid - tet->centroid;
            rhs(i) -= face->area / (d.DotProduct(face->normal)) * _faceBC[face->index].value;
            continue;
        }

        if (_faceBC[face->index].type == PoissonBCType::Neumann)
        {
            rhs(i) -= _faceBC[face->index].normalGrad * face->area;
            continue;
        }
    }
}

Vector3d PoissonSolver::_WeightedGradient(Tet* tet, int f) const
{
    Face* face = tet->faces[f];
    PoissonBCType bcType = _faceBC[face->index].type;

    Point d = face->centroid - tet->centroid;

    if (bcType == PoissonBCType::NonBoundary || bcType == PoissonBCType::Periodic)
    {
        Tet* adjTet = tet->adjTets[f];
        Point adjD = face->centroid - adjTet->centroid;

        double g = adjD.Abs() / (adjD.Abs() + d.Abs()); 
        double adjG = d.Abs() / (adjD.Abs() + d.Abs());

        Vector3d grad = _gradient[tet->index];
        Vector3d adjGrad = _gradient[adjTet->index];

        Vector3d weighted;
        for (int i = 0; i < 3; i++)
            weighted[i] = grad[i] * g + adjGrad[i] * adjG;
        return weighted;
    }

    if (bcType == PoissonBCType::Dirichlet)
    {
        return _gradient[tet->index];
    }
}

void PoissonSolver::_CorrectRHS(Eigen::VectorXd& rhs, int i) const
{
    if (!_solutionIsUnique && i == 0)
        return;

    for (int j = 0; j < 4; j++)
    {
        Tet* tet = _mesh->tets[i];

        Face* face = tet->faces[j];
        PoissonBCType bcType = _faceBC[face->index].type;

        Point d = face->centroid - tet->centroid;

        double crossDiffusion = 0;
        if (bcType == PoissonBCType::NonBoundary)
        {
            Tet* adjTet = tet->adjTets[j];
            Point e = adjTet->centroid - tet->centroid;
            e = e / e.Abs();

            Point weightedGrad = _WeightedGradient(tet, j);
            crossDiffusion = face->area * weightedGrad.DotProduct(face->normal -
                e / e.DotProduct(face->normal));
        }
        else if (bcType == PoissonBCType::Periodic)
        {
            Tet* adjTet = tet->adjTets[j];
            Point e = adjTet->centroid - tet->centroid;

            int k = 0;
            while (adjTet->adjTets[k] != tet)
                k++;

            e = e + face->centroid - adjTet->faces[k]->centroid;
            e = e / e.Abs();

            Point weightedGrad = _WeightedGradient(tet, j);
            crossDiffusion = face->area * weightedGrad.DotProduct(face->normal -
                e / e.DotProduct(face->normal));
        }
        else if (bcType == PoissonBCType::Dirichlet)
        {
            Point e = face->centroid - tet->centroid;
            e = e / e.Abs();

            Point weightedGrad = _WeightedGradient(tet, j);
            crossDiffusion = face->area * weightedGrad.DotProduct(face->normal -
                e / e.DotProduct(face->normal));
        }

        rhs(i) -= crossDiffusion;
    }
}

Vector3d PoissonSolver::_TetLSG(Tet* tet) const
{
    double val = _solution[tet->index];

    array<double, 4> adjVal;
    array<Point, 4> dist;
    
    for (int i = 0; i < 4; i++)
    {
        Tet* adjTet = tet->adjTets[i];

        Face* face = tet->faces[i];
        PoissonBCType bcType = _faceBC[face->index].type;

        if (bcType == PoissonBCType::NonBoundary)
        {
            adjVal[i] = _solution[tet->adjTets[i]->index];
            dist[i] = adjTet->centroid - tet->centroid;
            continue;
        }
        
        if (bcType == PoissonBCType::Dirichlet)
        {
            adjVal[i] = _faceBC[face->index].value;
            dist[i] = face->centroid - tet->centroid;
            continue;
        }

        if (bcType == PoissonBCType::Neumann)
        {
            Point d = face->centroid - tet->centroid;
            Point x = face->normal * face->normal.DotProduct(d);
            
            dist[i] = x;
            adjVal[i] = val + x.Abs() * _faceBC[face->index].normalGrad;

            continue;
        }

        if (bcType == PoissonBCType::Periodic)
        {
            adjVal[i] = _solution[tet->adjTets[i]->index];
            Point d = (adjTet->centroid - tet->centroid);

            int k = 0;
            while (adjTet->adjTets[k] != tet)
                k++;

            dist[i] = d + face->centroid - adjTet->faces[k]->centroid;
            continue;
        }
    }

    array<double, 4> weights;
    for (int i = 0; i < 4; i++)
        weights[i] = 1 / dist[i].Abs();

    Eigen::Matrix3d m;
    for (int k = 0; k < 3; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            m(k, i) = 0;
            for (int j = 0; j < 4; j++)
                m(k, i) += 2 * weights[j] * dist[j][k] * dist[j][i];
        }
    }

    if (m.determinant() == 0)
        throw runtime_error("Degenerate matrix in LSG.");

    Eigen::Vector3d rhs;
    for (int k = 0; k < 3; k++)
    {
        rhs(k) = 0;
        for (int j = 0; j < 4; j++)
            rhs(k) -= 2 * weights[j] * dist[j][k] * (val - adjVal[j]);
    }

    Eigen::Vector3d grad = m.fullPivLu().solve(rhs);
    return {grad(0), grad(1), grad(2)};
}

vector<Vector3d> PoissonSolver::_Gradient() const
{
    assert(!_solution.empty());

    vector<Vector3d> gradient(_mesh->tets.size());
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        Tet* tet = _mesh->tets[i];
        gradient[i] = _TetLSG(tet);
    }
    return gradient;
}

vector<double> EigenVectorToStdVector(const Eigen::VectorXd& eigenVector)
{
    int nElements = eigenVector.size();
    vector<double> vector(nElements);
    for (int i = 0; i < nElements; i++)
        vector[i] = eigenVector(i);
    return vector;
}
}