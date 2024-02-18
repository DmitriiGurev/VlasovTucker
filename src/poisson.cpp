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
    _faceBC = std::vector<PoissonBC>(_mesh->faces.size());

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

    // TODO: Use a map from entities to faces
    for (int i = 0; i < _mesh->faces.size(); i++)
    {
        Face* face = _mesh->faces[i];
        if (face->entity == boundaryInd)
            _faceBC[i] = bc;
    }
}

void PoissonSolver::Initialize()
{
    int nEquations = _mesh->tets.size();

    // Assemble the system using Triplets (line, row, value)
    typedef Eigen::Triplet<double> Triplet;
    vector<Triplet> coeffs;

    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        if (!_solutionIsUnique && i == 0)
        {
            coeffs.push_back(Triplet(0, 0, 1));
            continue;
        }

        // Over-relaxed correction
        for (int j = 0; j < 4; j++)
        {
            Tet* tet = _mesh->tets[i];
            Face* face = tet->faces[j];

            if (_faceBC[face->index].type == PoissonBCType::NonBoundary)
            {
                Tet* adjTet = tet->adjTets[j];
                Point d = adjTet->centroid - tet->centroid;

                coeffs.push_back(Triplet(i, adjTet->index, 
                    face->area / (d.DotProduct(face->normal))));

                coeffs.push_back(Triplet(i, i,
                    -face->area / (d.DotProduct(face->normal))));
            }
            else if (_faceBC[face->index].type == PoissonBCType::Periodic)
            {
                Tet* adjTet = tet->adjTets[j];
                Point d = (adjTet->centroid - tet->centroid);

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                d = d + face->centroid - adjTet->faces[k]->centroid;

                coeffs.push_back(Triplet(i, adjTet->index, 
                    face->area / (d.DotProduct(face->normal))));
                    
                coeffs.push_back(Triplet(i, i,
                    -face->area / (d.DotProduct(face->normal))));
            }
            else if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
            {
                Point d = face->centroid - tet->centroid;
                
                coeffs.push_back(Triplet(i, i,
                    -face->area / (d.DotProduct(face->normal))));
            }
        }
    }

    _system = Eigen::SparseMatrix<double>(nEquations, nEquations);
	_system.setFromTriplets(coeffs.begin(), coeffs.end());
	_system.makeCompressed();

	// _solver.analyzePattern(_system);
	// _solver.factorize(_system);

    _solver.compute(_system);
}

void PoissonSolver::Solve(vector<double> rho)
{
    // Make the total charge equal to zero
    if (!_solutionIsUnique)
    {
        cout << Indent(2) << "Solution is not unique\n";
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

    Eigen::VectorXd rhs(_mesh->tets.size());
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        if (!_solutionIsUnique && i == 0)
        {
            rhs(0) = 0;
            continue;       
        }

        Tet* tet = _mesh->tets[i];

        rhs(i) = (-rho[i] / epsilon0) * tet->volume;

        // Over-relaxed correction
        for (int j = 0; j < 4; j++)
        {
            Face* face = tet->faces[j];

            if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
            {
                Point d = face->centroid - tet->centroid;

                rhs(i) -= face->area / (d.DotProduct(face->normal)) * _faceBC[face->index].value;
            }
            else if (_faceBC[face->index].type == PoissonBCType::Neumann)
            {
                rhs(i) -= _faceBC[face->index].normalGrad * face->area;
            }
        }
    }

    if (_gradient.empty())
    {
        cout << Indent(2) << "Calculate the initial approximation of the gradient\n";

        Eigen::VectorXd solution = _solver.solve(rhs);

        _solution.resize(_mesh->tets.size());
        for (int i = 0; i < _mesh->tets.size(); i++)
            _solution[i] = solution(i);
        
        _gradient = move(Gradient());
    }

    // Cross-diffusion
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        if (!_solutionIsUnique && i == 0)
            continue;       

        for (int j = 0; j < 4; j++)
        {
            Tet* tet = _mesh->tets[i];
            Face* face = tet->faces[j];
            Point d = face->centroid - tet->centroid;

            double crossDiff;
            
            if (_faceBC[face->index].type == PoissonBCType::NonBoundary)
            {
                Tet* adjTet = tet->adjTets[j];
                Point adjD = face->centroid - adjTet->centroid;

                Point e = adjTet->centroid - tet->centroid;
                e = e / e.Abs();

                double g = adjD.Abs() / (adjD.Abs() + d.Abs()); 
                double adjG = d.Abs() / (adjD.Abs() + d.Abs());
                Point grad = _gradient[i];
                Point adjGrad = _gradient[adjTet->index];

                crossDiff = face->area * (grad * g + adjGrad * adjG).
                    DotProduct(face->normal - e / e.DotProduct(face->normal));
            }
            else if (_faceBC[face->index].type == PoissonBCType::Periodic)
            {
                Tet* adjTet = tet->adjTets[j];
                Point adjD = face->centroid - adjTet->centroid;

                Point e = adjTet->centroid - tet->centroid;

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                e = e + face->centroid - adjTet->faces[k]->centroid;
                e = e / e.Abs();

                double g = adjD.Abs() / (adjD.Abs() + d.Abs()); 
                double adjG = d.Abs() / (adjD.Abs() + d.Abs());
                Point grad = _gradient[i];
                Point adjGrad = _gradient[adjTet->index];

                crossDiff = face->area * (grad * g + adjGrad * adjG).
                    DotProduct(face->normal - e / e.DotProduct(face->normal));
            }
            else if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
            {
                Point e = face->centroid - tet->centroid;
                e = e / e.Abs();

                crossDiff = face->area * _gradient[i].DotProduct(face->normal -
                    e / e.DotProduct(face->normal));
            }
            else if (_faceBC[face->index].type == PoissonBCType::Neumann)
            {
                crossDiff = 0;
            }

            rhs(i) -= crossDiff;
        }
    }

    Eigen::VectorXd solution = _solver.solve(rhs);

    for (int i = 0; i < _mesh->tets.size(); i++)
        _solution[i] = solution(i);

    _gradient = move(Gradient());
}

// Least-Square Gradient
// TODO: Rewrite to make it easier to read
vector<Point> PoissonSolver::Gradient()
{
    assert(!_solution.empty());
    vector<Point> gradient(_mesh->tets.size());

    for (auto tet : _mesh->tets)
    {
        double val = _solution[tet->index];
        array<double, 4> adjVal;
        array<Point, 4> dist;
        
        bool neumann = false;

        // TODO: Rename it
        int iN;
        double normGrad;
        Point nN;

        for (int i = 0; i < 4; i++)
        {
            Face* face = tet->faces[i];
            Tet* adjTet = tet->adjTets[i];
            if (_faceBC[face->index].type == PoissonBCType::NonBoundary)
            {
                adjVal[i] = _solution[tet->adjTets[i]->index];
                dist[i] = adjTet->centroid - tet->centroid;
            }
            else if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
            {
                adjVal[i] = _faceBC[face->index].value;
                dist[i] = face->centroid - tet->centroid;
            }
            else if (_faceBC[face->index].type == PoissonBCType::Neumann)
            {
                neumann = true;

                iN = i;
                nN = face->normal;
                normGrad = _faceBC[face->index].normalGrad;
            }
            else if (_faceBC[face->index].type == PoissonBCType::Periodic)
            {
                adjVal[i] = _solution[tet->adjTets[i]->index];
                Point d = (adjTet->centroid - tet->centroid);

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                dist[i] = d + face->centroid - adjTet->faces[k]->centroid;
            }
        }

        array<double, 4> weights;
        for (int i = 0; i < 4; i++)
            weights[i] = 1 / dist[i].Abs();

        if (!neumann)
        {
            // TODO: Store these matrices in memory
            Eigen::Matrix3d m;
            for (int k = 0; k < 3; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    m(k, i) = 0;
                    for (int j = 0; j < 4; j++)
                        m(k, i) += 2 * weights[j] * dist[j].coords[k] * dist[j].coords[i];
                }
            }

            Eigen::Vector3d rhs;
            for (int k = 0; k < 3; k++)
            {
                rhs(k) = 0;
                for (int j = 0; j < 4; j++)
                    rhs(k) += -2 * weights[j] * dist[j].coords[k] * (val - adjVal[j]);
            }

            Eigen::Vector3d grad = m.colPivHouseholderQr().solve(rhs);
            gradient[tet->index] = Point({grad(0), grad(1), grad(2)});
        }
        else
        {
            // TODO: What if there are more than one Neumann face?

            Eigen::Matrix3d m;
            for (int k = 0; k < 2; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    m(k, i) = 0;
                    for (int j = 0; j < 4; j++) {
                        if (j == iN)
                            continue;

                        m(k, i) += 2 * weights[j] * dist[j].coords[k] * dist[j].coords[i];
                    }
                }
            }

            for (int i = 0; i < 3; i++)
                m(2, i) = nN.coords[i];

            Eigen::Vector3d rhs;
            for (int k = 0; k < 2; k++)
            {
                rhs(k) = 0;
                for (int j = 0; j < 4; j++)
                {
                    if (j == iN)
                        continue;

                    rhs(k) += -2 * weights[j] * dist[j].coords[k] * (val - adjVal[j]);
                }
            }

            rhs(2) = normGrad;

            Eigen::Vector3d grad = m.colPivHouseholderQr().solve(rhs);
            gradient[tet->index] = Point({grad(0), grad(1), grad(2)});
        }
    }

    // Smoothing
    Smoother smoother(_mesh);
    smoother.factor = 0.7;
    smoother.nRounds = 5;
    // smoother.SmoothField(gradient);

    return gradient;
}

const vector<double>& PoissonSolver::Potential() const
{
    return _solution;
}

vector<array<double, 3>> PoissonSolver::ElectricField() const
{
    vector<array<double, 3>> field(_mesh->tets.size());
    for (int i = 0; i < _mesh->tets.size(); i++)
        field[i] = (_gradient[i] * (-1)).coords;

    return field;
}
}