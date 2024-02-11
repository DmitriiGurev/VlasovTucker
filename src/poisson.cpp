#include "poisson.h"

#include "log.h"
#include "smoother.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

namespace VlasovTucker
{
PoissonSolver::PoissonSolver(const Mesh* mesh) :
    _mesh(mesh)
{
    int nEquations = _mesh->tets.size();
    _solutionIsUnique = false;

    for (int i = 0; i < _mesh->faces.size(); i++)
    {
        if (_mesh->faces[i]->type == FaceType::Internal)
            _faceBC[i].type == PoissonBCType::NonBoundary;

        if (_mesh->faces[i]->type == FaceType::Boundary)
        {
            auto FaceTypeIs = [this](int faceInd, string type)
            {
                Face* face = _mesh->faces[faceInd]; 
                auto& bcTypes = face->bcTypes;
                return find(bcTypes.begin(), bcTypes.end(), type) != bcTypes.end();
            };

            // TODO: Set boundary values here
            if (FaceTypeIs(i, "Poisson: Dirichlet"))
            {
                _faceBC[i].type = PoissonBCType::Dirichlet;
                _solutionIsUnique = true;
            }
            else if (FaceTypeIs(i, "Poisson: Neumann"))
            {
                _faceBC[i].type = PoissonBCType::Neumann;
            }
            else if (FaceTypeIs(i, "Poisson: Periodic"))
            {
                _faceBC[i].type = PoissonBCType::Periodic;
            }
            else
            {
                throw invalid_argument("Invalid BC type");
            }
        }
    }
    
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
            // else if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
            // {
            //     Point d = face->centroid - tet->centroid;

            //     coeffs.push_back(Triplet(i, i,
            //         -(d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()));

            //     _rhs(i) -= (d / d.Abs()).DotProduct(face->normal) *
            //         face->area * _faceBC[face->index].value / d.Abs();
            // }
            // else if (_faceBC[face->index].type == PoissonBCType::Neumann)
            // {
            //     // TODO: Move this to Solve 
            //     _rhs(i) -= _faceBC[face->index].gradient.DotProduct(face->normal) * face->area;
            // }
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
    auto rho2 = rho;

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

        rhs(i) = (-rho[i] / epsilon0) * _mesh->tets[i]->volume;
    }

    if (_gradient.empty())
    {
        cout << Indent(2) << "_gradient.empty()\n";

        // Initial approximation of the gradient
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
            Tet* adjTet = tet->adjTets[j];

            Point d = face->centroid - tet->centroid;
            Point adjD = face->centroid - adjTet->centroid;

            double g = adjD.Abs() / (adjD.Abs() + d.Abs()); 
            double adjG = d.Abs() / (adjD.Abs() + d.Abs());
            Point grad = _gradient[i];
            Point adjGrad = _gradient[adjTet->index];

            Point e;
            if (_faceBC[face->index].type == PoissonBCType::NonBoundary)
            {
                e = adjTet->centroid - tet->centroid;
                e = e / e.Abs();
            }
            else if (_faceBC[face->index].type == PoissonBCType::Periodic)
            {
                e = adjTet->centroid - tet->centroid;

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                e = e + face->centroid - adjTet->faces[k]->centroid;
                e = e / e.Abs();
            }

            // TODO: Move this up?
            // else if (_faceBC[face->index].type == PoissonBCType::Dirichlet)
            // {
            //     Point d = face->centroid - tet->centroid;

            //     coeffs.push_back(Triplet(i, i,
            //         -(d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()));

            //     _rhs(i) -= (d / d.Abs()).DotProduct(face->normal) *
            //         face->area * _faceBC[face->index].value / d.Abs();
            // }
            // else if (_faceBC[face->index].type == PoissonBCType::Neumann)
            // {
            //     _rhs(i) -= _faceBC[face->index].gradient.DotProduct(face->normal) * face->area;
            // }

            double crossDiff = face->area * (grad * g + adjGrad * adjG).
                DotProduct(face->normal - e / e.DotProduct(face->normal));

            rhs(i) -= crossDiff;
        }
    }

    Eigen::VectorXd solution = _solver.solve(rhs);

    for (int i = 0; i < _mesh->tets.size(); i++)
        _solution[i] = solution(i);

    _gradient = move(Gradient());
}

// Least-Square Gradient
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
        for (int i = 0; i < 4; i++)
        {
            Tet* adjTet = tet->adjTets[i];
            if (_faceBC[tet->faces[i]->index].type == PoissonBCType::NonBoundary)
            {
                adjVal[i] = _solution[tet->adjTets[i]->index];
                dist[i] = adjTet->centroid - tet->centroid;
            }
            else if (_faceBC[tet->faces[i]->index].type == PoissonBCType::Dirichlet)
            {
                adjVal[i] = _faceBC[tet->faces[i]->index].value;
                dist[i] = tet->faces[i]->centroid - tet->centroid;
            }
            else if (_faceBC[tet->faces[i]->index].type == PoissonBCType::Neumann)
            {
                neumann = true;
                gradient[tet->index] = _faceBC[tet->faces[i]->index].gradient;
            }
            else if (_faceBC[tet->faces[i]->index].type == PoissonBCType::Periodic)
            {
                adjVal[i] = _solution[tet->adjTets[i]->index];
                Point d = (adjTet->centroid - tet->centroid);

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                dist[i] = d + tet->faces[i]->centroid - adjTet->faces[k]->centroid;
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