#include "poisson.h"

#include <iostream>
#include <algorithm>

using namespace std;

PoissonSolver::PoissonSolver(const Mesh& mesh) :
    _mesh(mesh)
{
    int nEquations = _mesh.tets.size();
    _solutionIsUnique = false;

    for (int i = 0; i < _mesh.faces.size(); i++)
    {
        if (_mesh.faces[i]->type == Internal)
        {
            _faceTypes[i] == BCType::NonBoundary;
        }

        if (_mesh.faces[i]->type == Boundary)
        {
            if (find(_mesh.faces[i]->bcTypes.begin(),
                     _mesh.faces[i]->bcTypes.end(),
                     "Poisson: Dirichlet") !=
                     _mesh.faces[i]->bcTypes.end())
            {
                _faceTypes[i] = BCType::Dirichlet;
                _solutionIsUnique = true;
            }

            if (find(_mesh.faces[i]->bcTypes.begin(),
                        _mesh.faces[i]->bcTypes.end(),
                        "Poisson: Neumann") !=
                        _mesh.faces[i]->bcTypes.end())
            {
                _faceTypes[i] = BCType::Neumann;
            }

            if (find(_mesh.faces[i]->bcTypes.begin(),
                        _mesh.faces[i]->bcTypes.end(),
                        "Poisson: Periodic") !=
                        _mesh.faces[i]->bcTypes.end())
            {
                _faceTypes[i] = BCType::Periodic;
            }
        }
    }
    
    // Assemble the system using Triplets (line, row, value)
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> coeffs;

    _rhs = Eigen::VectorXd::Constant(nEquations, 0.0);

    for (int i = 0; i < _mesh.tets.size(); i++)
    {
        if (!_solutionIsUnique && i == 0)
        {
            coeffs.push_back(Triplet(0, 0, 1));
            continue;
        }

        for (int j = 0; j < 4; j++)
        {
            Tet* tet = _mesh.tets[i];
            Face* face = tet->faces[j];

            if (_faceTypes[face->index] == BCType::NonBoundary)
            {
                Tet* adjTet = tet->adjTets[j];
                Point d = adjTet->centroid - tet->centroid;

                coeffs.push_back(Triplet(
                    i,
                    adjTet->index,
                    (d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()
                    ));
                
                coeffs.push_back(Triplet(
                    i,
                    i,
                    -(d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()
                    ));
            }
            if (_faceTypes[face->index] == BCType::Periodic)
            {
                Tet* adjTet = tet->adjTets[j];
                Point d = (adjTet->centroid - tet->centroid);

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                d = d + face->centroid - adjTet->faces[k]->centroid;

                coeffs.push_back(Triplet(
                    i,
                    adjTet->index,
                    (d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()
                    ));
                
                coeffs.push_back(Triplet(
                    i,
                    i,
                    -(d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()
                    ));
            }
            if (_faceTypes[face->index] == BCType::Dirichlet)
            {
                // TODO
                double dirichletVal = 0.0;

                Point d = face->centroid - tet->centroid;

                coeffs.push_back(Triplet(
                    i,
                    i,
                    -(d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()
                    ));

                _rhs(i) -= (d / d.Abs()).DotProduct(face->normal) * face->area * dirichletVal / d.Abs();
            }
            if (_faceTypes[face->index] == BCType::Neumann)
            {
                // TODO
                Point neumannGrad = {0.0, 0.0, 0.0};
                
                _rhs(i) -= neumannGrad.DotProduct(face->normal) * face->area;
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

vector<double> PoissonSolver::Solve(std::vector<double> rho) const
{
    Eigen::VectorXd rhs = _rhs;
    for (int i = 0; i < _mesh.tets.size(); i++)
        rhs(i) += (-rho[i] / eps0) * _mesh.tets[i]->volume;

    if (!_solutionIsUnique)
        rhs(0) = 0;

    Eigen::VectorXd solution = _solver.solve(rhs);

    vector<double> solutionVec(_mesh.tets.size());
    for (int i = 0; i < _mesh.tets.size(); i++)
        solutionVec[i] = solution(i);

    return solutionVec;
}