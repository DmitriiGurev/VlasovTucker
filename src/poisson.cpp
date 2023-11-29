#include "poisson.h"

#include "log.h"
#include "smoother.h"

#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

PoissonSolver::PoissonSolver(const Mesh* mesh) :
    _mesh(mesh)
{
    int nEquations = _mesh->tets.size();
    _solutionIsUnique = false;

    for (int i = 0; i < _mesh->faces.size(); i++)
    {
        if (_mesh->faces[i]->type == Internal)
        {
            _faceTypes[i] == BCType::NonBoundary;
        }

        if (_mesh->faces[i]->type == Boundary)
        {
            if (find(_mesh->faces[i]->bcTypes.begin(),
                     _mesh->faces[i]->bcTypes.end(),
                     "Poisson: Dirichlet") !=
                     _mesh->faces[i]->bcTypes.end())
            {
                _faceTypes[i] = BCType::Dirichlet;
                _solutionIsUnique = true;
            }

            if (find(_mesh->faces[i]->bcTypes.begin(),
                        _mesh->faces[i]->bcTypes.end(),
                        "Poisson: Neumann") !=
                        _mesh->faces[i]->bcTypes.end())
            {
                _faceTypes[i] = BCType::Neumann;
            }

            if (find(_mesh->faces[i]->bcTypes.begin(),
                        _mesh->faces[i]->bcTypes.end(),
                        "Poisson: Periodic") !=
                        _mesh->faces[i]->bcTypes.end())
            {
                _faceTypes[i] = BCType::Periodic;
            }
        }
    }
    
    // Assemble the system using Triplets (line, row, value)
    typedef Eigen::Triplet<double> Triplet;
    vector<Triplet> coeffs;

    _rhs = Eigen::VectorXd::Constant(nEquations, 0.0);

    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        if (!_solutionIsUnique && i == 0)
        {
            coeffs.push_back(Triplet(0, 0, 1));
            continue;
        }

        for (int j = 0; j < 4; j++)
        {
            Tet* tet = _mesh->tets[i];
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
            else if (_faceTypes[face->index] == BCType::Periodic)
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
            else if (_faceTypes[face->index] == BCType::Dirichlet)
            {
                // TODO
                double dirichletVal = 0.0;

                Point d = face->centroid - tet->centroid;

                coeffs.push_back(Triplet(
                    i,
                    i,
                    -(d / d.Abs()).DotProduct(face->normal) * face->area / d.Abs()
                    ));

                _rhs(i) -= (d / d.Abs()).DotProduct(face->normal) *
                    face->area * dirichletVal / d.Abs();
            }
            else if (_faceTypes[face->index] == BCType::Neumann)
            {
                // TODO
                Point neumannGrad = Point({0.0, 0.0, 0.0});
                
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

vector<double> PoissonSolver::Solve(vector<double> rho) const
{
    // (?)
    if (!_solutionIsUnique)
    {
        // Make the total charge equal to zero
        double sum = 0;
        double vol = 0;
        for (int i = 0; i < _mesh->tets.size(); i++)
        {
            sum += rho[i] * _mesh->tets[i]->volume;
            vol += _mesh->tets[i]->volume;
        }
        // cout << "Sum = " << sum << "\n";

        // double sum2 = 0;
        for (int i = 0; i < _mesh->tets.size(); i++)
        {
            rho[i] -= sum / vol;
            // sum2 += rho[i] * _mesh->tets[i]->volume;
        }
        // cout << "Sum 2 = " << sum2 << "\n";
        // rho[1] -= sum2 / _mesh->tets[1]->volume;
    }

    Eigen::VectorXd rhs = _rhs;
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        rhs(i) += (-rho[i] / eps0) * _mesh->tets[i]->volume;
    }

    if (!_solutionIsUnique)
    {
        cout << Indent(2) << "Solution is not unique\n";
        rhs(0) = 0;        
    }

    Eigen::VectorXd solution = _solver.solve(rhs);

    vector<double> solutionVec(_mesh->tets.size());
    for (int i = 0; i < _mesh->tets.size(); i++)
        solutionVec[i] = solution(i);

    return solutionVec;
}

vector<array<double, 3>> PoissonSolver::ElectricField(const vector<double>& potential) const
{
    vector<array<double, 3>> result(_mesh->tets.size());

    for (auto tet : _mesh->tets)
    {
        // TODO: Add other bc types
        array<Point, 4> dist;
        for (int i = 0; i < 4; i++)
        {
            Tet* adjTet = tet->adjTets[i];
            if (_faceTypes[tet->faces[i]->index] == BCType::NonBoundary)
            {
                dist[i] = adjTet->centroid - tet->centroid;
            }
            if (_faceTypes[tet->faces[i]->index] == BCType::Periodic)
            {
                Point d = (adjTet->centroid - tet->centroid);

                int k = 0;
                while (adjTet->adjTets[k] != tet)
                    k++;

                dist[i] = d + tet->faces[i]->centroid - adjTet->faces[k]->centroid;
            }
        }

        double val = potential[tet->index];
        array<double, 4> adjVal;
        for (int i = 0; i < 4; i++)
            adjVal[i] = potential[tet->adjTets[i]->index];

        Eigen::Matrix3d m;
        for (int k = 0; k < 3; k++)
        {
            for (int i = 0; i < 3; i++)
            {
                m(k, i) = 0;
                for (int j = 0; j < 4; j++)
                    m(k, i) += 2 * dist[j].coords[k] * dist[j].coords[i];
            }
        }
        Eigen::Vector3d rhs;
        for (int k = 0; k < 3; k++)
        {
            rhs(k) = 0;
            for (int j = 0; j < 4; j++)
                rhs(k) += -2 * dist[j].coords[k] * (val - adjVal[j]);
        }

        Eigen::Vector3d grad = m.colPivHouseholderQr().solve(rhs);
        result[tet->index] = {-grad(0), -grad(1), -grad(2)};
    }

    // Smoothing (?)
    Smoother smoother(_mesh);
    smoother.factor = 0.7;
    smoother.nRounds = 5;
    smoother.SmoothField(result);

    return result;
}