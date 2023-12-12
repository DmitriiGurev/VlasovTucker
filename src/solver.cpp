#include "solver.h"
#include "poisson.h"
#include "smoother.h"
#include "vtk.h"
#include "timer.h"

using namespace std;

Solver::Solver(const Mesh* mesh,
               const VelocityGrid<Tensor>* velocityGrid,
               PlasmaParameters* plasmaParameters) :
    _mesh(mesh),
    _vGrid(velocityGrid),
    _plParams(plasmaParameters)
{
    // TODO: Set LogLevel in CMake
    _log = Log(LogLevel::AllText, "solver_");
    
    _PrecomputeNormalTensors();
    _PrecomputeGradCoeffs();

    int memorySpent = _vGrid->nCellsTotal * _mesh->tets.size() * 4 * 8 / pow(10, 6);
    _log << "The normal velocity tensors take up " << memorySpent << " MB of RAM\n"; 
}

void Solver::Solve(int nIterations)
{
    // TODO: Calculate it using the Courant number
    double timeStep = 1.0 * 1.0e-4;

    _log << "Initialize the Poisson solver\n";

    Timer timer;
    PoissonSolver pSolver(_mesh);
    timer.PrintSectionTime("Poisson solver initialization");

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        Timer timer;
        _log << "\n" << "Iteration #" << it << "\n";
        
        _log << Indent(1) << "Compute the electric field\n";
        vector<double> rho = _plParams->Density();
        for (int i = 0; i < rho.size(); i++) 
            rho[i] *= _plParams->charge;

        pSolver.Solve(rho);
        vector<double> phi = pSolver.Potential(); // For debug
        vector<array<double, 3>> field = move(pSolver.ElectricField());
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Update the right-hand side\n";
        vector<Tensor> rhs(_mesh->tets.size());
        for (auto tet : _mesh->tets)
        {
            rhs[tet->index] = Tensor(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
            rhs[tet->index].setZero();
        }

        _log << Indent(2) << "Boltzmann part\n";
        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];

            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;
                
                Tensor upwindPDF = 0.5 * (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] +
                    _plParams->pdf[tetInd]) - _vNormalAbs[tetInd][i] * (_plParams->pdf[adjTetInd] -
                    _plParams->pdf[tetInd]));

                if (fluxScheme == FluxScheme::Upwind)
                {
                    // Upwind reconstruction: X_f = X_C
                    rhs[tetInd] -= tet->faces[i]->area / tet->volume * upwindPDF;
                }
                else 
                {
                    // Compute gradients for a second-order scheme
                    array<Tensor, 3> grad = _Gradient(tet);
                    array<Tensor, 3> adjGrad = _Gradient(tet->adjTets[i]);
                    Point dist = tet->faces[i]->centroid - tet->centroid;
                    Point adjDist = tet->faces[i]->centroid - tet->adjTets[i]->centroid;

                    Tensor correction(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
                    correction.setZero();

                    if (fluxScheme == FluxScheme::UpwindFROMM)
                    {
                        // FROMM reconstruction: X_f = X_C + (grad X)_C * d_Cf
                        for (int k = 0; k < 3; k++)
                        {
                            correction += 0.5 * (_vNormal[tetInd][i] * ((adjDist.coords[k] *
                                adjGrad[k]) + (dist.coords[k] * grad[k])) - _vNormalAbs[tetInd][i] *
                                ((adjDist.coords[k] * adjGrad[k]) - (dist.coords[k] * grad[k])));
                        }

                        rhs[tetInd] -= tet->faces[i]->area / tet->volume * (upwindPDF + correction);
                    }
                    else if (fluxScheme == FluxScheme::UpwindQUICK)
                    {
                        // QUICK reconstruction: X_f = X_C + 0.5 ((grad X)_C + (grad X)_f) * d_Cf
                        array<Tensor, 3> gradFace;
                        double w = adjDist.Abs() / (adjDist.Abs() + dist.Abs());
                        double adjW = dist.Abs() / (adjDist.Abs() + dist.Abs());
                        for (int k = 0; k < 3; k++)
                            gradFace[k] = w * grad[k] + adjW * adjGrad[k];

                        for (int k = 0; k < 3; k++)
                        {
                            correction += 0.5 * (_vNormal[tetInd][i] * ((adjDist.coords[k] *
                                0.5 * (adjGrad[k] + gradFace[k])) + (dist.coords[k] *
                                0.5 * (grad[k] + gradFace[k]))) - _vNormalAbs[tetInd][i] *
                                ((adjDist.coords[k] * 0.5 * (adjGrad[k] + gradFace[k])) -
                                (dist.coords[k] * 0.5 * (grad[k] + gradFace[k]))));
                        }

                        rhs[tetInd] -= tet->faces[i]->area / tet->volume * (upwindPDF + correction);
                    }
                }
            }
        }

        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(2) << "Vlasov part\n";

        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];

            // TODO: Add a constant electric field
            for (int k = 0; k < 3; k++)
            {
                double forceComponent = _plParams->charge * field[tetInd][k];
                rhs[tetInd] -= forceComponent * _PDFDerivative(tet, k);
            }
        }
        
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time integration\n";

        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];
            _plParams->pdf[tetInd] += timeStep * rhs[tetInd];
        }

        timer.PrintSectionTime(Indent(2) + "Done");

        double nSumm = 0.0;
        vector<double> density = _plParams->Density();
        for (auto tet : _mesh->tets)
            nSumm += density[tet->index];

        _log << Indent(1) << "Time: " << it * timeStep << "\n";
        _log << Indent(1) << "Total density: " << nSumm << "\n";

        if (abs(nSumm) > 1e8)
            throw runtime_error("Solution diverged");

        if (it % writeStep == 0)
        {
            string densityFile = "solution/density/density_" + to_string(it / writeStep);
            VTK::WriteCellScalarData(densityFile, *_mesh, density);

            string phiFile = "solution/phi/phi_" + to_string(it / writeStep);
            VTK::WriteCellScalarData(phiFile, *_mesh, phi);

            string fieldFile = "solution/field/e_" + to_string(it / writeStep);
            VTK::WriteCellVectorData(fieldFile, *_mesh, field);

            string distrFile = "solution/distribution/distribution_" + to_string(it / writeStep);
            VTK::WriteDistribution(distrFile, *_vGrid, _plParams->pdf[_mesh->tets.size() / 2]);
        }
    }
}

void Solver::_PrecomputeNormalTensors()
{
    _vNormal = vector<array<Tensor, 4>>(_mesh->tets.size());
    _vNormalAbs = vector<array<Tensor, 4>>(_mesh->tets.size());

    for (auto tet : _mesh->tets)
    {
        int tetInd = tet->index;
        for (int i = 0; i < 4; i++)
        {
            Point normal = tet->faces[i]->normal;
            _vNormal[tetInd][i] = normal.coords[0] * _vGrid->v[0] +
                                  normal.coords[1] * _vGrid->v[1] +
                                  normal.coords[2] * _vGrid->v[2]; 

            _vNormalAbs[tetInd][i] = _vNormal[tetInd][i].abs();
        }
    }
}

Tensor Solver::_PDFDerivative(const Tet* tet, int ind) const
{
    int tetInd = tet->index;

    Tensor pdfDer(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);

    array<int, 3> i;
    for (i[0] = 0; i[0] < _vGrid->nCells[0]; i[0]++)
    {
        for (i[1] = 0; i[1] < _vGrid->nCells[1]; i[1]++)
        {
            for (i[2] = 0; i[2] < _vGrid->nCells[2]; i[2]++)
            {
                array<int, 3> iPlus = {i[0], i[1], i[2]};
                array<int, 3> iMinus = {i[0], i[1], i[2]};

                if (i[ind] == 0)
                {
                    iPlus[ind] = 1;
                    iMinus[ind] = _vGrid->nCells[ind] - 1;
                }
                else if (i[ind] == _vGrid->nCells[ind] - 1)
                {
                    iPlus[ind] = 0;
                    iMinus[ind] = _vGrid->nCells[ind] - 2;
                }
                else
                {
                    iPlus[ind] += 1;
                    iMinus[ind] -= 1;
                }

                pdfDer(i[0], i[1], i[2]) = (
                    _plParams->pdf[tetInd](iPlus[0], iPlus[1], iPlus[2]) - 
                    _plParams->pdf[tetInd](iMinus[0], iMinus[1], iMinus[2])
                    ) / (2 * _vGrid->step[ind]);
            }
        }
    }

    return pdfDer;
}

void Solver::_PrecomputeGradCoeffs()
{
    for (auto tet : _mesh->tets)
    {
        array<Point, 4> dist;
        
        for (int i = 0; i < 4; i++)
        {
            // Remark: Works only for periodic cases
            
            Tet* adjTet = tet->adjTets[i];
            Point d = (adjTet->centroid - tet->centroid);

            int k = 0;
            while (adjTet->adjTets[k] != tet)
                k++;

            dist[i] = d + tet->faces[i]->centroid - adjTet->faces[k]->centroid;

            _distances[{tet, adjTet}] = dist[i];
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
                    m(k, i) += 2 * weights[j] * dist[j].coords[k] * dist[j].coords[i];
            }
        }
        _gradMatrices.push_back(m);
    }
}

array<Tensor, 3> Solver::_Gradient(Tet* tet) const {
    array<Tensor, 3> grad;

    int tetInd = tet->index;
    
    double a1 = _gradMatrices[tetInd](0, 0);
    double a2 = _gradMatrices[tetInd](0, 1);
    double a3 = _gradMatrices[tetInd](0, 2);
    double b1 = _gradMatrices[tetInd](1, 0);
    double b2 = _gradMatrices[tetInd](1, 1);
    double b3 = _gradMatrices[tetInd](1, 2);
    double c1 = _gradMatrices[tetInd](2, 0);
    double c2 = _gradMatrices[tetInd](2, 1);
    double c3 = _gradMatrices[tetInd](2, 2);

    double e1 = a1 * b2 - a2 * b1;
    double e2 = a1 * b3 - a3 * b1;
    double f1 = a1 * c2 - a2 * c1;
    double f2 = a1 * c3 - a3 * c1;

    array<Point, 4> dist;
    array<double, 4> weights;
    for (int i = 0; i < 4; i++)
    {
        dist[i] = _distances.at({tet, tet->adjTets[i]}); 
        weights[i] = 1 / dist[i].Abs();
    }

    array<Tensor, 3> rhs;
    for (int k = 0; k < 3; k++)
    {
        rhs[k] = Tensor(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
        rhs[k].setZero();
    }

    for (int k = 0; k < 3; k++)
    {
        for (int i = 0; i < 4; i++)
        {
            rhs[k] += 2 * weights[i] * dist[i].coords[k] *
                (_plParams->pdf[tet->adjTets[i]->index] - _plParams->pdf[tetInd]);
        }
    }

    Tensor h1 = a1 * rhs[1] - b1 * rhs[0];
    Tensor h2 = a1 * rhs[2] - c1 * rhs[0];

    grad[2] = (h2 * e1 - h1 * f1) / (f2 * e1 - f1 * e2); 
    grad[1] = h1 / e1 - e2 * grad[2] / e1;
    grad[0] = rhs[0] / a1 - a2 * grad[1] / a1 - a3 * grad[2] / a1;

    return grad;
}