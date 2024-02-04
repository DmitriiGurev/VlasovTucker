#include "solver.h"
#include "poisson.h"
#include "smoother.h"
#include "vtk.h"
#include "timer.h"

using namespace std;

namespace VlasovTucker
{
Solver::Solver(const Mesh* mesh,
               const VelocityGrid<Tensor3d>* velocityGrid, // TODO: Change to Full/Tucker
               PlasmaParameters* plasmaParameters) :
    _mesh(mesh),
    _vGrid(velocityGrid),
    _plParams(plasmaParameters)
{
    // TODO: Set LogLevel in CMake
    _log = Log(LogLevel::AllText, "solver_");

    _maxRank = *max_element(_vGrid->nCells.begin(), _vGrid->nCells.end());
    cout << "The maximum rank is " << _maxRank << "\n";

    // Set boundary conditions
    // _boundaryConditions = std::vector<ParticleBC>(_mesh->faces.size());
    
    _PrecomputeNormalTensors();

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

    vector<Tucker> comprPDF;
    for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        comprPDF.push_back(Tucker(_plParams->pdf[tetInd], comprPrecision, _maxRank));

    // Print compression info
    cout << "PDF size reduction: " << _plParams->pdf[0].size() << " vs " << comprPDF[0].Size();
    double reduction = 0;
    for (auto tet : _mesh->tets)
        reduction += _plParams->pdf[0].size() / (double)comprPDF[0].Size();
    reduction /= _mesh->tets.size();
    cout << " (" << reduction << " times on average)\n";
    cout << "Normal speed size reduction: " << _vNormal[0][0].size() << " vs " <<
        _comprVNormal[0][0].Size();
    reduction = 0;
    for (auto tet : _mesh->tets)
        reduction += _vNormal[0][0].size() / (double)_comprVNormal[0][0].Size();
    reduction /= _mesh->tets.size();
    cout << " (" << reduction << " times on average)\n";
    cout << "Absolute normal speed size reduction: " << _vNormalAbs[0][0].size() << " vs " <<
        _comprVNormalAbs[0][0].Size();
    reduction = 0;
    for (auto tet : _mesh->tets)
        reduction += _vNormalAbs[0][0].size() / (double)_comprVNormalAbs[0][0].Size();
    reduction /= _mesh->tets.size();
    cout << " (" << reduction << " times on average)\n";

    _log << "Start the main loop\n";
    for (int it = 0; it < nIterations; it++)
    {
        Timer timer;
        _log << "\n" << "Iteration #" << it << "\n";
        
        _log << Indent(1) << "Compute the electric field\n";

        // Remark: Temporary
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
            _plParams->pdf[tetInd] = comprPDF[tetInd].Reconstructed();

        vector<double> rho = _plParams->Density();
        for (int i = 0; i < rho.size(); i++) 
            rho[i] *= _plParams->charge;

        pSolver.Solve(rho);
        // vector<double> phi = pSolver.Potential(); // For debug
        vector<array<double, 3>> field = move(pSolver.ElectricField());
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Update the right-hand side\n";
        vector<Tensor3d> rhs(_mesh->tets.size());
        vector<Tucker> comprRHS;
        for (auto tet : _mesh->tets)
        {
            rhs[tet->index] = Tensor3d(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);
            rhs[tet->index].setZero();

            comprRHS.push_back(Tucker(rhs[tet->index], comprPrecision, _maxRank));
        }

        _log << Indent(2) << "Boltzmann part\n";
        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tet* tet = _mesh->tets[tetInd];

            for (int i = 0; i < 4; i++)
            {
                int adjTetInd = tet->adjTets[i]->index;
                
                // Tensor upwindPDF = 0.5 * (_vNormal[tetInd][i] * (_plParams->pdf[adjTetInd] +
                //     _plParams->pdf[tetInd]) - _vNormalAbs[tetInd][i] * (_plParams->pdf[adjTetInd] -
                //     _plParams->pdf[tetInd]));

                // Upwind reconstruction: X_f = X_C
                // rhs[tetInd] -= tet->faces[i]->area / tet->volume * upwindPDF;

                comprRHS[tetInd] -= tet->faces[i]->area / tet->volume *
                    0.5 * (_comprVNormal[tetInd][i] * (comprPDF[adjTetInd] +
                    comprPDF[tetInd]) - _comprVNormalAbs[tetInd][i] *
                    (comprPDF[adjTetInd] - comprPDF[tetInd]));

                comprRHS[tetInd].Compress(comprPrecision, _maxRank);
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
                double forceComponent = (_plParams->charge / _plParams->mass) * field[tetInd][k];
                // rhs[tetInd] -= forceComponent * _PDFDerivative(tet, k); // TODO: /m?

                comprRHS[tetInd] -= forceComponent * _PDFDerivative(tet, k);
            }
            comprRHS[tetInd].Compress(comprPrecision, _maxRank);
        }
        
        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time integration\n";

        #pragma omp parallel for 
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            // Tet* tet = _mesh->tets[tetInd];
            // _plParams->pdf[tetInd] += timeStep * rhs[tetInd];

            comprPDF[tetInd] += timeStep * comprRHS[tetInd];
            comprPDF[tetInd].Compress(comprPrecision, _maxRank);
        }

        // _log << Indent(1) << "Difference: " <<
        //     (_plParams->pdf[0] - comprPDF[0].Reconstructed()).square().sum().sqrt() << "\n";

        timer.PrintSectionTime(Indent(2) + "Done");

        _log << Indent(1) << "Time: " << it * timeStep << "\n";

        array<double, 3> averageR = {0, 0, 0};
        for (int i = 0; i < _mesh->tets.size(); i++)
        {
            for (int j = 0; j < 3; j++)
                averageR[j] += comprPDF[i].Ranks()[j] / (double)_mesh->tets.size();
        }
        cout << "Average ranks = {" <<
            averageR[0] << ", " <<
            averageR[1] << ", " <<
            averageR[2] << "}\n";

        if (it % writeStep == 0)
        {
            for (int i = 0; i < _mesh->tets.size(); i++)
            {
                _plParams->pdf[i] = comprPDF[i].Reconstructed();
            }
            vector<double> density = _plParams->Density();

            for (auto d : density)
                assert(d == d);

            string densityFile = "solution/density/density_" + to_string(it / writeStep);
            WriteCellScalarDataVTK(densityFile, *_mesh, density);

            // string phiFile = "solution/phi/phi_" + to_string(it / writeStep);
            // WriteCellScalarDataVTK(phiFile, *_mesh, phi);

            // string fieldFile = "solution/field/e_" + to_string(it / writeStep);
            // WriteCellVectorDataVTK(fieldFile, *_mesh, field);

            string distrFile = "solution/distribution/distribution_" + to_string(it / writeStep);
            // WriteDistributionVTK(distrFile, *_vGrid, _plParams->pdf[_mesh->tets.size() / 2]);
            WriteDistributionVTK(distrFile, *_vGrid, comprPDF[_mesh->tets.size() / 2].Reconstructed());

            // string analytFile = "solution/analytical/analytical_" + to_string(it / writeStep);
            // auto rhoFunc = [it, timeStep](const Point& p)
            // { 
            //     return 10 + 0.2 * sin(1 * (p.coords[0] - it * timeStep) * (2 * pi));
            //     // return 10 + 0.2 * (((p.coords[0] - it * timeStep < 0.2) &&
            //     //     (p.coords[0] - it * timeStep > 0.0)) ? 1 : 0);
            // };
            // WriteCellScalarDataVTK(analytFile, *_mesh, ScalarField(_mesh, rhoFunc));
        }
    }
}

void Solver::_PrecomputeNormalTensors()
{
    _vNormal = vector<array<Tensor3d, 4>>(_mesh->tets.size());
    _vNormalAbs = vector<array<Tensor3d, 4>>(_mesh->tets.size());

    _comprVNormal = vector<array<Tucker, 4>>(_mesh->tets.size());
    _comprVNormalAbs = vector<array<Tucker, 4>>(_mesh->tets.size());

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

            _comprVNormal[tetInd][i] = Tucker(_vNormal[tetInd][i], comprPrecision);
            _comprVNormalAbs[tetInd][i] = Tucker(_vNormalAbs[tetInd][i], comprPrecision, 6);

            // WriteDistributionVTK("uncompr", *_vGrid, _vNormalAbs[tetInd][i]);
            // WriteDistributionVTK("compr", *_vGrid, _comprVNormalAbs[tetInd][i].Reconstructed());
            // exit(1);
        }
    }
}

// Tensor Solver::_PDFDerivative(const Tet* tet, int ind) const
Tucker Solver::_PDFDerivative(const Tet* tet, int ind) const
{
    int tetInd = tet->index;

    // Tensor pdfDer(_vGrid->nCells[0], _vGrid->nCells[1], _vGrid->nCells[2]);

    // array<int, 3> i;
    // for (i[0] = 0; i[0] < _vGrid->nCells[0]; i[0]++)
    // {
    //     for (i[1] = 0; i[1] < _vGrid->nCells[1]; i[1]++)
    //     {
    //         for (i[2] = 0; i[2] < _vGrid->nCells[2]; i[2]++)
    //         {
    //             array<int, 3> iPlus = {i[0], i[1], i[2]};
    //             array<int, 3> iMinus = {i[0], i[1], i[2]};

    //             if (i[ind] == 0)
    //             {
    //                 iPlus[ind] = 1;
    //                 iMinus[ind] = _vGrid->nCells[ind] - 1;
    //             }
    //             else if (i[ind] == _vGrid->nCells[ind] - 1)
    //             {
    //                 iPlus[ind] = 0;
    //                 iMinus[ind] = _vGrid->nCells[ind] - 2;
    //             }
    //             else
    //             {
    //                 iPlus[ind] += 1;
    //                 iMinus[ind] -= 1;
    //             }

    //             pdfDer(i[0], i[1], i[2]) = (
    //                 _plParams->pdf[tetInd](iPlus[0], iPlus[1], iPlus[2]) - 
    //                 _plParams->pdf[tetInd](iMinus[0], iMinus[1], iMinus[2])
    //                 ) / (2 * _vGrid->step[ind]);
    //         }
    //     }
    // }
 
    Tucker pdfCompr(_plParams->pdf[tetInd]);

    auto core = pdfCompr.Core();
    auto u = pdfCompr.U();

    u[ind] = (_vGrid->d[ind] * u[ind]).eval();

    return Tucker(core, u);
}
}