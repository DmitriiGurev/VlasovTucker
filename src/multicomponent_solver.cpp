#include "multicomponent_solver.h"
#include "timer.h"
#include "full.h"
#include "tucker.h"

using namespace std;

namespace VlasovTucker
{
template class MulticomponentSolver<Full>;
template class MulticomponentSolver<Tucker>;

template <typename TensorType>
MulticomponentSolver<TensorType>::MulticomponentSolver(Solver<TensorType>* base) : 
    _solvers({base})
{
    _log = Log(LogLevel::Console);
}

template <typename TensorType>
void MulticomponentSolver<TensorType>::AddSolver(Solver<TensorType>* solver)
{
    _solvers.push_back(solver);
}

template <typename TensorType>
void MulticomponentSolver<TensorType>::Solve()
{
    Solver<TensorType>* base = _solvers[0];

    // Set missing step multipliers equal to unity
    for (auto solver : _solvers)
    {
        if (!stepMultipliers.count(solver))
            stepMultipliers[solver] = 1;
    }

    for (auto solver : _solvers)
    {
        solver->timeStep = timeStep * stepMultipliers[solver];
        solver->writeStep = writeStep;
    }

    // Calculate the area of absorbing walls, if any, and set the collected charge to zero
    for (auto solver : _solvers)
        solver->_InitializeWallCharge();

    _log << "Initialize the Poisson solver\n";
    Timer timer;
    base->_poissonSolver.Initialize();
    timer.PrintSectionTime("Poisson solver initialization");

    _log << "Start the main loop\n";
    for (int iteration = 0; iteration < nIterations; iteration++)
    {
        _log << "\n" << "Iteration #" << iteration << "\n";
        _log << "Time: " << iteration * timeStep << "\n";

        _log << Indent(1) << "Compute the electric field\n";
        // Calculate the total charge density
        vector<double> rho(base->_mesh->tets.size(), 0); 
        for (auto solver : _solvers)
        {
            vector<double> density = solver->_pData->Density();
            for (int i = 0; i < rho.size(); i++)
                rho[i] += solver->_pData->charge * density[i];
        }

        // Add the background charge
        if (!base->backgroundChargeDensity.empty())
        {
            for (int i = 0; i < rho.size(); i++)
                rho[i] += base->backgroundChargeDensity[i];
        }

        // Solve the Poisson equation
        base->_poissonSolver.Solve(rho);
        for (auto solver : _solvers)
        {
            solver->_rho = rho;
            solver->_phi = base->_poissonSolver.Potential();
            solver->_field = base->_poissonSolver.ElectricField();
        }
        timer.PrintSectionTime(Indent(1) + "Done");

        for (auto solver : _solvers)
        {
            // Update this PDF once in several iterations
            if (iteration % stepMultipliers[solver])
                continue;

            _log << Indent(1) << "Update the PDF: " << solver->_pData->species << "\n"; 
            solver->_UpdatePDF();
        }
        timer.PrintSectionTime(Indent(1) + "Done");
        
        _log << Indent(1) << "Update the boundary conditions\n";
        // Update the boundary conditions
        unordered_map<int, double> wallCharge = base->_wallCharge;
        for (auto solver : _solvers)
        {
            if (solver == base)
                continue;

            for (auto pair : solver->_wallCharge)
            {
                int boundaryInd = pair.first;
                double charge = pair.second;
                wallCharge[boundaryInd] += charge;
            }
        }

        // Set the updated Neumann condition on the absorbing boundaries 
        for (auto pair : wallCharge)
        {
            int boundaryInd = pair.first;
            double charge = pair.second;
            
            double area = base->_wallArea[boundaryInd];

            FieldBC fieldBC;
            fieldBC.type = FieldBCType::ChargedPlane;
            fieldBC.chargeDensity = charge / area;

            base->SetFieldBC(boundaryInd, fieldBC);
        }

        // Write the results
        if (iteration % writeStep == 0)
        {
            for (auto solver : _solvers)
                solver->_WriteResults(iteration);
        }
    }
}
}