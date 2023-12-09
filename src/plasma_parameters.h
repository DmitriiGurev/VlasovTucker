#pragma once

#include <string>
#include <vector>
#include <functional>

#include "constants.h"
#include "mesh.h"
#include "typedefs.h"
#include "velocity_grid.h"

enum class ParticleType
{
    Electron,
    Ion,
    Neutral,
    Custom
};

class PlasmaParameters
{
public:
    PlasmaParameters(const Mesh* mesh, const VelocityGrid<Tensor>* vGrid) :
        _mesh(mesh), _vGrid(vGrid) {}

    struct MaxwellPDF
    {
        std::vector<double>   physDensity;
        double                temperature;
        std::array<double, 3> mostProbableV;
    };

    template <typename ParamsPDF>
    void SetPDF(const ParamsPDF& parameters);

    std::vector<double> Density() const;

public:
    ParticleType species;
    double       mass;
    double       charge;

    std::vector<Tensor> pdf;

private:
    const Mesh* _mesh;
    const VelocityGrid<Tensor>* _vGrid;
};

// TODO: Move it somewhere else?
std::vector<double> ScalarField(const Mesh* mesh,
    std::function<double(const Point&)> densityFunc);