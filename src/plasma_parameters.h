#pragma once

#include <string>
#include <vector>

#include "constants.h"
#include "mesh.h"
#include "tensor_type.h"
#include "velocity_grid.h"

enum class ParticleType
{
    Electron,
    Ion,
    Neutral
};

class PlasmaParameters
{
public:
    // PlasmaParameters(int nSpecies) :
    //     nSpecies(nSpecies) {}
    PlasmaParameters() {}

    struct UnifomMaxwell
    {
        double physDensity;
        double temperature;
        std::array<double, 3> averageV;
    };

    template <typename ParamsPDF>
    void SetPDF(const Mesh& mesh,
                const VelocityGrid<Tensor>& velocityGrid,
                const ParamsPDF& parameters);

public:
    // const int nSpecies;
    
    // std::vector<ParticleType> pSpecies;
    // std::vector<double>       pMass;
    // std::vector<double>       pCharge;
    // std::vector<std::vector<Tensor>> pPDF;

    ParticleType        pSpecies;
    double              pMass;
    double              pCharge;
    std::vector<Tensor> pPDF;
};