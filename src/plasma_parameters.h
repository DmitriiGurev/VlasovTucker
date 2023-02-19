#pragma once

#include <string>
#include <vector>

#include "constants.h"

enum class ParticleType
{
    Electron,
    Ion,
    Neutral
};

struct PlasmaParameters
{
    void ReadParameters(std::string fileName);

    int nSpecies;

    std::vector<ParticleType> pSpecies;
    std::vector<double>       pMass;
    std::vector<double>       pCharge;
    std::vector<double>       pDensity;
    std::vector<double>       pTemperature;
};

void PlasmaParameters::ReadParameters(std::string fileName)
{
    // TODO: Read from a file

    int nSpecies = 1;
    pSpecies     = {ParticleType::Electron};
    pMass        = {elMass, 131.29};
    pCharge      = {-elCharge};
    pDensity     = {1.0e17};
    pTemperature = {1.0 * electronvolt};

    // int nSpecies = 2;
    // pSpecies     = {ParticleType::Electron, ParticleType::Ion};
    // pMass        = {elMass, 131.29 * atomicMass};
    // pCharge      = {-elCharge, elCharge};
    // pDensity     = {1.0e17, 1.0e17};
    // pTemperature = {1.0 * electronvolt, 1.0 * electronvolt};
}