#pragma once

#include <string>
#include <vector>
#include <functional>

#include "constants.h"
#include "mesh.h"
#include "typedefs.h"
#include "velocity_grid.h"

namespace VlasovTucker
{
enum class ParticleType { Electron, Ion, Neutral, Custom };

struct MaxwellPDF
{
    std::vector<double> physDensity;
    double temperature;
    Vector3d mostProbableV;
};

template <typename TensorType>
class PlasmaParameters
{
public:
    PlasmaParameters(const Mesh* mesh, const VelocityGrid* vGrid);

    void SetMaxwellPDF(const MaxwellPDF& paramsPDF);

    std::vector<double> Density() const;
    std::vector<Vector3d> Velocity() const;

    void SetCompressionError(double error);
    double CompressionError() const;
    int MaxRank() const;

public:
    ParticleType species;
    double mass;
    double charge;

    // Particle distribution function
    std::vector<TensorType> pdf;

private:
    const Mesh* _mesh;
    const VelocityGrid* _vGrid;

    // Compression error
    double _comprErr = 1e-10;
    // Maximum tensor rank 
    int _maxRank;
};

// TODO: Move it somewhere else?
std::vector<double> ScalarField(const Mesh* mesh,
    std::function<double(const Point&)> densityFunc);

double DebyeLength(double temperature, double density, double charge);
double PlasmaFrequency(double density, double charge, double mass);
}