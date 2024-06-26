#pragma once

#include "primitives.h"

#include <string>
#include <vector>
#include <unordered_map>

#include <mshio/mshio.h>

namespace VlasovTucker
{
class Mesh
{
public:
    Mesh(std::string mshFile);

    void Reconstruct(double scaleFactor = 1);

    std::unordered_map<int, std::vector<std::string>> BoundaryLabels() const;
    void PrintBoundaryLabels() const;

    void SetPeriodicBounaries(const std::vector<std::array<int, 2>>& periodicPairs);
    std::vector<std::array<int, 2>> PeriodicBoundaries() const;

    const std::unordered_map<int, std::vector<Face*>>& EntityToFaces() const;

    double AverageCellSize() const;

    ~Mesh();

private:
    void _ExtractPoints(double scaleFactor);
    void _ExtractTetsAndFaces();
    void _LabelBoundaryFaces();
    void _FillAdjacencyInfo();
    void _ConfigurePeriodicity();

public:
    std::vector<Point*> points;
    std::vector<Face*> faces;
    std::vector<Tet*> tets;

private:
    mshio::MshSpec _mshSpec;

    std::vector<std::array<int, 2>> _periodicPairs;
    
    // Lookup tables
    std::unordered_map<int, std::vector<std::string>> _entityToPhysGroups;
    std::unordered_map<int, std::vector<Face*>> _entityToFaces;
    std::unordered_map<Triple, Face*, HashTriple> _pointsToFaces;
};
}