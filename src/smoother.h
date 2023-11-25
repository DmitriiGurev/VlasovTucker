#pragma once

#include "mesh.h"

#include <vector>
#include <array>

class Smoother
{
public:
    Smoother(const Mesh* mesh);

    void SmoothField(std::vector<double>& field);
    void SmoothField(std::vector<std::array<double, 3>>& field);

public:
    double factor = 0.5;
    double nRounds = 1;

private:
    const Mesh* _mesh;
};