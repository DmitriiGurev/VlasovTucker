#pragma once

#include "mesh.h"

#include <vector>
#include <array>

class Smoother
{
public:
    Smoother(const Mesh* mesh) : _mesh(mesh) {}

    template <typename T>
    void SmoothField(std::vector<T>& field);

public:
    double factor = 0.5;
    double nRounds = 1;

private:
    const Mesh* _mesh;
};

template <typename T>
void Smoother::SmoothField(std::vector<T>& field)
{
    for (int round = 0; round < nRounds; round++) 
    {
        for (auto tet : _mesh->tets)
        {
            T smoothVal;
            smoothVal = field[tet->index] * (1 - factor);

            Tet* adjTet;
            double coeff;
            T adjVal;

            std::vector<double> coeffs;

            adjTet = tet->adjTets[0];
            if (adjTet)
            {
                coeff = 1 / (adjTet->centroid - tet->centroid).Abs();
                coeffs.push_back(coeff);
                adjVal = field[adjTet->index] * coeff;
            }
            for (int j = 1; j < 4; j++)
            {
                Tet* adjTet = tet->adjTets[j];
                if (adjTet)
                {
                    coeff = 1 / (adjTet->centroid - tet->centroid).Abs();
                    coeffs.push_back(coeff);
                    adjVal = adjVal + field[adjTet->index] * coeff;
                }
            }

            double coeffSum = 0;
            for (auto coeff : coeffs)
                coeffSum += coeff;
            
            smoothVal = smoothVal + adjVal * factor / coeffSum;

            field[tet->index] = smoothVal;
        }
    } 
}