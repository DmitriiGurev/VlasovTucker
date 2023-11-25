#include "smoother.h"

using namespace std;

Smoother::Smoother(const Mesh* mesh) : _mesh(mesh) {}

void Smoother::SmoothField(vector<double>& field)
{
    for (int round = 0; round < nRounds; round++) 
    {
        for (auto tet : _mesh->tets)
        {
            double smoothVal;
            smoothVal = field[tet->index] * (1 - factor);

            vector<double> coeffs;
            double adjVal = 0;
            for (int j = 0; j < 4; j++)
            {
                Tet* adjTet = tet->adjTets[j];
                double coeff = 1 / (adjTet->centroid - tet->centroid).Abs();
                adjVal += field[adjTet->index] * coeff;
                coeffs.push_back(coeff);
            }

            double coeffSum = 0;
            for (auto coeff : coeffs)
                coeffSum += coeff;
        
            smoothVal += factor * adjVal / coeffSum;
            field[tet->index] = smoothVal;
        }
    }
}

void Smoother::SmoothField(vector<array<double, 3>>& field)
{
    for (int round = 0; round < nRounds; round++) 
    {
        for (auto tet : _mesh->tets)
        {
            array<double, 3> smoothVec;
            for (int i = 0; i < 3; i++)
                smoothVec[i] = field[tet->index][i] * (1 - factor);

            vector<double> coeffs;
            array<double, 3> adjVec = {0, 0, 0};
            for (int j = 0; j < 4; j++)
            {
                Tet* adjTet = tet->adjTets[j];

                double coeff = 1 / (adjTet->centroid - tet->centroid).Abs();
                for (int i = 0; i < 3; i++)
                    adjVec[i] += field[adjTet->index][i] * coeff;

                coeffs.push_back(coeff);
            }

            double coeffSum = 0;
            for (auto coeff : coeffs)
                coeffSum += coeff;
            
            for (int i = 0; i < 3; i++)
                smoothVec[i] += factor * adjVec[i] / coeffSum;

            field[tet->index] = smoothVec;
        }
    } 
}