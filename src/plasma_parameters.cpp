#include "plasma_parameters.h"

#include <iostream>

using namespace std;

template<>
void PlasmaParameters::SetPDF<PlasmaParameters::Maxwell>(
    const PlasmaParameters::Maxwell& paramsPDF)
{
    int n0 = _vGrid->nCells[0];
    int n1 = _vGrid->nCells[1];
    int n2 = _vGrid->nCells[2];

    Tensor v(n0, n1, n2);
    for (auto tet : _mesh->tets)
    {
        if (paramsPDF.temperature != 0.0)
        {
            double normConst = pow(mass / (2 * pi * boltzConst * 
                                paramsPDF.temperature), 1.5);

            for (int i0 = 0; i0 < n0; i0++)
            {
                for (int i1 = 0; i1 < n1; i1++)
                {
                    for (int i2 = 0; i2 < n2; i2++)
                    {
                        double velSquared = 0;
                        for (int j = 0; j < 3; j++)
                        {
                            double velJ = _vGrid->At(i0, i1, i2)[j] - paramsPDF.averageV[j];
                            velSquared += velJ * velJ;
                        }
                        v(i0, i1, i2) = paramsPDF.physDensity[tet->index] * normConst *
                                        exp(-mass * velSquared /
                                        (2 * boltzConst * paramsPDF.temperature));
                    }
                } 
            }
        }
        
        if (paramsPDF.temperature == 0.0)
        {
            v.setZero();

            int i0 = (paramsPDF.averageV[0] - _vGrid->minV[0]) / _vGrid->step[0];
            int i1 = (paramsPDF.averageV[1] - _vGrid->minV[1]) / _vGrid->step[1];
            int i2 = (paramsPDF.averageV[2] - _vGrid->minV[2]) / _vGrid->step[2];
            v(i0, i1, i2) = paramsPDF.physDensity[tet->index] / _vGrid->cellVolume;
        }

        pdf.push_back(v);
    }
}

vector<double> ConstDensity(const Mesh* mesh, double density)
{
    return vector<double>(mesh->tets.size(), density);
}

vector<double> PlasmaParameters::Density() const
{
    vector<double> result(_mesh->tets.size());

    int n0 = _vGrid->nCells[0];
    int n1 = _vGrid->nCells[1];
    int n2 = _vGrid->nCells[2];

    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        double density = 0;
        for (int i0 = 0; i0 < n0; i0++)
        {
            for (int i1 = 0; i1 < n1; i1++)
            {
                for (int i2 = 0; i2 < n2; i2++)
                {
                    density += pdf[i](i0, i1, i2) * _vGrid->cellVolume;
                }
            }
        }
        result[i] = density;
    }
    return result;
}