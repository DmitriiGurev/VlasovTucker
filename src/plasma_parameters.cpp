#include "plasma_parameters.h"

#include <iostream>

using namespace std;

template<>
void PlasmaParameters::SetPDF<PlasmaParameters::UnifomMaxwell>(
    const Mesh& mesh,
    const VelocityGrid<Tensor>& velocityGrid,
    const PlasmaParameters::UnifomMaxwell& paramsPDF)
{
    int n0 = velocityGrid.nCells[0];
    int n1 = velocityGrid.nCells[1];
    int n2 = velocityGrid.nCells[2];

    Tensor v(n0, n1, n2);

    double normConst = pow(pMass / (2 * pi * boltzConst * 
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
                    double velJ = velocityGrid.At(i0, i1, i2)[j] - paramsPDF.averageV[j];
                    velSquared += velJ * velJ;
                }
                v(i0, i1, i2) = paramsPDF.physDensity * normConst *
                                exp(-pMass * velSquared /
                                (2 * boltzConst * paramsPDF.temperature));
            }
        } 
    }
    pPDF = vector<Tensor>(mesh.tets.size(), v);
}