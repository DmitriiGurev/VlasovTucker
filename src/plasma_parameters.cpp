#include "plasma_parameters.h"
#include "full.h"
#include "tucker.h"

#include <iostream>

using namespace std;

namespace VlasovTucker
{
template class PlasmaParameters<Full>;
template class PlasmaParameters<Tucker>;

template <typename TensorType>
PlasmaParameters<TensorType>::PlasmaParameters(const Mesh* mesh, const VelocityGrid* vGrid) :
    _mesh(mesh), _vGrid(vGrid)
{
    _maxRank = *max_element(_vGrid->nCells.begin(), _vGrid->nCells.end());
    cout << "The maximum rank is " << _maxRank << "\n";
}

template <typename TensorType>
void PlasmaParameters<TensorType>::SetMaxwellPDF(const MaxwellPDF& paramsPDF)
{
    int n0 = _vGrid->nCells[0];
    int n1 = _vGrid->nCells[1];
    int n2 = _vGrid->nCells[2];

    Tensor3d t3d(n0, n1, n2);

    double reduction = 0;

    // TODO: Do it in parallel
    for (auto tet : _mesh->tets)
    {
        if (paramsPDF.temperature != 0.0)
        {
            double normConst = pow(mass / (2 * pi * boltzConst * paramsPDF.temperature), 1.5);
            for (int i0 = 0; i0 < n0; i0++)
            {
                for (int i1 = 0; i1 < n1; i1++)
                {
                    for (int i2 = 0; i2 < n2; i2++)
                    {
                        double velSquared = 0;
                        for (int j = 0; j < 3; j++)
                        {
                            double velJ = _vGrid->At(i0, i1, i2)[j] - paramsPDF.mostProbableV[j];
                            velSquared += velJ * velJ;
                        }
                        t3d(i0, i1, i2) = paramsPDF.physDensity[tet->index] * normConst *
                                        exp(-mass * velSquared /
                                        (2 * boltzConst * paramsPDF.temperature));
                    }
                } 
            }
        }
        
        if (paramsPDF.temperature == 0.0)
        {
            t3d.setZero();
            int i0 = (paramsPDF.mostProbableV[0] - _vGrid->minV[0]) / _vGrid->step[0];
            int i1 = (paramsPDF.mostProbableV[1] - _vGrid->minV[1]) / _vGrid->step[1];
            int i2 = (paramsPDF.mostProbableV[2] - _vGrid->minV[2]) / _vGrid->step[2];
            t3d(i0, i1, i2) = paramsPDF.physDensity[tet->index] / _vGrid->cellVolume;
        }

        TensorType tensor(t3d);
        tensor.Compress(_comprErr, _maxRank);
        pdf.push_back(tensor);

        reduction += t3d.size() / (double)tensor.Size();
        
    }
    reduction /= (double)_mesh->tets.size();

    // Print compression info
    cout << "PDF size reduction: " << reduction << " times on average\n";
}

template <typename TensorType>
vector<double> PlasmaParameters<TensorType>::Density() const
{
    vector<double> result(_mesh->tets.size());

    // TODO: Do it in parallel
    for (int i = 0; i < _mesh->tets.size(); i++)
        result[i] = pdf[i].Sum() * _vGrid->cellVolume;

    return result;
}

template <typename TensorType>
void PlasmaParameters<TensorType>::SetCompressionError(double error)
{
    _comprErr = error;
}

template <typename TensorType>
double PlasmaParameters<TensorType>::CompressionError() const
{
    return _comprErr;
}

template <typename TensorType>
int PlasmaParameters<TensorType>::MaxRank() const
{
    return _maxRank;
}

vector<double> ScalarField(const Mesh* mesh, function<double(const Point&)> densityFunc)
{
    vector<double> density(mesh->tets.size());
    for (int i = 0; i < mesh->tets.size(); i++)
        density[i] = densityFunc(mesh->tets[i]->centroid);

    return density;
}
}