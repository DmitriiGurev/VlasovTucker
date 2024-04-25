#include "particle_data.h"
#include "full.h"
#include "tucker.h"

#include <iostream>

using namespace std;

namespace VlasovTucker
{
template class ParticleData<Full>;
template class ParticleData<Tucker>;

template <typename TensorType>
ParticleData<TensorType>::ParticleData(const Mesh* mesh, const VelocityGrid* vGrid) :
    _mesh(mesh), _vGrid(vGrid)
{
    _maxRank = *max_element(_vGrid->nCells.begin(), _vGrid->nCells.end());
    cout << "The maximum rank is " << _maxRank << "\n";
}

template <typename TensorType>
void ParticleData<TensorType>::SetMaxwellPDF(const MaxwellPDF& paramsPDF)
{
    int n0 = _vGrid->nCells[0];
    int n1 = _vGrid->nCells[1];
    int n2 = _vGrid->nCells[2];

    #pragma omp parallel
    {
        vector<TensorType> pdfPrivate;

        #pragma omp for nowait
        for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        {
            Tensor3d t3d(n0, n1, n2);
            if (paramsPDF.temperature != 0.0)
            {
                double normConst = 0;
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
                            t3d(i0, i1, i2) = /* velSquared **/ exp(-mass * velSquared / (2 * boltzConst *
                                paramsPDF.temperature));
                            normConst += t3d(i0, i1, i2) * _vGrid->cellVolume;
                        }
                    } 
                }

                for (int i0 = 0; i0 < n0; i0++)
                {
                    for (int i1 = 0; i1 < n1; i1++)
                    {
                        for (int i2 = 0; i2 < n2; i2++)
                        {
                            t3d(i0, i1, i2) *= paramsPDF.physDensity[tetInd] / normConst;
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
                t3d(i0, i1, i2) = paramsPDF.physDensity[tetInd] / _vGrid->cellVolume;
            }

            TensorType tensor(t3d);
            tensor.Compress(_comprErr, _maxRank);
            pdfPrivate.push_back(tensor);
        }

        #pragma omp critical
        pdf.insert(pdf.end(), pdfPrivate.begin(), pdfPrivate.end());
    }

    double reduction = 0;
    for (int tetInd = 0; tetInd < _mesh->tets.size(); tetInd++)
        reduction += (n0 * n1 * n2) / (double)pdf[tetInd].Size();
    reduction /= (double)_mesh->tets.size();

    // Print compression info
    cout << "PDF size reduction: " << reduction << " times on average\n";
}

template <typename TensorType>
vector<double> ParticleData<TensorType>::Density() const
{
    vector<double> result(_mesh->tets.size());

    // TODO: Do it in parallel
    for (int i = 0; i < _mesh->tets.size(); i++)
        result[i] = pdf[i].Sum() * _vGrid->cellVolume;

    return result;
}

template <typename TensorType>
vector<Vector3d> ParticleData<TensorType>::Velocity() const {
    vector<Vector3d> result(_mesh->tets.size());
    vector<double> density = Density();

    // TODO: Do it in parallel
    for (int i = 0; i < _mesh->tets.size(); i++)
    {
        for (int k = 0; k < 3; k++)
        {
            TensorType vPDF = TensorType(_vGrid->v[k]) * pdf[i];
            if (density[i] != 0) {
                result[i][k] = vPDF.Sum() * _vGrid->cellVolume / density[i];
            }
            else
            {
                result[i][k] = 0;
            }
        }
    }
    return result;
}

template <typename TensorType>
void ParticleData<TensorType>::SetCompressionError(double error)
{
    _comprErr = error;
}

template <typename TensorType>
double ParticleData<TensorType>::CompressionError() const
{
    return _comprErr;
}

template <typename TensorType>
int ParticleData<TensorType>::MaxRank() const
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

double DebyeLength(double temperature, double density, double charge)
{
    return sqrt(epsilon0 * boltzConst * temperature / density) / charge;
}

double PlasmaFrequency(double density, double charge, double mass)
{
    return charge * sqrt(density / (mass * epsilon0));
}

}