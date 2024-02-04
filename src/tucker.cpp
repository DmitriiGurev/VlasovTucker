#include "tucker.h"

#include <iostream>
#include <cmath>

#include <Eigen/SVD>
#include <unsupported/Eigen/KroneckerProduct>

using namespace Eigen;
using std::vector;

namespace VlasovTucker
{
Tucker::Tucker() :
    _n({0, 0, 0}),
    _r({0, 0, 0})
{}

// Zero tensor with given ranks
Tucker::Tucker(int n0, int n1, int n2,
               int r0, int r1, int r2) :
    _n({n0, n1, n2}),
    _r({r0, r1, r2})
{
    _u[0] = MatrixXd::Zero(n0, r0);
    _u[1] = MatrixXd::Zero(n1, r1);
    _u[2] = MatrixXd::Zero(n2, r2);

    _core = Tensor3d(r0, r1, r2);
    _core.setZero();
}

// Compress a tensor with a given accuracy
Tucker::Tucker(const Tensor3d& tensor,
               double precision,
               int maxRank)
{
    _n = {(int)tensor.dimension(0),
          (int)tensor.dimension(1),
          (int)tensor.dimension(2)};

    _ComputeU(tensor, precision, maxRank);

    _r = {(int)_u[0].cols(),
          (int)_u[1].cols(),
          (int)_u[2].cols()};

    MatrixXd S0 = _u[0].transpose() * Unfolding(tensor, 0) * kroneckerProduct(_u[1], _u[2]);
    _core = Folding(_r[0], _r[1], _r[2], S0, 0);
}

// Create a tensor with given core and matrices U
Tucker::Tucker(const Tensor3d& core,
               const std::array<MatrixXd, 3>& u) :
    _core(core), _u(u) 
{
    _n = {(int)_u[0].rows(),
          (int)_u[1].rows(),
          (int)_u[2].rows()};

    _r = {(int)_core.dimension(0),
          (int)_core.dimension(1),
          (int)_core.dimension(2)};
}

Tucker& Tucker::Compress(double precision, int maxRank)
{
    vector<MatrixXd> Q(3);
    vector<MatrixXd> R(3);
    for (int i : {0, 1, 2})
    {
        MatrixXd thinQ(_u[i].rows(), _u[i].cols());
        MatrixXd q(_u[i].rows(), _u[i].rows());

        HouseholderQR<MatrixXd> householderQR(_u[i]);
        q = householderQR.householderQ();
        thinQ.setIdentity();

        Q[i] = householderQR.householderQ() * thinQ;
        R[i] = Q[i].transpose() * _u[i];
    }

    MatrixXd A0 = R[0] * Unfolding(_core, 0) * kroneckerProduct(R[1], R[2]).transpose();
    Tensor3d aux = Folding(_r[0], _r[1], _r[2], A0, 0);

    Tucker auxTucker(aux, precision, maxRank);

    _core = auxTucker.Core();

    for (int i : {0, 1, 2})
        _u[i] = Q[i] * auxTucker.U()[i];

    _r = {(int)_u[0].cols(),
          (int)_u[1].cols(),
          (int)_u[2].cols()};

    return *this;
}

int Tucker::Size() const
{
    return _core.size() + _u[0].size() + _u[1].size() + _u[2].size(); 
}

array<int, 3> Tucker::Dimensions() const
{
    return _n;
}

array<int, 3> Tucker::Ranks() const
{
    return _r;
}

array<MatrixXd, 3> Tucker::U() const
{
    return _u;
}

Tensor3d Tucker::Core() const
{
    return _core;
}

double Tucker::operator()(int i0, int i1, int i2) const
{
    double el = 0;

    for (int j0 = 0; j0 < _r[0]; j0++)
    {
        for (int j1 = 0; j1 < _r[1]; j1++)
        {
            for (int j2 = 0; j2 < _r[2]; j2++)
            {
                el += _core(j0, j1, j2) * _u[0](i0, j0) * _u[1](i1, j1) * _u[2](i2, j2);
            }
        }
    }
    return el;
}

Tensor3d Tucker::Reconstructed() const
{
    return Folding(_n[0], _n[1], _n[2],
                   _u[0] * Unfolding(_core, 0) * kroneckerProduct(_u[1], _u[2]).transpose(),
                   0);
}

double Tucker::Sum() const
{
    double sum = 0;
    for (int j0 = 0; j0 < _n[0]; j0++)
    {
        for (int j1 = 0; j1 < _n[1]; j1++)
        {
            for (int j2 = 0; j2 < _n[2]; j2++)
            {
                sum += (*this)(j0, j1, j2);
            }
        }
    }
    return sum;
}

double Tucker::Norm() const
{
    double sumSq = 0;
    for (int j0 = 0; j0 < _n[0]; j0++)
    {
        for (int j1 = 0; j1 < _n[1]; j1++)
        {
            for (int j2 = 0; j2 < _n[2]; j2++)
            {
                double el = (*this)(j0, j1, j2);
                sumSq += el * el;
            }
        }
    }
    return sqrt(sumSq);
}

std::ostream& operator<<(std::ostream& out, const Tucker& tucker)
{
    out << "This is a 3D tensor in the Tucker format with \n";
    out << "r0 = " << tucker._r[0] << ", n0 = " << tucker._n[0] << "\n";
    out << "r1 = " << tucker._r[1] << ", n1 = " << tucker._n[1] << "\n";
    out << "r2 = " << tucker._r[2] << ", n2 = " << tucker._n[2];
    return out;
}

Tucker operator+(const Tucker& t1, const Tucker& t2)
{
    if (t1.Dimensions() != t2.Dimensions())
        throw std::invalid_argument("Different shapes in sum");

    Tucker result(t1._n[0],
                  t1._n[1],
                  t1._n[2],
                  t1._r[0] + t2._r[0],
                  t1._r[1] + t2._r[1],
                  t1._r[2] + t2._r[2]);

    for (int k0 = 0; k0 < t1._r[0] + t2._r[0]; k0++)
    {
        for (int k1 = 0; k1 < t1._r[1] + t2._r[1]; k1++)
        {
            for (int k2 = 0; k2 < t1._r[2] + t2._r[2]; k2++)
            {
                if (k0 < t1._r[0] && 
                    k1 < t1._r[1] &&
                    k2 < t1._r[2])
                    result._core(k0, k1, k2) = t1._core(k0, k1, k2);
                    
                if (k0 >= t1._r[0] &&
                    k1 >= t1._r[1] &&
                    k2 >= t1._r[2])
                    result._core(k0, k1, k2) = t2._core(k0 - t1._r[0], k1 - t1._r[1], k2 - t1._r[2]);
            }
        }
    }

    for (int i : {0, 1, 2})
    {
        result._u[i].resize(t1._n[i], t1._r[i] + t2._r[i]);
        result._u[i].block(0, 0, t1._n[i], t1._r[i]) = t1._u[i];
        result._u[i].block(0, t1._r[i], t1._n[i], t2._r[i]) = t2._u[i];
    }
    return result;
}

Tucker& Tucker::operator+=(const Tucker& t)
{
    *this = *this + t;
    return *this;
}

Tucker& Tucker::operator-=(const Tucker& t)
{
    *this = *this - t;
    return *this;
}

Tucker& Tucker::operator*=(const Tucker& t)
{
    *this = *this * t;
    return *this;
}

Tucker& Tucker::operator*=(double d)
{
    *this = *this * d;
    return *this;
}

Tucker operator-(const Tucker& t1, const Tucker& t2)
{
    return t1 + (-1.0) * t2;
}

Tucker operator*(const Tucker& t1, const Tucker& t2) 
{
    if (t1.Dimensions() != t2.Dimensions())
        throw std::invalid_argument("Different shapes in mult");

    Tucker result(t1._n[0],
                  t1._n[1],
                  t1._n[2],
                  t1._r[0] * t2._r[0],
                  t1._r[1] * t2._r[1],
                  t1._r[2] * t2._r[2]);

    int kA0, kA1, kA2;
    int kB0, kB1, kB2;
    for (int k0 = 0; k0 < t1._r[0] * t2._r[0]; k0++)
    {
        for (int k1 = 0; k1 < t1._r[1] * t2._r[1]; k1++)
        {
            for (int k2 = 0; k2 < t1._r[2] * t2._r[2]; k2++)
            {
                kA0 = k0 / t2._r[0];
                kA1 = k1 / t2._r[1];
                kA2 = k2 / t2._r[2];

                kB0 = k0 % t2._r[0];
                kB1 = k1 % t2._r[1];
                kB2 = k2 % t2._r[2];

                result._core(k0, k1, k2) = t1._core(kA0, kA1, kA2) * t2._core(kB0, kB1, kB2);
            }
        }
    }

    for (int i : {0, 1, 2})
    {
        for (int r = 0; r < t1._n[i]; r++)
        {
            result._u[i].row(r) = kroneckerProduct(t1._u[i].row(r), t2._u[i].row(r));
        }
    }
    return result;
}

Tucker operator*(double d, const Tucker& t)
{
    Tucker result = t;
    for (int i0 = 0; i0 < t._r[0]; i0++)
    {
        for (int i1 = 0; i1 < t._r[1]; i1++)
        {
            for (int i2 = 0; i2 < t._r[2]; i2++)
            {
                result._core(i0, i1, i2) *= d;
            }
        }
    }
    return result;
}

Tucker operator*(const Tucker& t, double d)
{
    return d * t;
}

Tucker operator-(const Tucker& t)
{
    return (-1.0) * t;
}

template<typename Scalar, int rank, typename sizeType>
Map<const Matrix<Scalar, Dynamic, Dynamic>>
TensorToMatrix(const Tensor<Scalar, rank>& tensor,
               const sizeType rows,
               const sizeType cols)
{
    return Map<const Matrix<Scalar, Dynamic, Dynamic>>(tensor.data(), rows, cols);
}

MatrixXd Unfolding(const Tensor3d& tensor,
                   int index)
{
    int I0 = tensor.dimension(0);
    int I1 = tensor.dimension(1);
    int I2 = tensor.dimension(2);

    MatrixXd unfolding;

    switch (index)
    {
    case 0:
        unfolding.resize(I0, I1 * I2);

        for (int j = 0; j < I1; j++)
        {
            for (int i = 0; i < I2; i++)
            {
                Tensor<double, 3> thread = tensor.slice(array<Index, 3>{0, j, i},
                                                        array<Index, 3>{I0, 1, 1});

                unfolding.block(0, i + j * I2, I0, 1) = TensorToMatrix(thread, I0, 1);
            }
        }
        break;
    case 1:
        unfolding.resize(I1, I2 * I0);

        for (int j = 0; j < I2; j++)
        {
            for (int i = 0; i < I0; i++)
            {
                Tensor<double, 3> thread = tensor.slice(array<Index, 3>{i, 0, j},
                                                        array<Index, 3>{1, I1, 1});

                unfolding.block(0, i + j * I0, I1, 1) = TensorToMatrix(thread, I1, 1);
            }
        }
        break;
    case 2:
        unfolding.resize(I2, I0 * I1);

        for (int j = 0; j < I0; j++)
        {
            for (int i = 0; i < I1; i++)
            {
                Tensor3d thread = tensor.slice(array<Index, 3>{j, i, 0},
                                               array<Index, 3>{1, 1, I2});
                                                        
                unfolding.block(0, i + j * I1, I2, 1) = TensorToMatrix(thread, I2, 1);
            }
        }
        break;
    }
    return unfolding;
}

Tensor3d Folding(int I0, int I1, int I2,
                 const MatrixXd& unfolding,
                 int index)
{
    Tensor3d folding(I0, I1, I2);

    switch (index)
    {
    case 0:
        for (int i0 = 0; i0 < I0; i0++)
        {
            for (int i1 = 0; i1 < I1; i1++)
            {
                for (int i2 = 0; i2 < I2; i2++)
                {
                    folding(i0, i1, i2) = unfolding(i0, i2 + i1 * I2);
                }
            }
        }
        break;
    case 1:
        for (int i1 = 0; i1 < I1; i1++)
        {
            for (int i2 = 0; i2 < I2; i2++)
            {
                for (int i0 = 0; i0 < I0; i0++)
                {
                    folding(i0, i1, i2) = unfolding(i1, i0 + i2 * I0);
                }
            }
        }
        break;
    case 2:
        for (int i2 = 0; i2 < I2; i2++)
        {
            for (int i0 = 0; i0 < I0; i0++)
            {
                for (int i1 = 0; i1 < I1; i1++)
                {
                    folding(i0, i1, i2) = unfolding(i2, i1 + i0 * I1);
                }
            }
        }
        break;
    }
    return folding;
}

void Tucker::_ComputeU(const Tensor3d& tensor, double precision, int rmax)
{
    for (int i : {0, 1, 2})
    {
        BDCSVD<MatrixXd> SVD(Unfolding(tensor, i), ComputeThinU | ComputeThinV);
        VectorXd SV = SVD.singularValues();
        _u[i] = SVD.matrixU();

        double threshold = precision * SV.norm() / sqrt(3);

        vector<int> colsToKeep;
        int r = 0;
        for (int j = 0; j < _u[i].cols(); j++)
        {
            if ((r == 0) || (SV(j) > threshold && r < rmax))
            {
                colsToKeep.push_back(j);
                r++;
            }
        }
        VectorXi rowsToKeep = VectorXi::LinSpaced(_u[i].rows(), 0, _u[i].rows());
        _u[i] = _u[i](rowsToKeep, colsToKeep).eval();
    }
}
}