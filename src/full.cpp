#include "full.h"

namespace VlasovTucker
{
Full::Full(const Tensor3d& tensor) : _tensor(tensor) {}

int Full::Size() const
{
    return _tensor.size();
}

std::array<int, 3> Full::Dimensions() const
{
    return {(int)_tensor.dimension(0),
            (int)_tensor.dimension(1),
            (int)_tensor.dimension(2)};
}

double Full::operator()(int i0, int i1, int i2) const
{
    return _tensor(i0, i1, i2);
}

Tensor3d Full::Reconstructed() const
{
    return _tensor;
}

double Full::Sum() const {
    return Eigen::Tensor<double, 0>(_tensor.sum())(0);
}

Full& Full::Compress(double precision, int maxRank)
{
    // Do nothing
    return *this;
}

std::ostream& operator<<(std::ostream& out, const Full& t)
{
    out << t._tensor;
    return out;
}

Full& Full::operator+=(const Full& t)
{
    *this = *this + t; 
    return *this;
}

Full& Full::operator-=(const Full& t)
{
    *this = *this - t; 
    return *this;
}

Full& Full::operator*=(const Full& t)
{
    *this = *this * t; 
    return *this;
}

Full& Full::operator*=(double d)
{
    *this = *this * d; 
    return *this;
}

Full operator+(const Full& t1, const Full& t2)
{
    Full result(t1._tensor + t2._tensor);
    return result;
}

Full operator-(const Full& t1, const Full& t2)
{
    Full result(t1._tensor - t2._tensor);
    return result;
}

Full operator*(const Full& t1, const Full& t2)
{
    Full result(t1._tensor * t2._tensor);
    return result;
}

Full operator*(double d, const Full& t)
{
    Full result(d * t._tensor);
    return result;
}

Full operator*(const Full& t, double d)
{
    return d * t;
}

Full operator-(const Full& t)
{
    return (-1) * t;
}
}