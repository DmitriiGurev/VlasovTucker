#pragma once

#include "typedefs.h"

#include <unsupported/Eigen/CXX11/Tensor>

namespace VlasovTucker
{
class Full
{
public:
	Full() {}

    Full(const Tensor3d& tensor);

    int Size() const;

	std::array<int, 3> Dimensions() const;

    double operator()(int i0, int i1, int i2) const;

    Tensor3d Reconstructed() const;

    double Sum() const;

    Full& Compress(double precision = 0, int maxRank = 1e+6);

    friend std::ostream& operator<<(std::ostream& out, const Full& t);

    Full& operator+=(const Full& t);
	Full& operator-=(const Full& t);
	Full& operator*=(const Full& t);
	Full& operator*=(double d);

	friend Full operator+(const Full& t1, const Full& t2);
	friend Full operator-(const Full& t1, const Full& t2);
	friend Full operator*(const Full& t1, const Full& t2);
	friend Full operator*(double d, const Full& t);
	friend Full operator*(const Full& t, double d);
	friend Full operator-(const Full& t);

private:
    Eigen::Tensor<double, 3> _tensor;
};
}