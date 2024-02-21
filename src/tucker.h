#pragma once

#include "typedefs.h"

#include <Eigen/Dense>

namespace VlasovTucker
{
class Tucker
{
public:
	Tucker();

	// Zero tensor with given ranks
	Tucker(int n0, int n1, int n2,
		   int r0, int r1, int r2);

	// Compress a tensor with a given accuracy
	Tucker(const Tensor3d& tensor,
		   double precision = 0,
		   int maxRank = 1e+6);

	// Create a tensor with given core tensor and matrices U
	Tucker(const Tensor3d& core,
		   const std::array<Eigen::MatrixXd, 3>& u);

	Tucker& Compress(double precision = 0, int maxRank = 1e+6);

	Tensor3d Reconstructed() const;

	int Size() const;
	std::array<int, 3> Dimensions() const;
	std::array<int, 3> Ranks() const;

	std::array<Eigen::MatrixXd, 3> U() const;
	Tensor3d Core() const;

	// Reductions
	double Sum() const;
	double Norm() const;

	Tucker Abs() const;

public:
	// Element access
	double operator()(int i0, int i1, int i2) const;

	// Element-wise operations
	Tucker& operator+=(const Tucker& t);
	Tucker& operator-=(const Tucker& t);
	Tucker& operator*=(const Tucker& t);
	Tucker& operator*=(double d);

	friend Tucker operator+(const Tucker& t1, const Tucker& t2);
	friend Tucker operator-(const Tucker& t1, const Tucker& t2);
	friend Tucker operator*(const Tucker& t1, const Tucker& t2);
	friend Tucker operator*(double d, const Tucker& t);
	friend Tucker operator*(const Tucker& t, double d);
	friend Tucker operator-(const Tucker& t);

	friend std::ostream& operator<<(std::ostream& out, const Tucker& t);

private:
	void _ComputeU(const Tensor3d& tensor,
				   double precision,
				   int maxRank);

private:
	std::array<int, 3> _n;
	std::array<int, 3> _r;
	std::array<Eigen::MatrixXd, 3> _u;
	Tensor3d _core;
};

// Tensor-matrix conversions
template<typename Scalar, int rank, typename sizeType>
Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
TensorToMatrix(const Eigen::Tensor<Scalar, rank>& tensor,
			   const sizeType rows,
			   const sizeType cols);

Eigen::MatrixXd Unfolding(const Tensor3d& tensor,
                          int index);

Tensor3d Folding(int I0, int I1, int I2,
			     const Eigen::MatrixXd& unfolding,
			     int index);
}