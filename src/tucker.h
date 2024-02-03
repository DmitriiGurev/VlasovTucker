#pragma once

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

class Tucker
{
public:
	Tucker();

	// Zero tensor with given ranks
	Tucker(int n0, int n1, int n2,
		   int r0, int r1, int r2);

	// Compress a tensor with a given accuracy
	Tucker(const Eigen::Tensor<double, 3>& tensor,
		   double eps = 0,
		   int rmax = 1e+6);

	// Create a tensor with given core tensor and matrices U
	Tucker(const Eigen::Tensor<double, 3>& core,
		   const std::array<Eigen::MatrixXd, 3>& u);

	int Size() const;

	std::array<int, 3> Dimensions() const;
	std::array<int, 3> Ranks() const;
	std::array<Eigen::MatrixXd, 3> U() const;
	Eigen::Tensor<double, 3> Core() const;

	double At(int i0, int i1, int i2) const;
	Eigen::Tensor<double, 3> Reconstructed() const;

	double Sum() const;
	double Norm() const;

	Tucker Abs() const;

	Tucker& Recompress(double eps = 0, int rmax = 1e+6);

	friend std::ostream& operator<<(std::ostream& out, const Tucker& t);

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

	friend Tucker Reflection(Tucker t, int axis);

private:
	void _ComputeU(const Eigen::Tensor<double, 3>& tensor,
				   double eps = 1e-14,
				   int rmax = 1e+6);

private:
	std::array<int, 3> _n;
	std::array<int, 3> _r;
	std::array<Eigen::MatrixXd, 3> _u;
	Eigen::Tensor<double, 3> _core;
};

// Tensor-matrix conversions
template<typename Scalar, int rank, typename sizeType>
Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
TensorToMatrix(const Eigen::Tensor<Scalar, rank>& tensor,
			   const sizeType rows,
			   const sizeType cols);

Eigen::MatrixXd Unfolding(const Eigen::Tensor<double, 3>& tensor,
                          int index);

Eigen::Tensor<double, 3> Folding(int I0, int I1, int I2,
								 const Eigen::MatrixXd& unfolding,
								 int index);