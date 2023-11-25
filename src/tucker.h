#pragma once

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

template<typename Scalar, int rank, typename sizeType>
Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>
TensorToMatrix(const Eigen::Tensor<Scalar, rank>& tensor,
				const sizeType rows,
				const sizeType cols);

Eigen::MatrixXd Unfolding(const Eigen::Tensor<double, 3>& tensor, int index);

Eigen::Tensor<double, 3> Folding(int I0, int I1, int I2,
								 const Eigen::MatrixXd& unfolding,
								 int index);

class Tucker
{
public:
	Tucker();
	// Zero tensor with given ranks
	Tucker(int n0, int n1, int n2,
		   int r0, int r1, int r2);

	// Compress a tensor with a given accuracy
	Tucker(const Eigen::Tensor<double, 3>& tensor,
		   double eps = 1e-14,
		   int rmax = 1e+6);

	// Create a rank-one tensor from given factors
	Tucker(const std::vector<Eigen::VectorXd>& u);

	int Size() const;

	std::vector<int> Dimensions() const;
	std::vector<int> Ranks() const;
	std::vector<Eigen::MatrixXd> U() const;
	Eigen::Tensor<double, 3> Core() const;

	double At(int i0, int i1, int i2) const;
	Eigen::Tensor<double, 3> Reconstructed() const;

	double Sum() const;
	double Norm() const;

	Tucker Abs() const;

	void Recompress(double eps = 1e-14, int rmax = 1e+6);

	friend std::ostream& operator <<(std::ostream& out, const Tucker& t);

	Tucker& operator +=(const Tucker& t);

	friend Tucker operator +(const Tucker& t1, const Tucker& t2);
	friend Tucker operator -(const Tucker& t1, const Tucker& t2);
	friend Tucker operator *(const Tucker& t1, const Tucker& t2);
	friend Tucker operator *(const double alpha, const Tucker& t);
	friend Tucker operator *(const Tucker& t, const double alpha);
	friend Tucker operator /(Tucker t1, const Tucker& t2);
	friend Tucker operator -(const Tucker& t);

	friend Tucker Reflection(Tucker t, int axis);

private:
	// TODO: Remove u
	void _ComputeU(const Eigen::Tensor<double, 3>& tensor,
				  double eps = 1e-14,
				  int rmax = 1e+6);

private:
	std::vector<int> _n = std::vector<int>(3);
	std::vector<int> _r = std::vector<int>(3);
	std::vector<Eigen::MatrixXd> _u = std::vector<Eigen::MatrixXd>(3);
	Eigen::Tensor<double, 3> _core;
};