#pragma once

#include <array>

#include <unsupported/Eigen/CXX11/Tensor>

namespace VlasovTucker
{
typedef Eigen::Tensor<double, 3> Tensor3d;
typedef std::array<double, 3> Vector3d;
}