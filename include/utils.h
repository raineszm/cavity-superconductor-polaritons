#pragma once

#include <Eigen/Core>
#include <cmath>

using Eigen::Matrix2cd;
using Eigen::Matrix3cd;

inline double
tanh_over(double x, double T)
{
  if (std::fabs(x) > 1e-5 * T) {
    return std::tanh(x / (2 * T)) / x;
  } else {
    return 1 / (2 * T) - x * x / (24 * std::pow(T, 3));
  }
}

inline double
c(double x1, double x2, double T)
{
  if (std::fabs(x1) > 1e-5 * T) {
    return (std::tanh((x1 + x2) / (2 * T)) - std::tanh((-x1 + x2) / (2 * T))) /
           (4 * x1);
  } else {
    return (1. / (T * std::pow(std::cosh(x2 / (2 * T)), 2)) +
            x1 * x1 * (std::cosh(x2 / T) - 2) /
              (24 * std::pow(T, 3) * std::pow(std::cosh(x2 / (2 * T)), 4))) /
           4;
  }
}

inline Matrix2cd
adjugate(const Matrix2cd& m)
{
  return m.trace() * Matrix2cd::Identity() - m;
}

inline Matrix3cd
adjugate(const Matrix3cd& m)
{
  auto m2 = m * m;
  auto trm = m.trace();
  return (trm * trm - m2.trace()) * Matrix3cd::Identity() / 2 - m * trm + m2;
}
