#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <complex>

using boost::math::tools::bracket_and_solve_root, Eigen::Matrix2cd,
  Eigen::Vector2cd;

#include "bs.h"
#include "cavity.h"
#include "coupling.h"

class Polariton
{
public:
  const BS bs;
  const Cavity cav;
  const Coupling coupling;

  Polariton(const BS& bs_, const Cavity& cav_, const Coupling& c_)
    : bs(bs_)
    , cav(cav_)
    , coupling(c_)
  {}

  Matrix2cd action(double omega, double qx, double qy) const
  {
    Matrix2cd mat;
    double c = coupling.ImDA(omega);
    mat << bs.action(omega), std::complex<double>(0, c),
      -std::complex<double>(0, c), cav.action(omega, qx, qy);
    return mat;
  }

  Vector2cd eigen(double omega, double qx, double qy) const
  {
    return action(omega, qx, qy).eigenvalues();
  }

  double find_mode(double qx, double qy) const
  {
    boost::uintmax_t max = 1e5;
    auto [a, b] = bracket_and_solve_root(
      [this, qx, qy](double omega) {
        return std::real(action(omega, qx, qy).determinant());
      },
      bs.root(),
      2.,
      true,
      [this, qx, qy](double a, double b) {
        double x = (a + b) / 2;
        return std::abs(action(x, qx, qy).determinant()) < 1e-6;
      },
      max);
    return (a + b) / 2;
  }
};