#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <complex>

using boost::math::tools::bracket_and_solve_root;
using Eigen::Matrix3cd;
using Eigen::Vector3cd;

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

  Matrix3cd action(double omega, double qx, double qy) const
  {
    Matrix3cd mat;
    std::complex<double> c(0., coupling.ImDA(omega));
    auto cs = std::cos(coupling.sys.theta_v);
    auto sn = std::sin(coupling.sys.theta_v);
    auto se00 = coupling.photon_se(omega, qx, qy, 0, 0);
    auto se01 = coupling.photon_se(omega, qx, qy, 0, 1);
    auto se10 = coupling.photon_se(omega, qx, qy, 1, 0);
    auto se11 = coupling.photon_se(omega, qx, qy, 1, 1);
    mat << bs.action(omega), c * cs, c * sn, -c * cs,
      cav.action(omega, qx, qy) + se00, se01, -c * sn, se10,
      cav.action(omega, qx, qy) + se11;
    return mat;
  }

  Vector3cd eigen(double omega, double qx, double qy) const
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