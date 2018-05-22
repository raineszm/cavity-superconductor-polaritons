//! \file cavity.h
//!  We use Gaussian Hartree Units
#pragma once
#include <Eigen/Core>
#include <cmath>

using Eigen::Matrix2d;

//! The fine structure constant \f$\alpha\f$
const double ALPHA = 7.2973525664e-3;

//! The speed of light
const double C = 1 / ALPHA;

//! The paramagnetic coupling strength
const double GPAR = std::sqrt(ALPHA / C);

//! The physics of the cavity modes
class Cavity
{
public:
  //! Cavity frequency
  double omega0;

  explicit Cavity(double omega0_)
    : omega0(omega0_)
  {}

  //! The dispersion of the cavity
  //! \f[\omega^2 = \omega_0^2 + c^2 q^2\f]
  double omega(double qx, double qy) const
  {
    return std::sqrt(omega0 * omega0 + C * C * (qx * qx + qy * qy));
  }

  /** The inverse GF of the photon modes


  \f[
    S_A = \frac{1}{16 \pi c^2}\sum_q \mathbf{A}(-q) \left[ (i \omega_m)^2 -
  \omega_\mathbf{q}^2\right] \left[ \left(1 +
  \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
  - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin
  2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
  \right]
  \mathbf{A}(q)
  \f]
  **/
  Matrix2d action(double omega, double qx, double qy, double theta_s) const
  {
    auto wq = this->omega(qx, qy);
    auto q = std::hypot(qx, qy);
    auto theta_q = q == 0 ? 0. : std::atan2(qy, qx);
    auto prefactor =
      (omega * omega - this->omega(qx, qy)) / (16 * M_PI * C * C);

    auto P0 = 1 + wq * wq / (omega0 * omega0);
    auto P1 =
      C * C * q * q / (omega0 * omega0) * std::sin(2 * (theta_q - theta_s));
    auto P3 =
      -C * C * q * q / (omega0 * omega0) * std::cos(2 * (theta_q - theta_s));

    return prefactor * (Matrix2d() << P0 + P3, P1, P1, P0 - P3).finished();
  }
};