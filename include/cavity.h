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
const double GPAR = 1 / C;

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

  double L() const { return M_PI * C / omega0; }

  Matrix2d matrix_structure(double qx, double qy, double theta_s) const
  {
    auto wq = this->omega(qx, qy);
    auto q = std::hypot(qx, qy);
    auto theta_q = q == 0 ? 0. : std::atan2(qy, qx);

    auto P0 = 1 + wq * wq / (omega0 * omega0);
    auto P1 =
      C * C * q * q / (omega0 * omega0) * std::sin(2 * (theta_q - theta_s));
    auto P3 =
      -C * C * q * q / (omega0 * omega0) * std::cos(2 * (theta_q - theta_s));

    return (Matrix2d() << P0 + P3, P1, P1, P0 - P3).finished();
  }

  /** The inverse GF of the photon modes

  \f[
    S_A = \frac{L \alpha^2}{32 \pi (c\alpha)^3}\sum_q \mathbf{A}(-q) \left[ (i
  \omega_m)^2 - \omega_\mathbf{q}^2\right] \left[ \left(1 +
  \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
  - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin
  2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
  \right]
  \mathbf{A}(q)
  \f]
  * \sa matrix_structure()
  */
  Matrix2d action(double omega, double qx, double qy, double theta_s) const
  {
    auto omega_q = this->omega(qx, qy);
    auto prefactor = L() * std::pow(ALPHA, 2) *
                     (omega * omega - omega_q * omega_q) /
                     (32 * M_PI * std::pow(ALPHA * C, 3));

    return prefactor * matrix_structure(qx, qy, theta_s);
  }

  Matrix2d d_action(double omega, double qx, double qy, double theta_s) const
  {
    auto prefactor =
      L() * std::pow(ALPHA, 2) * omega / (16 * M_PI * std::pow(C * ALPHA, 3));

    return prefactor * matrix_structure(qx, qy, theta_s);
  }
};
