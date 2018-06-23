/**
 * @brief the physics of the photons
 *
 * @file cavity.h
 * @author Zach Raines
 * @date 2018-06-12
 *
 * We use Gaussian natural units
 */
#pragma once
#include "utils.h"
#include <Eigen/Core>
#include <cmath>

using Eigen::Matrix2d;

//! The fine structure constant \f$\alpha\f$
const double ALPHA = 7.2973525664e-3;

//! The speed of light
const double C = 1;

//! The paramagnetic coupling strength
const double GPAR = std::sqrt(ALPHA / C);

/**
 * @brief The physics of the cavity modes
 *
 */
class Cavity
{
public:
  /**
   * @brief The photon 'mass'
   *
   */
  double omega0;

  explicit Cavity(double omega0_)
    : omega0(omega0_)
  {}

  /**
   * @brief The cavity dispersion
   *
   * @param q momentum
   * @return double
   *
   * The bare photon cavity disperions \f$\omega_q^2 = \omega_0^2 + c^2q^2\f$.
   */
  double omega(double q) const
  {
    return std::sqrt(omega0 * omega0 + C * C * q * q);
  }

  /**
   * @brief The size of the cavity
   *
   * We determine this by relating it to the photon 'mass' \f$\omega_0\f$
   * as \f$L=\tfrac{\pi c}{\omega_0}\f$.
   *
   * @return double
   */
  double L() const { return M_PI * C / omega0; }

  /**
   * @brief The matrix structure of the inverse photon Green's function in the
   * \f$\parallel-\perp\f$ field basis.
   *
   * @param q photon momentum
   * @param theta_q angle of the momentum w.r.t to \f$v_s\f$
   * axis
   * @return Matrix2d
   * @sa inv_gf()
   *
   * \f[
   * \left(1 +\frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
   * - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin
   * 2\theta_q \sigma_1 - \cos 2\theta_q\sigma_3\right)
   * \f]
   */
  Matrix2d matrix_structure(double q, double theta_q) const
  {
    auto wq = this->omega(q);

    auto P0 = 1 + wq * wq / (omega0 * omega0);
    auto P1 = C * C * q * q / (omega0 * omega0) * std::sin(2 * theta_q);
    auto P3 = -C * C * q * q / (omega0 * omega0) * std::cos(2 * theta_q);

    return (Matrix2d() << P0 + P3, P1, P1, P0 - P3).finished();
  }

  /**
   * @brief The inverse Green's function of the photon fields in the
   * \f$\parallel-\perp\f$ basis
   * @verbatim embed:rst:leading-asterisk
   * As derived in :ref:`photon-action`
   * @endverbatim
   *
   * @param omega the frequency
   * @param q photon momentum
   * @param theta_q angle of the momentum w.r.t \f$v_s\f$
   * @return Matrix2d
   * @see matrix_structure()
   *
   * We take the inverse Green's function from the above action and analytically
   * continue it to real frequency
   *
   *
   * \f[
   * D^{-1}(\omega, \mathbf{q}) = \frac{L}{32 \pi c^2}\left[
   * (\omega + i0^+
   * )^2 - \omega_\mathbf{q}^2\right] \left[ \left(1 +
   * \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
   * - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin
   * 2\theta_q\sigma_1 - \cos 2\theta_q\sigma_3\right)
   * \right]
   * \f]
   */
  Matrix2d inv_gf(double omega, double q, double theta_q) const
  {
    auto omega_q = this->omega(q);
    auto prefactor =
      L() * (omega * omega - omega_q * omega_q) / (32 * M_PI * C * C);

    return prefactor * matrix_structure(q, theta_q);
  }

  /**
   * @brief The derivative of the cavity inverse photon green's function w.r.t
   * frequency
   *
   * @param omega the frequency
   * @param q photon momentum
   * @param theta_q angle of the momentum w.r.t \f$v_s\f$
   * @return Matrix2d
   * @see inv_gf()
   *
   * \f[ \frac{\partial D^{-1}(\omega, \mathbf q)}{d\omega}\f]
   */
  Matrix2d d_inv_gf(double omega, double q, double theta_q) const
  {
    auto prefactor = L() * omega / (16 * M_PI * C * C);

    return prefactor * matrix_structure(q, theta_q);
  }

  /** The polarization vectors of the problem at \f$z=L/2\f$ in the basis where
   * \f$\theta_s\f$ defines the \f$x\f$-axis.
   *
   * The vectors are along the columns
   * The factor \f$i\sqrt{\tfrac{2}{L}}\f$ is not included.
   *
   * \note this is not a unitary transformation.
   */
  Matrix2d polarizations(double q, double theta_q) const
  {
    double e2_factor = omega0 / omega(q);

    auto [q1, q2] = gsl::polar_to_rect(1, theta_q);

    return (Matrix2d() << -q2, e2_factor * q1, -q1, e2_factor * -q2).finished();
  }

  using pickle_type = double;

  pickle_type pickle() const { return omega0; }
  static inline Cavity unpickle(pickle_type omega0) { return Cavity(omega0); }
};
