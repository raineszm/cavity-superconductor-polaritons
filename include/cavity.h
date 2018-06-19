/**
 * @brief the physics of the photons
 *
 * @file cavity.h
 * @author Zach Raines
 * @date 2018-06-12
 *
 * We use Gaussian Hartree Units
 */
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
   * @param qx momentum along \f$v_s\f$
   * @param qy momentum along \f$\hat z \times v_s\f$
   * @return double
   *
   * The bare photon cavity disperions \f$\omega_q^2 = \omega_0^2 + c^2q^2\f$.
   */
  double omega(double qx, double qy) const
  {
    return std::sqrt(omega0 * omega0 + C * C * (qx * qx + qy * qy));
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
   * @param qx momentum along \f$v_s\f$
   * @param qy momentum along \f$\hat z \times v_s\f$
   * @param theta_s the angle that the supercurrent makes with system \f$x\f$
   * axis
   * @return Matrix2d
   * @sa inv_gf()
   *
   * \f[
   * \left(1 +\frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
   * - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin
   * 2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
   * \f]
   */
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

  /**
   * @brief The inverse Green's function of the photon fields in the
   * \f$\parallel-\perp\f$ basis
   * @verbatim embed:rst:leading-asterisk
   * As derived in :ref:`photon-action`
   * @endverbatim
   *
   * @param omega the frequency
   * @param qx momentum along \f$v_s\f$
   * @param qy momentum along \f$\hat z \times v_s\f$
   * @param theta_s the angle that the supercurrent makes with system \f$x\f$
   * @return Matrix2d
   * @see matrix_structure()
   *
   * We take the inverse Green's function from the above action and analytically
   * continue it to real frequency
   *
   *
   * \f[
   * D^{-1}(\omega, \mathbf{q}) = \frac{L \alpha^2}{32 \pi (c\alpha)^3}\left[
   * (\omega + i0^+
   * )^2 - \omega_\mathbf{q}^2\right] \left[ \left(1 +
   * \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right)\sigma_0
   * - \left(1 - \frac{\omega_\mathbf{q}^2}{\omega_0^2}\right) \left(\sin
   * 2(\theta_q - \theta_s)\sigma_1 - \cos 2(\theta_q - \theta_s)\sigma_3\right)
   * \right]
   * \f]
   */
  Matrix2d inv_gf(double omega, double qx, double qy, double theta_s) const
  {
    auto omega_q = this->omega(qx, qy);
    auto prefactor = L() * std::pow(ALPHA, 2) *
                     (omega * omega - omega_q * omega_q) /
                     (32 * M_PI * std::pow(ALPHA * C, 3));

    return prefactor * matrix_structure(qx, qy, theta_s);
  }

  /**
   * @brief The derivative of the cavity inverse photon green's function w.r.t
   * frequency
   *
   * @param omega the frequency
   * @param qx momentum along \f$v_s\f$
   * @param qy momentum along \f$\hat z \times v_s\f$
   * @param theta_s the angle that the supercurrent makes with system \f$x\f$
   * @return Matrix2d
   * @see inv_gf()
   *
   * \f[ \frac{\partial D^{-1}(\omega, \mathbf q)}{d\omega}\f]
   */
  Matrix2d d_inv_gf(double omega, double qx, double qy, double theta_s) const
  {
    auto prefactor =
      L() * std::pow(ALPHA, 2) * omega / (16 * M_PI * std::pow(C * ALPHA, 3));

    return prefactor * matrix_structure(qx, qy, theta_s);
  }

  /** The polarization vectors of the problem at \f$z=L/2\f$ in the basis where
   * \f$\theta_s\f$ defines the \f$x\f$-axis.
   *
   * The vectors are along the columns
   * The factor \f$i\sqrt{\tfrac{2}{L}}\f$ is not includes.
   *
   * \note this is not a unitary transformation.
   */
  Matrix2d polarizations(double qx, double qy, double theta_s) const
  {
    double e2_factor = omega0 / omega(qx, qy);
    double theta_q = std::atan2(qy, qx);

    auto q1 = std::cos(theta_q - theta_s);
    auto q2 = std::sin(theta_q - theta_s);

    return (Matrix2d() << -q2, e2_factor * q1, q1, e2_factor * q2).finished();
  }
};
