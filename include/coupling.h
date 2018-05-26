#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>

#include "cavity.h"
#include "state.h"

using Eigen::Matrix2d;

//! The coupling between the Superconductor and Cavity

//! \note All values returned from this class are missing the factor of
//! \f$\sqrt{2/L}\f$ in the effective paramagnetic coupling. This factor is
//! reinstanted in Polariton::action()
class Coupling
{
public:
  //! The associated mean field state

  //! This also holds the related System
  const State state;

  explicit Coupling(const State& state_)
    : state(state_)
  {}

  /** The integrand for the imaginary part of the coupling
   *
   * In evaluating the integral we pull out the constant prefactor
   * \f$-2 \frac{e}{c} \nu v_s \Omega\Delta\f$
   * leaving the integrand
   * \f[
   * \frac{1}{\sqrt{\lambda^2 - \Delta^2}}
   * \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(\Omega + i0^+)^2 -
   * (2\lambda)^2}f_d(\theta) \f]
   *
   * \sa ImDA()
   */
  double ImDA_int(double l, double theta, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(2 * theta);
    double drift = state.sys.drift_theta(state.sys.kf(), theta);
    double Ep = drift + l;
    double Em = drift - l;

    return fd *
           (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
           ((omega * omega - 4 * l * l) *
            std::sqrt(l * l - state.delta * state.delta));
  }

  double d_ImDA_int(double l, double theta, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(2 * theta);
    double drift = state.sys.drift_theta(state.sys.kf(), theta);
    double Ep = drift + l;
    double Em = drift - l;

    return -2 * omega * fd *
           (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
           (std::pow(omega * omega - 4 * l * l, 2) *
            std::sqrt(l * l - state.delta * state.delta));
  }

  /** The coupling between the Bardasis Schrieffer mode and photons
   *
   * This enters the action
   *  \f[ i\sum_q g_q \left(A^\parallel_q \ \bar{d}_{\perp,q} -
   * A^{\parallel\ast}_q d_{\perp,q}\right)\f]
   * with
   * \f[
   *  g_q = -\frac{e}{c} v_s \Omega\Delta \sum_{\mathbf{k}}
   * \frac{n_F(E^-_\mathbf{k})-n_F(E^+_\mathbf{k}) }{(\Omega + i0^+)^2 -
   * (2\lambda_k)^2}\frac{f_d(\mathbf{k})}{\lambda_k} \f]
   * where here we use unrationalized (Gaussian) units
   *
   * We rewrite the integral in the \f$\xi\f$ approximation
   * \f[
   *  g_q \approx -2 \frac{e}{c} \nu v_s \Omega\Delta \int_\Delta^\infty
   * \frac{\lambda d\lambda}{\sqrt{\lambda^2 - \Delta^2}}
   * \int_0^{2\pi}\frac{d\theta}{2\pi}
   * \frac{n_F(E^-(\lambda))-n_F(E^+(\lambda))}{(\Omega + i0^+)^2 -
   * (2\lambda)^2}\frac{f_d(\theta)}{\lambda} \f]
   *
   * \note We take into account the angular factor arising from \f$v_s\cdot A\f$
   * when building the polariton matrix.
   *
   * \sa ImDA_int(), Polariton::action()
   */
  double ImDA(double omega) const
  {
    return -2 * omega * state.delta * state.sys.m * GPAR * state.sys.vs *
           state.sys.dos() *
           gsl_xi_integrate(
             [this, omega](double l, double theta) {
               return ImDA_int(l, theta, omega);
             },
             state.delta);
  }

  double d_ImDA(double omega) const
  {
    return -2 * state.delta * state.sys.m * GPAR * state.sys.vs *
           state.sys.dos() *
           (gsl_xi_integrate(
              [this, omega](double l, double theta) {
                return ImDA_int(l, theta, omega);
              },
              state.delta) +
            omega * gsl_xi_integrate(
                      [this, omega](double l, double theta) {
                        return d_ImDA_int(l, theta, omega);
                      },
                      state.delta));
  }

  /** The polarization bubble entering into the photon self energy
   *
   * \f[\pi_0(E_1, E_2, \omega) =
   * \frac{n_F(E_2) - n_F(E_1)}{\omega + i0 - E_1 + E_2}\f]
   *
   * For numerical convenience we make use of the relation
   * \f$n_f = \tfrac{1 - \tanh}{2}\f$ to rewrite this as
   * \f[\pi_0(E_1, E_2, \omega) =
   * \frac{1}{2}\frac{\tanh\frac{E_1}{2T}- \tanh\frac{E_2}{2T}}{\omega + i0 -
   * E_1 + E_2}\f]
   *
   * \sa photon_se_int(), photon_se()
   */
  double pi0(double E1, double E2, double omega, bool deriv) const
  {
    if (deriv) {
      return -0.5 *
             (std::tanh(E1 / (2 * state.T)) - std::tanh(E2 / (2 * state.T))) /
             std::pow(omega - E1 + E2, 2);

    } else {
      return 0.5 *
             (std::tanh(E1 / (2 * state.T)) - std::tanh(E2 / (2 * state.T))) /
             (omega - E1 + E2);
    }
  }

  /**
   * The traces of pi0() with the pauli matrices.
   * Specifically we return
   * \f[\{\operatorname{tr}\pi_0\tau_0, \operatorname{tr}\pi_0\tau_1,
   * -i\operatorname{tr}\pi_0\tau_2, \operatorname{tr}\pi_0\tau_3\}\f]
   *
   * \sa pi0(), photon_se()
   */
  std::array<double, 4> pi0_elems(double d1,
                                  double d2,
                                  double l1,
                                  double l2,
                                  double omega,
                                  bool deriv) const
  {
    auto p00 = pi0(d1 + l1, d2 + l2, omega, deriv);
    auto p01 = pi0(d1 + l1, d2 - l2, omega, deriv);
    auto p10 = pi0(d1 - l1, d2 + l2, omega, deriv);
    auto p11 = pi0(d1 - l1, d2 - l2, omega, deriv);

    return { { p00 + p11, p01 + p10, p01 - p10, p00 - p11 } };
  }

  double vs_comp(int i) const
  {
    assert(i < 2 and i >= 0);
    if (i == 0) {
      return state.sys.vs;
    } else {
      return 0.;
    }
  }
  double v_comp(double kx, double ky, int i) const
  {
    assert(i < 2 and i >= 0);
    auto v = std::hypot(kx, ky) / state.sys.m;
    auto theta = std::atan2(ky, kx);
    if (i == 0) {
      return v * std::cos(theta - state.sys.theta_v);
    } else {
      return v * std::sin(theta - state.sys.theta_v);
    }
  }

  /** The photon self-energy due to renormalization by the s-wave state
   *
   * This term can be written in the form
   * \f[
   * -\Pi^{ij}(\mathbf{q}, \Omega) =
   * \frac{e^2}{2c^2} \sum_\mathbf{k}
   * \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2,
   * \mathbf{k} - \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k},
   * \mathbf{q}) \f]
   *
   * In this notation \f[
   * T^{ij}_0 =\ell_\mathbf{k,q}^2 v^i v^j + n_\mathbf{k,q}^2\, v_s^i v_s^j\\
   * T^{ij}_1 = p_\mathbf{k,q}^2 v^i v^j + m_\mathbf{k,q}^2\, v_s^i v_s^j\\
   * T^{ij}_2 = -i p_\mathbf{k,q}m_\mathbf{k,q}\left(v_s^i v^j + v^i
   * v_s^j\right)\\
   * T^{ij}_3 = \ell_\mathbf{k,q}n_\mathbf{k,q}\left(v_s^i v^j +
   * v^i v_s^j\right) \f]
   *
   * For covenience we take the the components \f$i=\{0,1\}\f$ and be along and
   * perpendicular to \f$\mathbf{v}_s\f$ respectively. Specifically, if
   * \f$\theta_s\f$ is the angle that \f$\mathbf{v}_s\f$ makes with the
   * \f$x\f$-axis, we define the new photon operators \f[ \begin{pmatrix}
   * A_\parallel\\
   * A_\perp
   * \end{pmatrix}
   * = \hat{R}(\theta_s)
   * \begin{pmatrix}
   * A_x\\
   * A_y
   * \end{pmatrix}
   * \f]
   * where \f$\hat{R}\f$ is a rotation matrix.
   * This induces the transformation \f$\mathbf{v} \to
   * R^{-1}(\theta_s)\mathbf{v}\f$ on the velocities. Concretely then \f[
   * \mathbf{v}_s &\to (v_s, 0)^T\\
   * \mathbf{v} &\to \frac{k}{m} (\cos(\theta - \theta_s),
   * \sin(\theta-\theta_s))^T \f]
   *
   * \sa State::, photon_se(), Polarition::action()
   */
  double photon_se_int(double kx,
                       double ky,
                       double omega,
                       double qx,
                       double qy,
                       int i,
                       int j,
                       bool deriv) const
  {
    auto xp = state.sys.xi(kx + qx / 2, ky + qy / 2);
    auto xm = state.sys.xi(kx - qx / 2, ky - qy / 2);

    auto vi = v_comp(kx, ky, i);
    auto vj = v_comp(kx, ky, j);
    auto vsi = vs_comp(i);
    auto vsj = vs_comp(j);

    auto T0 = state.l2(xp, xm) * vi * vj + state.n2(xp, xm) * vsi * vsj;
    auto T1 = state.p2(xp, xm) * vi * vj + state.m2(xp, xm) * vsi * vsj;
    auto iT2 = state.mp(xp, xm) * (vsi * vj + vi * vsj);
    auto T3 = state.ln(xp, xm) * (vsi * vj + vi * vsj);

    auto lp = std::hypot(xp, state.delta);
    auto lm = std::hypot(xm, state.delta);

    auto dp = state.sys.drift(kx + qx / 2, ky + qy / 2);
    auto dm = state.sys.drift(kx - qx / 2, ky - qy / 2);

    auto pl = pi0_elems(dp, dm, lp, lm, omega, deriv);

    return T0 * pl[0] + T1 * pl[1] + iT2 * pl[2] + T3 * pl[3];
  }

  /** Evaluates the photon self-energy
   *
   * The self-energy is given by
   * \f[
   * -\Pi^{ij}(\mathbf{q}, \Omega) =
   * \frac{e^2}{2c^2} \sum_\mathbf{k}
   * \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2,
   * \mathbf{k} - \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k},
   * \mathbf{q}) \f]
   * as described in photon_se_int().
   *
   * Numerically it is convenient to evaluate this integral in polar coordinates
   * via the change of variables \f$k=\sqrt{2m(\xi + \mu - \frac{1}{2}m
   * v_s^2)}\f$ \f[ -\Pi^{ij}(\mathbf{q}, \Omega) =
   *  \frac{e^2}{2c^2} \nu \int_{-\mu + \frac{1}{2}m
   * v_s^2}^\infty d\xi \int_0^{2\pi}\frac{d\theta}{2\pi} \sum_l
   * \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2, \mathbf{k} -
   * \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k}, \mathbf{q}) \f]
   *
   * We then make use of the fact that \f$\mu\f$ is the largest energy scale in
   * the problem and take the lower limit of the \f$\xi\f$ integral to
   * \f$-\infty\f$ and symmetrize in \f$\xi\f$. Defining \f[ g(\xi, \theta,
   * \mathbf{q}) = \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} +
   * \mathbf{q}/2, \mathbf{k} - \mathbf{q}/2, \Omega)
   * \hat\tau_l\right]T^{ij}_l(\mathbf{k}, \mathbf{q}) \f] i.e. photon_se_int(),
   * we can express the integral as \f[ -\Pi^{ij}(\mathbf{q}, \Omega)\approx
   * \frac{e^2}{2c^2} \nu \int_0^\infty d\xi
   * \int_0^{2\pi}\frac{d\theta}{2\pi} \left[ g(\xi, \theta, \mathbf{q}) +
   * g(-\xi, \theta, \mathbf{q})\right] \f]
   *
   * \sa photon_se_int(), System::xi_k(), System::n()
   *
   * @{
   */
  double photon_se(double omega, double qx, double qy, int i, int j) const
  {
    return _photon_se(omega, qx, qy, i, j, false);
  }

  double _photon_se(double omega,
                    double qx,
                    double qy,
                    int i,
                    int j,
                    bool deriv) const
  {
    auto ret =
      state.sys.dos() *
      gsl_xi_integrate(
        [this, omega, qx, qy, i, j, deriv](double x, double theta) {
          auto k1 =
            std::sqrt(2 * state.sys.m *
                      (x + state.sys.mu -
                       0.5 * state.sys.m * state.sys.vs * state.sys.vs));
          auto k2 =
            std::sqrt(2 * state.sys.m *
                      (-x + state.sys.mu -
                       0.5 * state.sys.m * state.sys.vs * state.sys.vs));
          return photon_se_int(k1 * std::cos(theta),
                               k1 * std::sin(theta),
                               omega,
                               qx,
                               qy,
                               i,
                               j,
                               deriv) +
                 photon_se_int(k2 * std::cos(theta),
                               k2 * std::sin(theta),
                               omega,
                               qx,
                               qy,
                               i,
                               j,
                               deriv);
          ;
        },
        0);

    return 0.5 * GPAR * GPAR * ret;
  }

  Matrix2d photon_se(double omega, double qx, double qy) const
  {
    return (Matrix2d() << photon_se(omega, qx, qy, 0, 0),
            photon_se(omega, qx, qy, 0, 1),
            photon_se(omega, qx, qy, 1, 0),
            photon_se(omega, qx, qy, 1, 1))
      .finished();
  }

  Matrix2d d_photon_se(double omega, double qx, double qy) const
  {
    return (Matrix2d() << _photon_se(omega, qx, qy, 0, 0, true),
            _photon_se(omega, qx, qy, 0, 1, true),
            _photon_se(omega, qx, qy, 1, 0, true),
            _photon_se(omega, qx, qy, 1, 1, true))
      .finished();
  }
  //@}
};