#pragma once
#include <cmath>

#include "cavity.h"
#include "state.h"

inline double
angular(double theta, int i)
{
  if (i == 0) {
    return std::cos(theta);
  } else {
    return std::sin(theta);
  }
}

//! The coupling between the Superconductor and Cavity
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
  double pi0(double E1, double E2, double omega) const
  {
    return 0.5 *
           (std::tanh(E1 / (2 * state.T)) - std::tanh(E2 / (2 * state.T))) /
           (omega - E1 + E2);
  }

  /** The photon self-energy due to renormalization by the s-wave state
   *
   * This term can be written in the form
   * \f[
   * \Pi^{ij}(\mathbf{q}, \Omega) =
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
   * \sa State::, photon_se()
   */
  double photon_se_int(double kx,
                       double ky,
                       double omega,
                       double qx,
                       double qy,
                       int i,
                       int j) const
  {
    auto kpx = kx + qx / 2;
    auto kpy = ky + qy / 2;
    auto kmx = kx - qx / 2;
    auto kmy = ky - qy / 2;

    auto x1 = state.sys.xi(kpx, kpy);
    auto x2 = state.sys.xi(kmx, kmy);

    auto ki = i == 0 ? kx : ky;
    auto kj = j == 0 ? kx : ky;
    auto vsi = state.sys.vs * angular(state.sys.theta_v, i);
    auto vsj = state.sys.vs * angular(state.sys.theta_v, j);
  }

  /** Evaluates the photon self-energy
   *
   * The self-energy is given by
   * \f[
   * \Pi^{ij}(\mathbf{q}, \Omega) =
   * \frac{e^2}{2c^2} \sum_\mathbf{k}
   * \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2,
   * \mathbf{k} - \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k},
   * \mathbf{q}) \f]
   * as described in photon_se_int().
   *
   * Numerically it is convenient to evaluate this integral in polar coordinates
   * via the change of variables \f$k=\sqrt{2m(\xi + \mu - \frac{1}{2}m
   * v_s^2)}\f$ \f[ \Pi^{ij}(\mathbf{q}, \Omega) = \frac{e^2}{2c^2} \nu
   * \int_{-\mu + \frac{1}{2}m v_s^2}^\infty d\xi
   * \int_0^{2\pi}\frac{d\theta}{2\pi} \sum_l
   * \operatorname{tr}\left[\hat\pi_0(\mathbf{k} + \mathbf{q}/2, \mathbf{k} -
   * \mathbf{q}/2, \Omega) \hat\tau_l\right]T^{ij}_l(\mathbf{k}, \mathbf{q}) \f]
   *
   * We then make use of the fact that \f$\mu\f$ is the largest energy scale in
   * the problem and take the lower limit of the \f$\xi\f$ integral to
   * \f$-\infty\f$ and symmetrize in \f$\xi\f$. Defining \f[ g(\xi, \theta,
   * \mathbf{q}) = \sum_l \operatorname{tr}\left[\hat\pi_0(\mathbf{k} +
   * \mathbf{q}/2, \mathbf{k} - \mathbf{q}/2, \Omega)
   * \hat\tau_l\right]T^{ij}_l(\mathbf{k}, \mathbf{q}) \f] i.e. photon_se_int(),
   * we can express the integral as \f[ \Pi^{ij}(\mathbf{q}, \Omega)\approx
   * \frac{e^2}{2c^2} \nu \int_0^\infty d\xi \int_0^{2\pi}\frac{d\theta}{2\pi}
   * \left[ g(\xi, \theta, \mathbf{q}) + g(-\xi, \theta, \mathbf{q})\right]
   * \f]
   *
   * \sa photon_se_int(), System::xi_k()
   */
  double photon_se(double omega, double qx, double qy, int i, int j) const
  {
    return 0.5 * state.sys.dos() * GPAR * GPAR *
           gsl_xi_integrate(
             [this, omega, qx, qy, i, j](double x, double theta) {
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
                                    j) +
                      photon_se_int(k2 * std::cos(theta),
                                    k2 * std::sin(theta),
                                    omega,
                                    qx,
                                    qy,
                                    i,
                                    j);
               ;
             },
             0);
  }
};