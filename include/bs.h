#pragma once
#include "integrate.h"
#include "state.h"
#include "system.h"
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <functional>

using boost::math::tools::bracket_and_solve_root;

//! Describes the Bardasis-Schrieffer mode of the system
class BS
{
public:
  /** The bare Bardasis-Schriefer 'mass'
   *
   * The mass of the BS mode is given by solving
   * \f[
   * \underbrace{\frac{1}{g_d} - \frac{1}{g_s}}_{m_0} - I(\omega)\f]
   * where \f$I(\omega^2)\f$ is the integrand described by
   * action_int()
   */
  double mass;
  //! The associated mean field state
  const State state;

  BS(double mass_, const State& state_)
    : mass(mass_)
    , state(state_)
  {}

  /** The integrand of the function \f$I(\omega^2)\f$ appearing in
   * the BS mode mass eqn.
   *
   * \f[
   * f(\Omega^2) =
   * \Omega^2 \frac{1}{2\lambda_k} \frac{n_F(E^-_\mathbf{k}) -
   * n_F(E^+_\mathbf{k})}{(\Omega + i0^+)^2-(2\lambda_k)^2} + 2\lambda_k
   * \cos(4\theta_k) \frac{n_F(E^-_\mathbf{k}) - n_F(E^+_\mathbf{k})}{(\Omega +
   * i0^+)^2-(2\lambda_k)^2}
   * \f]
   *
   * \note This integral will be calculated within the \f$\xi\f$ approximation
   * as \f[ I(\omega^2) = \int \frac{d^2k}{(2\pi)^2} f(\mathbf{k}, \omega^2)
   * \approx 2\nu\int_\Delta^\infty \frac{\lambda d\lambda}{\sqrt{\lambda^2 -
   * \Delta^2}} \int_0^{2\pi}\frac{d\theta}{2\pi} f(\lambda, \theta,
   * \omega^2)\f] In particular we group the factor of the quasiparticle density
   * of states into the integrand.
   * \verbatim embed:rst:leading-asterisk
   * For more on the :math:`\xi` approx see :ref:`xi-approx`
   * \endverbatim
   *
   * \sa action()
   */
  double action_int(double l, double theta, double omega) const
  {
    double drift = state.sys.drift_theta(state.sys.kf(), theta);

    double Ep = drift + l;
    double Em = drift - l;

    double F1 =
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      (omega * omega - 4 * l * l);
    return F1 * (omega * omega / 2 + 2 * l * l * std::cos(4 * theta)) /
           std::sqrt(l * l - state.delta * state.delta);
  }

  double d_action_int(double l, double theta, double omega) const
  {
    double drift = state.sys.drift_theta(state.sys.kf(), theta);

    double Ep = drift + l;
    double Em = drift - l;

    double F1 =
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      (omega * omega - 4 * l * l);

    double dF1 =
      -(std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      std::pow(omega * omega - 4 * l * l, 2) * 2 * omega;

    return (F1 * omega +
            dF1 * (omega * omega / 2 + 2 * l * l * std::cos(4 * theta))) /
           std::sqrt(l * l - state.delta * state.delta);
  }

  //! The value of the BS Mode inverse GF at frequency \f$\Omega\f$
  double action(double omega) const
  {
    return mass + 2 * state.sys.dos() *
                    gsl_xi_integrate(
                      [this, omega](double l, double theta) {
                        return action_int(l, theta, omega);
                      },
                      state.delta);
  }

  double d_action(double omega) const
  {
    return 2 * state.sys.dos() *
           gsl_xi_integrate(
             [this, omega](double l, double theta) {
               return d_action_int(l, theta, omega);
             },
             state.delta);
  }

  //! The pole of the Bardasis Schrieffer mode GF

  //! Obtained by solving for where action() is zero
  double root() const
  {
    boost::uintmax_t max = 1e7;
    auto [a, b] = bracket_and_solve_root([this](double x) { return action(x); },
                                         mass,
                                         2.,
                                         false,
                                         [this](double a, double b) {
                                           double x = (a + b) / 2;
                                           return std::fabs(action(x)) < 1e-8;
                                         },
                                         max);
    return (a + b) / 2;
  }
};