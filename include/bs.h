#pragma once
#include "integrate.h"
#include "roots.h"
#include "state.h"
#include "system.h"
#include <cmath>
#include <functional>

//! Describes the Bardasis-Schrieffer mode of the system
class BS
{
private:
  mutable double _mass = NAN;
  mutable size_t _state_hash;
  mutable double _root;

public:
  /** The bare Bardasis-Schrieffer 'mass'
   *
   * The mass of the BS mode is given by solving
   * \f[
   * \underbrace{\frac{1}{g_d} - \frac{1}{g_s}}_{m_0} + I(\omega^2)\f]
   * where \f$I(\omega^2)\f$ is the integrand described by
   * action_int()
   *
   * In particular however, we divide the gap equation by \f$\nu\f$.
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
           (2 * std::sqrt(l * l - state.delta * state.delta));
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
      -2 * omega *
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      std::pow(omega * omega - 4 * l * l, 2);

    double F2 = (omega * omega / 2 + 2 * l * l * std::cos(4 * theta));

    double dF2 = omega;

    return (F1 * dF2 + dF1 * F2) /
           (2 * std::sqrt(l * l - state.delta * state.delta));
  }

  //! The value of the BS Mode inverse GF at frequency \f$\Omega\f$
  double action(double omega) const
  {
    return -state.sys.dos() *
           (mass + 2 * gsl_xi_integrate(
                         [this, omega](double l, double theta) {
                           return action_int(l, theta, omega);
                         },
                         state.delta));
  }

  //! Derivative of the BS mode inverse GF wrt \f$\Omega\f$
  double d_action(double omega) const
  {
    return -2 * state.sys.dos() *
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
    // If the cached value is still good return it
    if (mass == _mass && std::hash<State>{}(state) == _state_hash) {
      return _root;
    }

    // Keep track of parameters for when we cache new value
    _mass = mass;
    _state_hash = std::hash<State>{}(state);

    auto f = [this](double x) { return action(x); };

    auto gsl_f = gsl_function_pp<decltype(f)>(f);
    auto solver = FSolver::create<decltype(f)>(
      gsl_root_fsolver_brent, gsl_f, 1e-3 * state.delta, 1.99 * state.delta);

    return (_root = solver.solve(1e-8 * state.delta, 0));
  }

  /** The Hamiltonian of the BS mode operators
   *
   * We can define mode operators \f$b,b^\dagger\f$ such the BS action can be
   * rewritten \f[\frac{1}{\beta}\sum_q \bar{b}_q \left(-i \omega_m +
   * \Omega_\text{BS}\right)b_q\f]
   *
   * This defines the Hamiltonian of the BS mode.
   */
  double hamiltonian() const { return root(); }

  /** The kinetic mass of the BS mode
   */
  double M() const { return 2 * state.sys.dos() * mass / (root() * root()); }
};
