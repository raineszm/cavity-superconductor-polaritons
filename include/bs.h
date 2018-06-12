/**
 * @brief The Bardasis Schrieffer mode
 *
 * @file bs.h
 * @author Zach Raines
 * @date 2018-06-12
 */
#pragma once
#include "integrate.h"
#include "roots.h"
#include "state.h"
#include "system.h"
#include <cmath>
#include <functional>

/**
 * @brief Describes the Bardasis-Schrieffer mode of the system
 *
 */
class BS
{
private:
  mutable double _mass = NAN;
  mutable size_t _state_hash;
  mutable double _root;

public:
  /**
   * @brief The bare Bardasis-Schrieffer 'mass'
   *
   * The mass of the BS mode is given by solving
   * \f[
   * \underbrace{\frac{1}{g_d} - \frac{1}{g_s}}_{m_0} + I(\Omega^2)\f]
   * where \f$I(\Omega^2)\f$ is the integrand described by
   * gf_int()
   *
   * @note That we define this variable such that \f$m_0 = \nu M\f$
   * where \f$M\f$ is the mass defined here and \f$\nu\f$ is the
   * density of states @ref dos().
   *
   * @seealso gf_int()
   */
  double mass;
  //! The associated mean field state
  const State state;

  /**
   * @brief Construct a new BS object
   *
   * @param mass_
   * @param state_
   */
  BS(double mass_, const State& state_)
    : mass(mass_)
    , state(state_)
  {}

  /**
   * @brief The integrand of the function \f$I(\omega^2)\f$ appearing in
   * the BS mode mass eqn.
   *
   * @param l the BDG quasiparticle energy (without the dopplershift)
   * \f$\lambda\f$
   * @param theta the angle along the Gap minimum
   * @param omega the frequency
   * @return double
   * @seealso gf()
   *
   * This integrand defines the function \f$I((i\Omega_m)^2) = \sum_\mathbf{k}
   * f(i\Omega, \mathbf{k})\f$ appearing in the BS mode action \f[ f((i
   * \Omega_m)^2, \mathbf{k}) = (i\Omega_m)^2 \frac{1}{2\lambda_k}
   * \frac{n_F(E^-_\mathbf{k}) -
   * n_F(E^+_\mathbf{k})}{(i\Omega_m)^2-(2\lambda_k)^2} + 2\lambda_k
   * \cos(4\theta_k) \frac{n_F(E^-_\mathbf{k}) - n_F(E^+_\mathbf{k})}{(i\Omega_m
   * )^2-(2\lambda_k)^2}
   * \f]
   *
   * \note This integral will be calculated within the \f$\xi\f$ approximation
   * as \f[ I((i\Omega_m)^2) = \int \frac{d^2k}{(2\pi)^2} f((i\Omega_m)^2,
   * \mathbf{k}) \approx 2\nu\int_\Delta^\infty \frac{\lambda
   * d\lambda}{\sqrt{\lambda^2 - \Delta^2}} \int_0^{2\pi}\frac{d\theta}{2\pi}
   * f((i\Omega_m)^2, \lambda, \theta,
   * )
   * \f]
   * In particular we group the factor of the quasiparticle density
   * of states into the integrand.
   * \verbatim embed:rst:leading-asterisk
   * For more on the :math:`\xi` approx see :ref:`xi-approx`
   * \endverbatim
   *
   * @note We will be using the integral in the retarded inverse Green's
   * function for the BS mode, so we analytically continue \f$i\Omega_m \to
   * \omega + i0^+\f$.
   */
  double gf_int(double l, double theta, double omega) const
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

  /**
   * @brief The integrand of the derivative of \f$I((i\Omega)^2)\f$ w.r.t
   * frequency
   *
   * @param l the BDG quasiparticle energy (without the dopplershift)
   * \f$\lambda\f$
   * @param theta the angle along the Gap minimum
   * @param omega the frequency
   * @return double
   * @seealso d_gf(), gf_int()
   */
  double d_gf_int(double l, double theta, double omega) const
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

  /**
   * @brief The Bardasis-Schrieffer inverse Green's function
   *
   * @param omega the frequency \f$\omega\f$
   * @return double
   *
   * The BS action
   *
   * \f[
   * S = \frac{1}{\beta}\sum_q \left[
   *     \underbrace{\frac{1}{g_d} - \frac{1}{g_s}}_{m_0} +
   * I((i\Omega_\text{BS})^2) \right] \f]
   *
   * Defines a Green's function
   *
   * \f[
   *  G^{-1}(i\Omega_m) = -I((i\Omega_m)^2) - m_0
   * \f]
   *
   * @note We analytically continue \f$i\Omega_m \to \omega + i0^+\f$
   */
  double gf(double omega) const
  {
    return -state.sys.dos() *
           (mass + 2 * gsl_xi_integrate(
                         [this, omega](double l, double theta) {
                           return gf_int(l, theta, omega);
                         },
                         state.delta));
  }

  //! Derivative of the BS mode inverse GF wrt \f$\Omega\f$
  double d_gf(double omega) const
  {
    return -2 * state.sys.dos() *
           gsl_xi_integrate(
             [this, omega](double l, double theta) {
               return d_gf_int(l, theta, omega);
             },
             state.delta);
  }

  //! The pole of the Bardasis Schrieffer mode GF

  //! Obtained by solving for where gf() is zero
  double root() const
  {
    // If the cached value is still good return it
    if (mass == _mass && std::hash<State>{}(state) == _state_hash) {
      return _root;
    }

    // Keep track of parameters for when we cache new value
    _mass = mass;
    _state_hash = std::hash<State>{}(state);

    auto f = [this](double x) { return gf(x); };

    auto gsl_f = gsl_function_pp<decltype(f)>(f);
    auto solver = FSolver::create<decltype(f)>(
      gsl_root_fsolver_brent, gsl_f, 1e-3 * state.delta, 1.99 * state.delta);

    return (_root = solver.solve(1e-8 * state.delta, 0));
  }

  /**
   * @brief The Hamiltonian of the BS mode operators
   *
   * @return double
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
  double M() const { return 2 * d_gf(bs.root()); }
};
