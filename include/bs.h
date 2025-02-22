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
#include <array>
#include <cmath>
#include <functional>

/**
 * @brief Describes the Bardasis-Schrieffer mode of the system
 *
 */
class BS
{
private:
  mutable double _mass = std::numeric_limits<double>::quiet_NaN();
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
   * inv_gf_int()
   *
   * @note That we define this variable such that \f$m_0 = \nu \tilde{m}_0\f$
   * where \f$\tilde{m}_0\f$ is the mass defined here and \f$\nu\f$ is the
   * density of states @ref System::dos().
   *
   * @see inv_gf_int()
   */
  double mass;
  //! The associated mean field state
  const State state;

  /**
   * @brief Construct a new BS object
   *
   * @param mass_ #mass
   * @param state_ #state
   */
  BS(double mass_, const State& state_)
    : mass(mass_)
    , state(state_)
  {}

  /**
   * @brief The integrand of the function \f$I(\omega^2)\f$ appearing in
   * the BS mode mass eqn.
   *
   * @param l the BDG quasiparticle energy (without the Doppler shift)
   * \f$\lambda\f$
   * @param theta the angle along the Gap minimum
   * @param omega the frequency
   * @return double
   * @see inv_gf()
   *
   * This integrand defines the function \f$I((i\Omega_m)^2) = \sum_\mathbf{k}
   * f((i\Omega)^2, \mathbf{k})\f$ appearing in the BS mode action \f[ f((i
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
  double inv_gf_int(double l, double theta, double omega) const
  {
    double doppler = state.sys.doppler(state.sys.kf(), theta);

    double Ep = doppler + l;
    double Em = doppler - l;

    double F1 =
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      (omega * omega - 4 * l * l);
    return F1 *
           (omega * omega / 2 +
            2 * l * l * std::cos(4 * (theta + state.sys.theta_s))) /
           (2 * std::sqrt(l * l - state.delta * state.delta));
  }

  /**
   * @brief The integrand of the derivative of \f$I((i\Omega)^2)\f$ w.r.t
   * frequency
   *
   * @param l the BDG quasiparticle energy (without the Doppler shift)
   * \f$\lambda\f$
   * @param theta the angle along the Gap minimum
   * @param omega the frequency
   * @return double
   * @see d_inv_gf(), inv_gf_int()
   */
  double d_inv_gf_int(double l, double theta, double omega) const
  {
    double doppler = state.sys.doppler(state.sys.kf(), theta);

    double Ep = doppler + l;
    double Em = doppler - l;

    double F1 =
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      (omega * omega - 4 * l * l);

    double dF1 =
      -2 * omega *
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      std::pow(omega * omega - 4 * l * l, 2);

    double F2 = (omega * omega / 2 +
                 2 * l * l * std::cos(4 * (theta + state.sys.theta_s)));

    double dF2 = omega;

    return (F1 * dF2 + dF1 * F2) /
           (2 * std::sqrt(l * l - state.delta * state.delta));
  }

  double I(double omega) const
  {
    return -state.sys.dos() * 2 *
           gsl_xi_integrate(
             [this, omega](double l, double theta) {
               return inv_gf_int(l, theta, omega);
             },
             state.delta);
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
  double inv_gf(double omega) const
  {
    return -state.sys.dos() * mass + I(omega);
  }

  /**
   * @brief Derivative of the BS mode inverse inv_gf() wrt \f$\Omega\f$
   *
   * @param omega frequency
   * @return double
   */
  double d_inv_gf(double omega) const
  {
    return -2 * state.sys.dos() *
           gsl_xi_integrate(
             [this, omega](double l, double theta) {
               return d_inv_gf_int(l, theta, omega);
             },
             state.delta);
  }

  /**
   * @brief The pole of the Bardasis Schrieffer mode inv_gf
   *
   * @return double \f$\Omega_\text{BS}\f$
   *
   * Obtained by solving for where inv_gf() is zero
   */
  double root() const
  {
    // If the cached value is still good return it
    if (mass == _mass && std::hash<State>{}(state) == _state_hash) {
      return _root;
    }

    // Keep track of parameters for when we cache new value
    _mass = mass;
    _state_hash = std::hash<State>{}(state);

    auto f = [this](double x) { return inv_gf(x); };

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
   *
   * \f[\frac{1}{2} M = \left.\frac{\partial I(\Omega^2)}{\partial
   * \Omega^2}\right|_{\Omega=\Omega_\text{BS}}\f]
   *
   * Note that \f$\tfrac{\partial I(\Omega^2)}{\partial \Omega} = I'(\Omega^2)2
   * \Omega\f$.
   *
   * So \f$M = I'(\Omega_\text{BS}^2)/\Omega_\text{BS}\f$.
   */
  double M() const { return d_inv_gf(root()) / root(); }

  double inv_gf_int_q(double k,
                      double theta_k,
                      double omega,
                      double q,
                      double theta_q,
                      bool deriv) const
  {
    const System& sys = state.sys;
    auto theta_qk = theta_q - theta_k;
    auto kp = std::sqrt(k * k + q * q / 4 + k * q * std::cos(theta_qk));
    auto km = std::sqrt(k * k + q * q / 4 - k * q * std::cos(theta_qk));

    auto xp = sys.xi_k(kp);
    auto xm = sys.xi_k(km);
    auto lp = std::hypot(xp, state.delta);
    auto lm = std::hypot(xm, state.delta);

    auto theta_p = theta_k + std::atan2(q / 2 * std::sin(theta_qk),
                                        k + q / 2 * std::cos(theta_qk));
    auto theta_m = theta_k + std::atan2(-q / 2 * std::sin(theta_qk),
                                        k - q / 2 * std::cos(theta_qk));
    auto dp = sys.doppler(kp, theta_p);
    auto dm = sys.doppler(km, theta_m);

    auto p0 = state.pi0(dp + lp, dm + lm, omega, deriv) +
              state.pi0(dp - lp, dm - lm, omega, deriv);
    auto p1 = state.pi0(dp + lp, dm - lm, omega, deriv) +
              state.pi0(dp - lp, dm + lm, omega, deriv);

    auto lk = std::hypot(sys.xi_k(k), state.delta);
    return (p0 * state.p2(xp, xm) - p1 * state.l2(xp, xm)) *
             std::pow(std::cos(2 * (theta_k + sys.theta_s)), 2) -
           c(lk, sys.doppler(k, theta_k), state.T);
  }

  double Iq(double omega, double q, double theta_q, bool deriv) const
  {
    return gsl_angular_integrate(
      [this, omega, q, theta_q, deriv](double k, double theta) {
        return inv_gf_int_q(k, theta, omega, q, theta_q, deriv);
      },
      state.sys.kf() - 2 * state.delta / state.sys.vf(),
      state.sys.kf() + 2 * state.delta / state.sys.vf());
  }

  double dispersion(double q, double theta_q) const
  {
    auto f = [this, q, theta_q](double x) {
      return -state.sys.dos() * mass + Iq(x, q, theta_q, false);
    };

    auto gsl_f = gsl_function_pp<decltype(f)>(f);
    auto solver = FSolver::create<decltype(f)>(
      gsl_root_fsolver_brent, gsl_f, 1e-3 * state.delta, 1.99 * state.delta);

    return (_root = solver.solve(1e-8 * state.delta, 0));
  }

  using pickle_type = double;

  pickle_type pickle() const { return mass; }
  static inline BS unpickle(const State& state, pickle_type mass)
  {
    return BS(mass, state);
  }
};
