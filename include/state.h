/**
 * @brief Properties of the superconducting state
 *
 * @file state.h
 * @author Zach Raines
 * @date 2018-06-12
 */
#pragma once

#include "roots.h"
#include "system.h"
#include <array>
#include <cmath>

/**
 * @brief Mean field solution
 *
 */
class State
{
public:
  //! Temperature of the State we are describing
  double T;
  //! Associated System object
  const System sys;
  //! Order parameter
  double delta;

  /**
   * @brief Construct a new State object
   *
   * @param sys_ #sys
   * @param T_ #T
   * @param D #delta
   */
  State(const System& sys_, double T_, double D)
    : T(T_)
    , sys(sys_)
    , delta(D)
  {}

  bool operator==(const State& rhs) const
  {
    return T == rhs.T && sys == rhs.sys && delta == rhs.delta;
  }

  bool operator!=(const State& rhs) const { return !((*this) == rhs); }

  /**
   * @brief Solve for the gap given a system and a temperature
   *
   * @param sys the system being studied
   * @param T temperature
   * @return State
   *
   * This includes the effects of the superfluid velocity \f$v_s\f$
   */
  inline static State solve(const System& sys, double T)
  {
    if (T >= sys.Tc) {
      return State(sys, T, 0.);
    }
    auto f = [sys, T](double x) { return sys.gap_eq(T, x); };

    auto gsl_f = make_gsl_function(f);
    auto solver =
      FSolver::create(gsl_root_fsolver_brent, gsl_f, 0., 5 * sys.Tc);

    return State(sys, T, solver.solve(1e-4 * sys.Tc, 0));
  }

  /**
   * @name Coherence functions

    These methods evaluate the products of the coherence functions
    \f[
    l = u u' + v v'\\
    m = u v' + v u'\\
    n = u u' - v v'\\
    p = u v' - v u'
    \f]

    @param x1 \f$\xi\f$
    @param x2 \f$\xi'\f$

    @{
  */
  double l2(double x1, double x2) const
  {
    auto l1 = std::hypot(x1, delta);
    auto l2 = std::hypot(x2, delta);

    return 0.5 * (1 + (x1 * x2 + delta * delta) / (l1 * l2));
  }

  double m2(double x1, double x2) const
  {
    auto l1 = std::hypot(x1, delta);
    auto l2 = std::hypot(x2, delta);

    return 0.5 * (1 + (-x1 * x2 + delta * delta) / (l1 * l2));
  }

  double n2(double x1, double x2) const
  {
    auto l1 = std::hypot(x1, delta);
    auto l2 = std::hypot(x2, delta);

    return 0.5 * (1 + (x1 * x2 - delta * delta) / (l1 * l2));
  }

  double p2(double x1, double x2) const
  {
    auto l1 = std::hypot(x1, delta);
    auto l2 = std::hypot(x2, delta);

    return 0.5 * (1 - (x1 * x2 + delta * delta) / (l1 * l2));
  }

  double ln(double x1, double x2) const
  {
    auto l1 = std::hypot(x1, delta);
    auto l2 = std::hypot(x2, delta);

    return 0.5 * (x1 / l1 + x2 / l2);
  }

  double mp(double x1, double x2) const
  {
    auto l1 = std::hypot(x1, delta);
    auto l2 = std::hypot(x2, delta);

    return 0.5 * (x1 / l1 - x2 / l2);
  }

  //! @}

  using pickle_type = std::array<double, 2>;
  pickle_type pickle() const { return { { T, delta } }; }

  static inline State unpickle(const System& sys, const pickle_type& a)
  {
    return State(sys, a[0], a[1]);
  }

  /**
   * @brief The polarization bubble
   *
   * @param E1 quasiparticle energy
   * @param E2 quasiparticle energy
   * @param omega photon frequency
   * @param deriv is this the bubble or its derivative w.r.t Freq
   * @return double
   * @see Coupling::photon_se_int(), Coupling::photon_se()
   *
   * \f[\pi_0(E_1, E_2, \omega) =
   * \frac{n_F(E_2) - n_F(E_1)}{i\Omega_m - E_1 + E_2}\f]
   *
   * For numerical convenience we make use of the relation
   * \f$n_f = \tfrac{1 - \tanh}{2}\f$ to rewrite this as
   * \f[\pi_0(E_1, E_2, \omega) =
   * \frac{1}{2}\frac{\tanh\frac{E_1}{2T}- \tanh\frac{E_2}{2T}}{i\Omega_m-
   * E_1 + E_2}\f]
   *
   * @note This is analytically continued to real frequency
   */
  double pi0(double E1, double E2, double omega, bool deriv) const
  {
    if (deriv) {
      return -0.5 * (std::tanh(E1 / (2 * T)) - std::tanh(E2 / (2 * T))) /
             std::pow(omega - E1 + E2, 2);

    } else {
      return 0.5 * (std::tanh(E1 / (2 * T)) - std::tanh(E2 / (2 * T))) /
             (omega - E1 + E2);
    }
  }
};

namespace std {
template<>
struct hash<State>
{
  typedef State argument_type;
  typedef std::size_t result_type;
  result_type operator()(argument_type const& s) const noexcept
  {
    result_type const hT = std::hash<double>{}(s.T);
    result_type const hsys = std::hash<System>{}(s.sys);
    result_type const hdelta = std::hash<double>{}(s.delta);
    return hT ^ (hsys << 1) ^ (hdelta << 2);
  }
};
}