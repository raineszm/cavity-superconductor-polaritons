#pragma once

#include "roots.h"
#include "system.h"
#include <cmath>

//! Mean Field solution
class State
{
public:
  //! Temperature of the State we are describing
  double T;
  //! Associated System object
  const System sys;
  //! Order parameter
  double delta;

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

  //! Solve the mean field problem for System sys at temperature \f$T\f$

  //! This includes the effects of the superfluid velocity \f$v_s\f$
  inline static State solve(const System& sys, double T)
  {
    if (T >= sys.Tc) {
      return State(sys, T, 0.);
    }
    auto f = [sys, T](double x) { return sys.gap_eq(T, x); };

    auto gsl_f = gsl_function_pp(f);
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