#pragma once

#include "system.h"
#include <boost/math/tools/roots.hpp>
#include <cmath>

using boost::math::tools::bracket_and_solve_root;

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

  //! Solve the mean field problem for System sys at temperature \f$T\f$

  //! This includes the effects of the superfluid velocity \f$v_s\f$
  inline static State solve(const System& sys, double T)
  {
    if (T >= sys.Tc) {
      return State(sys, T, 0.);
    }
    boost::uintmax_t max = 1e5;
    auto [a, b] = bracket_and_solve_root(
      [sys, T](double x) { return sys.gap_eq(T, x); },
      sys.Tc,
      2.,
      false,
      [sys](double a, double b) { return b - a < 1e-4 * sys.Tc; },
      max);

    return State(sys, T, (a + b) / 2);
  }

  /** @name Coherence functions

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
