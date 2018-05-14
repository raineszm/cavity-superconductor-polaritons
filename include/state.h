#pragma once

#include "system.h"
#include <boost/math/tools/roots.hpp>
#include <cmath>

using boost::math::tools::bracket_and_solve_root;

class State
{
public:
  //! Temperature
  double T;
  //! associated system
  const System sys;
  //! Order parameter
  double delta;

  State(const System& sys_, double T_, double D)
    : T(T_)
    , sys(sys_)
    , delta(D)
  {}

  inline static State solve(const System& sys, double T)
  {
    if (T >= sys.Tc) {
      return State(sys, T, 0.);
    }
    boost::uintmax_t max = 1e5;
    auto [a, b] =
      bracket_and_solve_root([sys, T](double x) { return sys.gap_eq(T, x); },
                             sys.Tc,
                             2.,
                             false,
                             [](double a, double b) { return b - a < 1e-5; },
                             max);

    return State(sys, T, (a + b) / 2);
  }

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
};
