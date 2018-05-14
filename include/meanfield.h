#pragma once

#include <boost/math/tools/roots.hpp>
#include <cmath>

#include "integrate.h"
#include "system.h"
#include "utils.h"

using boost::math::tools::bracket_and_solve_root;

class State
{
public:
  double T;
  double delta;

  State(double T_, double delta_)
    : T(T_)
    , delta(delta_)
  {}

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

class MeanField
{
public:
  double Tc;
  double T;
  System sys;
  double delta;

  MeanField(double D, double Tc_, double T_, const System& sys_)
    : Tc(Tc_)
    , T(T_)
    , sys(sys_)
    , delta(D)
  {}

  MeanField(double Tc_, double T_, const System& sys_)
    : Tc(Tc_)
    , T(T_)
    , sys(sys_)
    , delta(0.)
  {
    solve();
  }

  void solve()
  {
    if (T >= Tc) {
      delta = 0.;
      return;
    }
    boost::uintmax_t max = 1e5;
    auto [a, b] = bracket_and_solve_root(
      [this](double x) {
        delta = x;
        return gap_eq();
      },
      Tc,
      2.,
      false,
      [](double a, double b) { return b - a < 1e-5; },
      max);

    delta = (a + b) / 2;
  }

  double gap_eq() const
  {
    return integrate(
      [this](double kx, double ky) { return gap_eq_int(kx, ky); });
  }

  double gap_eq_int(double kx, double ky) const
  {

    double x = sys.xi(kx, ky);
    double l = std::hypot(x, delta);
    double drift = sys.drift(kx, ky);
    return c(l, drift, T) - tanh_over(x, Tc) / 2;
  }

  // Coherence factors

  State to_state() const { return State(T, delta); }
};