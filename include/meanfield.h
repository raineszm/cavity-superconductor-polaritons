#pragma once

#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <functional>

#include "integrate.h"
#include "system.h"
#include "utils.h"

using namespace std::placeholders;
using boost::math::tools::bracket_and_solve_root;

class MeanField
{
public:
  double Tc;
  double T;
  const System sys;
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
    // solve();
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
    return integrate(std::bind(&MeanField::gap_eq_int, this, _1, _2));
  }

  double gap_eq_int(double kx, double ky) const
  {

    double x = sys.xi(kx, ky);
    double l = std::hypot(x, delta);
    double drift = sys.drift(kx, ky);
    // return (std::tanh((drift + l) / (2 * T)) -
    //         std::tanh((drift - l) / (2 * T))) /
    //          (4 * l) -
    return tanh_over(l, T) / 2 - tanh_over(x, Tc) / 2;
  }
};