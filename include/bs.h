#pragma once
#include "integrate.h"
#include "state.h"
#include "system.h"
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <functional>

using boost::math::tools::bracket_and_solve_root;

class BS
{
public:
  double mass; // 1/g_d - 1/g_s
  const State state;

  BS(double mass_, const State& state_)
    : mass(mass_)
    , state(state_)
  {}

  double action_int(double l, double theta, double omega) const
  {
    double drift = state.sys.drift_theta(state.sys.kf(), theta);

    double Ep = drift + l;
    double Em = drift - l;

    double F1 =
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      (omega * omega - 4 * l * l);
    return F1 * (omega * omega / 2 + 2 * l * l * std::cos(4 * theta));
  }

  double action(double omega) const
  {
    return mass + 2 * state.sys.m / M_PI *
                    gsl_xi_integrate(
                      [this, omega](double l, double theta) {
                        return action_int(l, theta, omega) /
                               std::sqrt(l * l - state.delta * state.delta);
                      },
                      state.delta);
  }

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