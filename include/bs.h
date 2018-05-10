#pragma once
#include "state.h"
#include "system.h"
#include <array>
#include <boost/math/tools/roots.hpp>
#include <cmath>
#include <tuple>

using boost::math::tools::bracket_and_solve_root;

class BS
{
public:
  double mass; // 1/g_d - 1/g_s
  const System sys;
  const State state;

  BS(double mass_, const System& sys_, const State& state_)
    : mass(mass_)
    , sys(sys_)
    , state(state_)
  {}

  double action_int(double kx, double ky, double omega) const
  {
    double x = sys.xi(kx, ky);
    double l = std::hypot(x, state.delta);
    double drift = sys.drift(kx, ky);

    double theta = std::atan2(ky, kx);

    double Ep = drift + l;
    double Em = drift - l;

    double F1 =
      (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
      (omega * omega - 4 * l * l);
    return F1 * (omega * omega / (2 * l) + 2 * l * std::cos(4 * theta));
  }

  double action(double omega) const
  {
    auto integrand = [this, omega](double v[1], const double k[2]) -> int {
      v[0] = action_int(M_PI * k[0], M_PI * k[1], omega);
      return 0;
    };
    std::array<double, 1> result, err;
    std::tie(result, err) = rzmcmt::integrate<2, 1>(integrand, 0.);

    return mass - 4 * result[0];
  }

  double root() const
  {
    double a, b;
    boost::uintmax_t max = 1e5;
    std::tie(a, b) =
      bracket_and_solve_root([this](double x) -> double { return action(x); },
                             state.delta,
                             2.,
                             true,
                             [](double a, double b) {
                               double x = (a + b) / 2;
                               return (b - a) <= 1e-6 * x;
                             },
                             max);
    return (a + b) / 2;
  }
};