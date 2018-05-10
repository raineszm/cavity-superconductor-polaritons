#pragma once
#include <cmath>

#include <rzmcmt/fermi.h>
#include <rzmcmt/integrate.h>

#include "state.h"
#include "system.h"

using rzmcmt::nf;

class Coupling
{
public:
  const System sys;
  const State state;

  double ImDA_int(double kx, double ky, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(std::atan2(ky, ky));
    double x = sys.xi(kx, ky);

    double l = std::hypot(x, state.delta);
    double drift = -sys.As / sys.m * kx; // TODO: fix angle of As

    return fd * (nf(drift - l, state.T) - nf(drift + l, state.T)) /
           ((omega * omega - 4 * l * l) * l);
  }

  double ImDA(double omega) const
  {
    auto integrand = [this, omega](double v[1], const double k[2]) -> int {
      v[0] = ImDA_int(M_PI * k[0], M_PI * k[1], omega);
      return 0;
    };
    auto [result, err] = rzmcmt::integrate<2, 1>(integrand, 0.);

    return 4 * state.delta * result[0] * omega / sys.m;
  }
};