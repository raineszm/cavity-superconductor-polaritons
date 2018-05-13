#pragma once
#include <cmath>

#include <rzmcmt/fermi.h>

#include "integrate.h"
#include "meanfield.h"
#include "system.h"

using rzmcmt::nf;

class Coupling
{
public:
  const System sys;
  const State state;

  double ImDA_int(double k, double theta, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(2 * theta);
    auto kx = k * std::cos(theta);
    auto ky = k * std::sin(theta);
    double x = sys.xi(kx, ky);

    double l = std::hypot(x, state.delta);
    double drift = sys.drift(kx, ky);

    return fd * (nf(drift - l, state.T) - nf(drift + l, state.T)) /
           ((omega * omega - 4 * l * l) * l);
  }

  double ImDA(double omega) const
  {
    return state.delta * omega / sys.m *
           angular_integrate([this, omega](double k, double theta) {
             return ImDA_int(k, theta, omega);
           });
  }
};