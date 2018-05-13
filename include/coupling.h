#pragma once
#include <cmath>
#include <functional>

#include <rzmcmt/fermi.h>

#include "integrate.h"
#include "state.h"
#include "system.h"

using rzmcmt::nf;
using namespace std::placeholders;

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
    double drift = sys.drift(kx, ky);

    return fd * (nf(drift - l, state.T) - nf(drift + l, state.T)) /
           ((omega * omega - 4 * l * l) * l);
  }

  double ImDA(double omega) const
  {
    return state.delta * omega / sys.m *
           integrate(std::bind(&Coupling::ImDA_int, this, _1, _2, omega));
  }
};