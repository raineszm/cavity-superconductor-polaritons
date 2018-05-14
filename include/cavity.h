#pragma once
#include <cmath>

const double C = 1.; // speed of light (FIXME)
const double ALPHA = 1 / 137.;
const double GPAR = std::sqrt(ALPHA / C);

class Cavity
{
public:
  double omega0; // Cavity Frequency

  explicit Cavity(double omega0_)
    : omega0(omega0_)
  {}

  double omega(double qx, double qy) const
  {
    return std::sqrt(omega0 * omega0 + C * C * qx * qx + C * C * qy * qy);
  }

  double action(double omega, double qx, double qy) const
  {
    return omega * omega - C * C * (qx * qx + qy * qy) - omega0 * omega0;
  }
};