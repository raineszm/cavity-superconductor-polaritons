#pragma once
#include <cmath>

// We use atomic units

//! The fine structure constant
const double ALPHA = 1 / 137.;

//! The speed of light
const double C = 1 / ALPHA;

//! The paramagnetic coupling strength
const double GPAR = std::sqrt(ALPHA / C);

class Cavity
{
public:
  //! Cavity frequency
  double omega0;

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