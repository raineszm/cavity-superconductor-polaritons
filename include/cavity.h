#pragma once
#include <cmath>

// We use Gaussian Hartree Units

//! The fine structure constant
const double ALPHA = 7.2973525664e-3;

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
    return std::sqrt(omega0 * omega0 + C * C * (qx * qx + qy * qy));
  }

  double action(double omega, double qx, double qy) const
  {
    return omega * omega - this->omega(qx, qy);
  }
};