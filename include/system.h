#pragma once

#include <cmath>

class System
{
public:
  double m;
  double mu;
  double vs;

  System(double m_, double mu_, double vs_)
    : m(m_)
    , mu(mu_)
    , vs(vs_)
  {}

  double xi(double kx, double ky) const
  {
    return (kx * kx + ky * ky) / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  double kf() const { return std::sqrt(2 * m * mu); }
  double vf() const { return kf() / m; }

  double drift(double kx, double ky) const { return kx * vs; }
};