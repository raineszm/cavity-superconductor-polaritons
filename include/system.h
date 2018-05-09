#pragma once

#include <cmath>

class System
{
public:
  double m;
  double mu;
  double As;

  System(double m_, double mu_, double As_) 
  : m(m_), mu(mu_), As(As_) {}

  double xi(double kx, double ky) const
  {
    return (kx * kx + ky * ky + As * As) / (2 * m) - mu;
  }

  double kf() const { return std::sqrt(2 * m * mu); }
  double vf() const { return kf() / m; }
};