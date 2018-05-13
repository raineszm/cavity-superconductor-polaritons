#pragma once

#include <cmath>

class System
{
public:
  double m;
  double mu;
  double vs;
  double theta_v;

  System(double m_, double mu_, double vs_, double theta_v_)
    : m(m_)
    , mu(mu_)
    , vs(vs_)
    , theta_v(theta_v_)
  {}

  double xi(double kx, double ky) const
  {
    return (kx * kx + ky * ky) / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  double kf() const { return std::sqrt(2 * m * mu); }
  double vf() const { return kf() / m; }

  double drift(double kx, double ky) const
  {
    auto theta = std::atan2(ky, kx);
    return vs * std::hypot(kx, ky) * std::cos(theta - theta_v);
  }
};