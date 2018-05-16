#pragma once

#include "integrate.h"
#include "utils.h"
#include <cmath>

class System
{
public:
  //! The electron effective mass in units of m_e
  double m;
  //! The chemical potential
  double mu;
  //! Critical temperature
  double Tc;
  //! The superfluid velocity
  double vs;
  //! The angle that the superfluid velocity makes with the x axis
  double theta_v;

  System(double m_, double mu_, double Tc_, double vs_, double theta_v_)
    : m(m_)
    , mu(mu_)
    , Tc(Tc_)
    , vs(vs_)
    , theta_v(theta_v_)
  {}

  double xi(double kx, double ky) const
  {
    return (kx * kx + ky * ky) / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  double xi_k(double k) const
  {
    return k * k / (2 * m) - mu + 0.5 * m * vs * vs;
  }

  double kf() const { return std::sqrt(2 * m * mu); }
  double vf() const { return kf() / m; }

  double drift(double kx, double ky) const
  {
    return drift_theta(std::hypot(kx, ky), std::atan2(ky, kx));
  }

  double drift_theta(double k, double theta) const
  {
    return vs * k * std::cos(theta - theta_v);
  }

  double gap_eq(double T, double delta) const
  {
    return 2 * m / M_PI *
           xi_integrate(
             [this, delta, T](double x, double theta) {
               return gap_eq_int(x, theta, T, delta);
             },
             0.);
  }

  double gap_eq_int(double x, double theta, double T, double delta) const
  {

    double l = std::hypot(x, delta);
    double d = drift_theta(kf(), theta);
    return c(l, d, T) - tanh_over(x, Tc) / 2;
  }
};