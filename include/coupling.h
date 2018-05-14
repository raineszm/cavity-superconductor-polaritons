#pragma once
#include <cmath>

#include <rzmcmt/fermi.h>

#include "cavity.h"
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

  double pi0(double E1, double E2, double omega) const
  {
    return 0.5 * (nf(E2, state.T) - nf(E1, state.T)) / (omega - E1 + E2);
  }

  double photon_se_int(double kx,
                       double ky,
                       double omega,
                       double qx,
                       double qy) const
  {
    auto kpx = kx + qx / 2;
    auto kpy = ky + qy / 2;
    auto kmx = kx - qx / 2;
    auto kmy = ky - qy / 2;

    auto x1 = sys.xi(kpx, kpy);
    auto x2 = sys.xi(kmx, kmy);

    auto k = std::hypot(kx, ky);
    // auto theta = std::atan2(ky, kx);

    // TODO: add alpha to copulings
    auto diagonal =
      (state.l2(x1, x2) * k * k + state.n2(x1, x2) * sys.vs * sys.vs) /
      (sys.m * sys.m); // INCLUDE ANGULAR FACTORS
    auto off = (state.p2(x1, x2) * k * k + state.m2(x1, x2) * sys.vs * sys.vs) *
               GPAR * GPAR / (sys.m * sys.m); // INCLUDE ANGULAR FACTORS

    auto diagonal_mix =
      state.ln(x1, x2) * k * sys.vs / sys.m; // INCLUDE ANGULAR FACTORS
    auto off_mix =
      state.mp(x1, x2) * k * sys.vs / sys.m; // INCLUDE ANGULAR FACTORS

    auto drift_p = sys.drift(kpx, kpy);
    auto drift_m = sys.drift(kmx, kmy);

    auto l1 = std::hypot(x1, state.delta);
    auto l2 = std::hypot(x2, state.delta);

    return diagonal * (pi0(drift_p + l1, drift_m + l2, omega) +
                       pi0(drift_p - l1, drift_m - l2, omega)) +
           off * (pi0(drift_p + l1, drift_m - l2, omega) +
                  pi0(drift_p - l1, drift_m + l2, omega)) +
           diagonal_mix * (pi0(drift_p + l1, drift_m + l2, omega) -
                           pi0(drift_p - l1, drift_m - l2, omega)) +
           off_mix * (pi0(drift_p + l1, drift_m - l2, omega) -
                      pi0(drift_p - l1, drift_m + l2, omega));
  }

  double photon_se(double omega, double qx, double qy) const
  {
    return 0.5 * GPAR * GPAR *
           integrate([this, omega, qx, qy](double kx, double ky) {
             return photon_se_int(kx, ky, omega, qx, qy);
           });
  }
};