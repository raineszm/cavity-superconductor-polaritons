#pragma once
#include <cmath>

#include <rzmcmt/fermi.h>

#include "cavity.h"
#include "integrate.h"
#include "meanfield.h"
#include "system.h"

using rzmcmt::nf;

inline double
angular(double theta, int i)
{
  if (i == 0) {
    return std::cos(theta);
  } else {
    return std::sin(theta);
  }
}

class Coupling
{
public:
  const System sys;
  const MeanField state;

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
                       double qy,
                       int i,
                       int j) const
  {
    auto kpx = kx + qx / 2;
    auto kpy = ky + qy / 2;
    auto kmx = kx - qx / 2;
    auto kmy = ky - qy / 2;

    auto x1 = sys.xi(kpx, kpy);
    auto x2 = sys.xi(kmx, kmy);

    auto ki = i == 0 ? kx : ky;
    auto kj = j == 0 ? kx : ky;
    auto vsi = sys.vs * angular(sys.theta_v, i);
    auto vsj = sys.vs * angular(sys.theta_v, j);

    auto diagonal =
      (state.l2(x1, x2) * ki * kj + state.n2(x1, x2) * vsi * vsj) /
      (sys.m * sys.m);
    auto off = (state.p2(x1, x2) * kj * kj + state.m2(x1, x2) * vsi * vsj) /
               (sys.m * sys.m);

    auto diagonal_mix = state.ln(x1, x2) * (ki * vsj + vsi * kj) / sys.m;
    auto off_mix = state.mp(x1, x2) * (ki * vsj + vsi * kj) / sys.m;

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

  double photon_se(double omega, double qx, double qy, int i, int j) const
  {
    return 0.5 * GPAR * GPAR *
           integrate([this, omega, qx, qy, i, j](double kx, double ky) {
             return photon_se_int(kx, ky, omega, qx, qy, i, j);
           });
  }
};