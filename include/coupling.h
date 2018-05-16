#pragma once
#include <cmath>

#include <rzmcmt/fermi.h>

#include "cavity.h"
#include "state.h"

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
  const State state;

  explicit Coupling(const State& state_)
    : state(state_)
  {}

  double ImDA_int(double l, double theta, double omega) const
  {
    double fd = std::sqrt(2) * std::cos(2 * theta);
    double drift = state.sys.drift_theta(state.sys.kf(), theta);
    double Ep = drift + l;
    double Em = drift - l;

    return fd *
           (std::tanh(Ep / (2 * state.T)) - std::tanh(Em / (2 * state.T))) /
           ((omega * omega - 4 * l * l) *
            std::sqrt(l * l - state.delta * state.delta));
  }

  double ImDA(double omega) const
  {
    return state.delta * omega / M_PI *
           xi_integrate(
             [this, omega](double l, double theta) {
               return ImDA_int(l, theta, omega);
             },
             state.delta);
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

    auto x1 = state.sys.xi(kpx, kpy);
    auto x2 = state.sys.xi(kmx, kmy);

    auto ki = i == 0 ? kx : ky;
    auto kj = j == 0 ? kx : ky;
    auto vsi = state.sys.vs * angular(state.sys.theta_v, i);
    auto vsj = state.sys.vs * angular(state.sys.theta_v, j);

    auto diagonal =
      (state.l2(x1, x2) * ki * kj + state.n2(x1, x2) * vsi * vsj) /
      (state.sys.m * state.sys.m);
    auto off = (state.p2(x1, x2) * kj * kj + state.m2(x1, x2) * vsi * vsj) /
               (state.sys.m * state.sys.m);

    auto diagonal_mix = state.ln(x1, x2) * (ki * vsj + vsi * kj) / state.sys.m;
    auto off_mix = state.mp(x1, x2) * (ki * vsj + vsi * kj) / state.sys.m;

    auto drift_p = state.sys.drift(kpx, kpy);
    auto drift_m = state.sys.drift(kmx, kmy);

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
    return state.sys.m / M_PI * GPAR * GPAR *
           xi_integrate(
             [this, omega, qx, qy, i, j](double x, double theta) {
               auto k = std::sqrt(2 * state.sys.m * x);
               return photon_se_int(
                 k * std::cos(theta), k * std::sin(theta), omega, qx, qy, i, j);
             },
             0);
  }
};