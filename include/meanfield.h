#pragma once

#include <cmath>
#include <functional>
#include <tuple>

#include <nlopt.hpp>

#include "optimize.h"
#include "rzmcmt/integrate.h"

#include "system.h"

class MeanField
{
public:
  double g;
  double T;
  const System sys;
  double delta;

  MeanField(double D, double g_, double T_, const System& sys_)
    : g(g_)
    , T(T_)
    , sys(sys_)
    , delta(D)
  {
  }

  MeanField(double g_, double T_, const System& sys_)
    : g(g_)
    , T(T_)
    , sys(sys_)
    , delta(0.)
  {
    solve();
  }

  OptResult solve()
  {
    nlopt::opt opt(nlopt::LN_SBPLX, 1);
    opt.set_lower_bounds(0.);
    opt.set_upper_bounds(1.);

    std::function<double(const std::vector<double>&)> objective =
      [this](const std::vector<double>& D) {
        delta = D[0];
        return F();
      };

    opt.set_min_objective(&wrap, static_cast<void*>(&objective));
    opt.set_ftol_abs(1e-6);
    opt.set_maxeval(200);

    std::vector<double> x = { g };
    double fval;

    auto result = opt.optimize(x, fval);

    if (F() > -1e-9) {
      delta = 0.;
    }

    return OptResult(result);
  }

  double FMF() const
  {
    auto integrand = [this](double v[1], const double k[2]) -> int {
      v[0] = FMF_int(M_PI * k[0], M_PI * k[1]);
      return 0;
    };
    std::array<double, 1> result, err;
    std::tie(result, err) = rzmcmt::integrate<2, 1>(integrand, 0.);

    return result[0];
  }

  double F() const
  {
    double delta_ = delta > 0 ? delta * delta_ratio() : 0.;
    return FMF() + (-delta_ * delta_ + 2 * delta * delta_) / g;
  }

  double delta_ratio() const
  {
    auto integrand = [this](double v[1], const double k[2]) -> int {
      v[0] = delta_ratio_int(M_PI * k[0], M_PI * k[1]);
      return 0;
    };
    std::array<double, 1> result, err;
    std::tie(result, err) = rzmcmt::integrate<2, 1>(integrand);

    return g * result[0];
  }

  double ddt_delta_ratio() const
  {
    auto integrand = [this](double v[1], const double k[2]) -> int {
      v[0] = ddt_delta_ratio_int(M_PI * k[0], M_PI * k[1]);
      return 0;
    };
    std::array<double, 1> result, err;
    std::tie(result, err) = rzmcmt::integrate<2, 1>(integrand);

    return g * result[0];
  }

  double xi(double kx, double ky) const { return sys.xi(kx, ky); }

  static double Tc(double g, const System& sys);

private:
  double FMF_int(double kx, double ky) const
  {
    double x = std::abs(xi(kx, ky));
    double E = std::hypot(x, delta);

    return -(E - x +
             2 * T *
               (std::log1p(std::exp(-E / T)) - std::log1p(std::exp(-x / T))));
  }

  double delta_ratio_int(double kx, double ky) const
  {

    double x = xi(kx, ky);
    double E = std::hypot(x, delta);
    return std::tanh(E / (2 * T)) / (2 * E);
  }

  double ddt_delta_ratio_int(double kx, double ky) const
  {

    double x = xi(kx, ky);
    double E = std::hypot(x, delta);
    double cosh = std::cosh(E / (2 * T));
    return 1 / (4 * E * T * cosh * cosh);
  }
};
