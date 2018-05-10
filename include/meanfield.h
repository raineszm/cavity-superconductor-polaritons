#pragma once

#include <cmath>
#include <functional>

#include <nlopt.hpp>

#include "optimize.h"
#include "integrate.h"
#include "system.h"

using namespace std::placeholders;

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
  {}

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
    return integrate(std::bind(&MeanField::FMF_int, this, _1, _2));
  }

  double F() const
  {
    double delta_ = delta > 0 ? delta * delta_ratio() : 0.;
    return FMF() + (-delta_ * delta_ + 2 * delta * delta_) / g;
  }

  double delta_ratio() const
  {
    return g*integrate(std::bind(&MeanField::delta_ratio_int, this, _1, _2));
  }

  double xi(double kx, double ky) const { return sys.xi(kx, ky); }

  static double Tc(double g, const System& sys);

private:
  // TODO:
  double FMF_int(double kx, double ky) const
  {
    double x = std::abs(xi(kx, ky));
    double l = std::hypot(x, delta);
    double drift = sys.drift(kx, ky);

    return -(l - drift - x +
             T * (std::log1p(std::exp(-(l + drift) / T)) +
                  std::log1p(std::exp(-(l - drift) / T)) -
                  std::log1p(std::exp(-x / T))));
  }

  double delta_ratio_int(double kx, double ky) const
  {

    double x = xi(kx, ky);
    double l = std::hypot(x, delta);
    double drift = sys.drift(kx, ky);
    return (std::tanh((drift + l) / (2 * T)) -
            std::tanh((drift - l) / (2 * T))) /
           (4 * l);
  }
};