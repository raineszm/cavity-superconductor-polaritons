#include "meanfield.h"

double
MeanField::Tc(double g, const System& sys)
{
  nlopt::opt opt(nlopt::LN_BOBYQA, 1);
  opt.set_lower_bounds(1e-6);
  opt.set_upper_bounds(sys.mu);

  std::function<double(const std::vector<double>&)> objective =
    [sys, g](const std::vector<double>& T) {
      auto mf = MeanField(0., g, T[0], sys);
      return std::pow(mf.delta_ratio() - 1,2);
    };

  opt.set_min_objective(&wrap, static_cast<void*>(&objective));
  opt.set_xtol_abs(1e-6);
  opt.set_maxeval(200);

  std::vector<double> T = { sys.mu * std::exp(-1 / g) };
  double fval;

  opt.optimize(T, fval);

  return T[0];
}
