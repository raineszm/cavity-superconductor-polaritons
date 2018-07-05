#include "exceptions.h"
#include "function.h"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <memory>
#include <tuple>

class FMinimizer
{
  std::unique_ptr<gsl_min_fminimizer, decltype(&gsl_min_fminimizer_free)>
    _minimizer;

public:
  FMinimizer(const gsl_min_fminimizer_type* t)
    : _minimizer(gsl_min_fminimizer_alloc(t), &gsl_min_fminimizer_free)
  {}

  template<typename F>
  static FMinimizer create(const gsl_min_fminimizer_type* t,
                           gsl_function_pp<F>& f,
                           double a,
                           double b)
  {
    auto minimizer = FMinimizer(t);
    minimizer.set(f, a, b);
    return minimizer;
  }

  template<typename F>
  int set(gsl_function_pp<F>& f, double a, double b)
  {
    return gsl_min_fminimizer_set(_minimizer.get(), &f, (a + b) / 2, a, b);
  }

  int step() { return gsl_min_fminimizer_iterate(_minimizer.get()); }

  double min() const { return (xu() + xl()) / 2; }

  double xl() const { return gsl_min_fminimizer_x_lower(_minimizer.get()); }
  double xu() const { return gsl_min_fminimizer_x_upper(_minimizer.get()); }

  bool found(double epsabs, double epsrel) const
  {
    return gsl_min_test_interval(xl(), xu(), epsabs, epsrel) == GSL_SUCCESS;
  }

  double solve(double epsabs, double epsrel, size_t niter = 100)
  {
    try {
      for (size_t i = 0; i < niter; i++) {
        step();

        if (found(epsabs, epsrel)) {
          return min();
        }
      }
      return std::numeric_limits<double>::quiet_NaN();
    } catch (const gsl::RootException&) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }
};
