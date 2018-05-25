#pragma once
#include "exceptions.h"
#include "function.h"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <memory>
#include <tuple>

class FSolver
{
  std::unique_ptr<gsl_root_fsolver, decltype(&gsl_root_fsolver_free)> _solver;

public:
  FSolver(const gsl_root_fsolver_type* t)
    : _solver(gsl_root_fsolver_alloc(t), &gsl_root_fsolver_free)
  {}

  template<typename F>
  static FSolver create(const gsl_root_fsolver_type* t,
                        gsl_function_pp<F>& f,
                        double a,
                        double b)
  {
    auto solver = FSolver(t);
    solver.set(f, a, b);
    return solver;
  }

  template<typename F>
  int set(gsl_function_pp<F>& f, double a, double b)
  {
    return gsl_root_fsolver_set(_solver.get(), &f, a, b);
  }

  int step() { return gsl_root_fsolver_iterate(_solver.get()); }

  double root() const { return gsl_root_fsolver_root(_solver.get()); }

  double xl() const { return gsl_root_fsolver_x_lower(_solver.get()); }
  double xu() const { return gsl_root_fsolver_x_upper(_solver.get()); }

  bool found(double epsabs, double epsrel) const
  {
    return gsl_root_test_interval(xl(), xu(), epsabs, epsrel) == GSL_SUCCESS;
  }

  double solve(double epsabs, double epsrel, size_t niter = 100)
  {
    try {
      for (size_t i = 0; i < niter; i++) {
        step();

        if (found(epsabs, epsrel)) {
          return root();
        }
      }
      return NAN;
    } catch (const gsl::RootException&) {
      return NAN;
    }
  }
};

class FDFSolver
{
  std::unique_ptr<gsl_root_fdfsolver, decltype(&gsl_root_fdfsolver_free)>
    _solver;

public:
  FDFSolver(const gsl_root_fdfsolver_type* t)
    : _solver(gsl_root_fdfsolver_alloc(t), &gsl_root_fdfsolver_free)
  {}

  template<typename F, typename U>
  int set(gsl_function_fdf_pp<F, U>& fdf, double guess)
  {
    return gsl_root_fdfsolver_set(_solver.get(), &fdf, guess);
  }

  int step() { return gsl_root_fdfsolver_iterate(_solver.get()); }

  double root() const { return gsl_root_fdfsolver_root(_solver.get()); }
};