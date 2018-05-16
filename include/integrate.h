#pragma once

#include <memory>

#include <gsl/gsl_integration.h>
#include <rzmcmt/integrate.h>

template<class F>
class Integrand
{
  const F _f;

  double jac(double t) const { return 1 / (t * t); }

  double z(double t) const { return (1 - t) / t; }

public:
  explicit Integrand(const F& f)
    : _f(f)
  {}

  int operator()(double v[1], const double t[2]) const
  {
    auto zx = z(t[0]);
    auto zy = z(t[1]);
    v[0] = (_f(zx, zy) + _f(-zx, zy) + _f(zx, -zy) + _f(-zx, -zy)) * jac(t[0]) *
           jac(t[1]);
    return 0;
  }
};

template<class F>
double
integrate(const F& f)
{
  auto i = Integrand(f);
  auto [result, err] = rzmcmt::integrate<2, 1>(i);
  return result[0] / (4 * M_PI * M_PI);
}

using workspace = std::unique_ptr<gsl_integration_workspace,
                                  decltype(&gsl_integration_workspace_free)>;

class IntegrationWorkspace
{
  workspace wsp;

public:
  IntegrationWorkspace(const size_t n = 1000)
    : wsp(gsl_integration_workspace_alloc(n), &gsl_integration_workspace_free)
  {}

  operator gsl_integration_workspace*() { return wsp.get(); }
};

// Build gsl_function from lambda
template<typename F>
class gsl_function_pp : public gsl_function
{
  const F func;
  static double invoke(double x, void* params)
  {
    return static_cast<gsl_function_pp*>(params)->func(x);
  }

public:
  gsl_function_pp(const F& f)
    : func(f)
  {
    function = &gsl_function_pp::invoke; // inherited from gsl_function
    params = this;                       // inherited from gsl_function
  }
  operator gsl_function*() { return this; }
};

// Helper function for template construction
template<typename F>
gsl_function_pp<F>
make_gsl_function(const F& func)
{
  return gsl_function_pp<F>(func);
}

template<typename F>
double
gsl_xi_integrate(const F& f, double a)
{
  double result, abserr, inner_result, inner_abserr;

  const size_t limit = 1000;

  IntegrationWorkspace wsp1(limit);
  IntegrationWorkspace wsp2(limit);

  auto outer = make_gsl_function([&](double x) {
    auto inner = make_gsl_function([&](double theta) { return f(x, theta); });
    gsl_integration_qag(inner,
                        0,
                        2 * M_PI,
                        EPSABS,
                        EPSREL,
                        limit,
                        GSL_INTEG_GAUSS21,
                        wsp1,
                        &inner_result,
                        &inner_abserr);
    return inner_result / (2 * M_PI);
  });
  gsl_integration_qagiu(
    outer, a, EPSABS, EPSREL, limit, wsp2, &result, &abserr);

  return result;
}

template<class F>
class XiIntegrand
{
  const F _f;
  const double _a;

public:
  XiIntegrand(const F& f, double a)
    : _f(f)
    , _a(a)
  {}

  int operator()(double v[1], const double t[2]) const
  {
    auto x = _a + (1 - t[0]) / t[0];
    auto jac = 1 / (t[0] * t[0]);
    auto theta = 2 * M_PI * t[1];
    v[0] = _f(x, theta) * jac;
    return 0;
  }
};

template<class F>
double
xi_integrate(const F& f, double a)
{
  auto i = XiIntegrand(f, a);
  auto [result, err] = rzmcmt::integrate<2, 1>(i);
  return result[0];
}