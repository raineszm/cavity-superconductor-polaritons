#pragma once

#include <memory>

#include <gsl/gsl_integration.h>

const double EPSABS = 1e-8;
const double EPSREL = 1e-5;

// https://stackoverflow.com/a/43636411/267610
class IntegrationWorkspace
{
  std::unique_ptr<gsl_integration_workspace,
                  decltype(&gsl_integration_workspace_free)>
    wsp;

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