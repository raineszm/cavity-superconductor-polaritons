#pragma once

#include "function.h"
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

class IntegrationCquadWorkspace
{
  std::unique_ptr<gsl_integration_cquad_workspace,
                  decltype(&gsl_integration_cquad_workspace_free)>
    wsp;

public:
  IntegrationCquadWorkspace(const size_t n = 1000)
    : wsp(gsl_integration_cquad_workspace_alloc(n),
          &gsl_integration_cquad_workspace_free)
  {}

  operator gsl_integration_cquad_workspace*() { return wsp.get(); }
};

class IntegrationQAWOTable
{
  std::unique_ptr<gsl_integration_qawo_table,
                  decltype(&gsl_integration_qawo_table_free)>
    table;

public:
  IntegrationQAWOTable(double omega,
                       double L,
                       enum gsl_integration_qawo_enum sine,
                       size_t n)
    : table(gsl_integration_qawo_table_alloc(omega, L, sine, n),
            &gsl_integration_qawo_table_free)
  {}

  operator gsl_integration_qawo_table*() { return table.get(); }
};

template<typename F>
double
gsl_xi_integrate(const F& f, double a, double epsrel = EPSREL)
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
                        epsrel,
                        limit,
                        GSL_INTEG_GAUSS21,
                        wsp1,
                        &inner_result,
                        &inner_abserr);
    return inner_result / (2 * M_PI);
  });
  gsl_integration_qagiu(
    outer, a, EPSABS, epsrel, limit, wsp2, &result, &abserr);

  return result;
}

template<typename F>
double
gsl_angular_integrate(const F& f, double k1, double k2)
{
  double total = 0;
  double result, abserr, inner_result, inner_abserr;

  const size_t limit = 1000;

  IntegrationWorkspace wsp1(limit);
  IntegrationWorkspace wsp2(limit);
  IntegrationCquadWorkspace cwsp(limit);

  auto outer = make_gsl_function([&](double k) {
    auto inner = make_gsl_function([&](double theta) { return f(k, theta); });
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
    return k * inner_result;
  });

  size_t nevals;
  gsl_integration_cquad(
    outer, k1, k2, EPSABS, EPSREL, cwsp, &result, &abserr, &nevals);
  total += result;
  gsl_integration_qag(outer,
                      0,
                      k1,
                      EPSABS,
                      EPSREL,
                      limit,
                      GSL_INTEG_GAUSS21,
                      wsp2,
                      &result,
                      &abserr);
  total += result;
  gsl_integration_qagiu(
    outer, k2, EPSABS, EPSREL, limit, wsp2, &result, &abserr);
  total += result;

  return total / (4 * M_PI * M_PI);
}

template<typename F>
double
gsl_xi_cos_integrate(const F& f,
                     double a,
                     double theta_s,
                     double epsrel = EPSREL)
{
  double result, abserr, inner_result, inner_abserr;

  const size_t limit = 1000;

  IntegrationWorkspace wsp1(limit);
  IntegrationWorkspace wsp2(limit);
  IntegrationQAWOTable tbl(2., 2 * M_PI, GSL_INTEG_COSINE, limit);

  auto outer = make_gsl_function([&](double x) {
    auto inner =
      make_gsl_function([&](double theta) { return f(x, theta - theta_s); });
    gsl_integration_qawo(
      inner, 0, EPSABS, epsrel, limit, wsp1, tbl, &inner_result, &inner_abserr);
    return inner_result / (2 * M_PI);
  });
  gsl_integration_qagiu(
    outer, a, EPSABS, epsrel, limit, wsp2, &result, &abserr);

  return result;
}