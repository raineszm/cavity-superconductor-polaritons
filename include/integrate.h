#pragma once

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

template<class F>
class AngularIntegrand
{
  const F _f;

public:
  explicit AngularIntegrand(const F& f)
    : _f(f)
  {}

  int operator()(double v[1], const double t[2]) const
  {
    auto x = (1 - t[0]) / t[0];
    auto jac = 1 / (t[0] * t[0]);
    auto theta = 2 * M_PI * t[1];
    v[0] = _f(x, theta) * x * jac;
    return 0;
  }
};

template<class F>
double
angular_integrate(const F& f)
{
  auto i = AngularIntegrand(f);
  auto [result, err] = rzmcmt::integrate<2, 1>(i);
  return result[0] / (2 * M_PI);
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