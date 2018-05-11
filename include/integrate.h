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
integrate(const Integrand<F>& i)
{
  auto [result, err] = rzmcmt::integrate<2, 1>(i);
  return result[0] / (4 * M_PI * M_PI);
}

template<class F>
double
integrate(const F& f)
{
  auto i = Integrand(f);
  return integrate(i);
}