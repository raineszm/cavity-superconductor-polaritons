#pragma once

#include <rzmcmt/integrate.h>

template <class F>
class Integrand {
    const F _f;

    double jac(double t) const {
        double t2 = t*t;
        return (1+t2)/((1-t2)*(1-t2));

    }

    double z(double t) const {
        return t/(1-t*t);
    }
    public:

    explicit Integrand(const F& f) : _f(f) {}

    int operator() (double v[1], const double t[2]) const {
      v[0] = _f(z(t[0]), z(t[1]))*jac(t[0])*jac(t[1]);
      return 0;
    }

};
template <class F>
double integrate(const Integrand<F>& i) {
    auto [result, err] = rzmcmt::integrate<2, 1>(i, 0.);
    return result[0]/(4*M_PI*M_PI);
}

template <class F>
double integrate(const F& f) {
    auto i = Integrand(f);
    return integrate(i);
}