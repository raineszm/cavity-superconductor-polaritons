#pragma once
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_trig.h>

using Eigen::Matrix2cd;
using Eigen::Matrix3cd;

inline double
tanh_over(double x, double T)
{
  if (std::fabs(x) > 1e-5 * T) {
    return std::tanh(x / (2 * T)) / x;
  } else {
    return 1 / (2 * T) - x * x / (24 * std::pow(T, 3));
  }
}

inline double
c(double x1, double x2, double T)
{
  if (std::fabs(x1) > 1e-5 * T) {
    return (std::tanh((x1 + x2) / (2 * T)) - std::tanh((-x1 + x2) / (2 * T))) /
           (4 * x1);
  } else {
    return (1. / (T * std::pow(std::cosh(x2 / (2 * T)), 2)) +
            x1 * x1 * (std::cosh(x2 / T) - 2) /
              (24 * std::pow(T, 3) * std::pow(std::cosh(x2 / (2 * T)), 4))) /
           4;
  }
}

inline double
diff_of_tanh(double E1, double E2, double T, double cutoff = 1e-4)
{
  if (std::abs(E1) < cutoff * T && std::abs(E2) < cutoff * T) {
    return (E1 - E2) / (2 * T);
  } else if (T < cutoff * std::abs(E1) && T < cutoff * std::abs(E2)) {
    if (E1 * E2 > 0) {
      return 0.;
    } else {
      return std::copysign(2., E1);
    }
  } else {
    return std::tanh(E1 / (2 * T)) - std::tanh(E2 / (2 * T));
  }
}

inline Matrix2cd
adjugate(const Matrix2cd& m)
{
  return m.trace() * Matrix2cd::Identity() - m;
}

inline Matrix3cd
adjugate(const Matrix3cd& m)
{
  auto m2 = m * m;
  auto trm = m.trace();
  return (trm * trm - m2.trace()) * Matrix3cd::Identity() / 2 - m * trm + m2;
}

template<typename F>
inline double
deriv(const F& f, double x, double h)
{
  auto gsl_f = make_gsl_function(f);
  return deriv_gsl(gsl_f, x, h);
}

inline double
deriv_gsl(const gsl_function& f, double x, double h)
{
  double deriv, err;
  gsl_deriv_central(&f, x, h, &deriv, &err);
  return deriv;
}

namespace gsl {

inline std::array<double, 2>
rect_to_polar(double x, double y)
{
  gsl_sf_result r, theta;
  gsl_sf_rect_to_polar(x, y, &r, &theta);
  return { { r.val, theta.val } };
}

inline std::array<double, 2>
polar_to_rect(double r, double theta)
{
  gsl_sf_result x, y;
  gsl_sf_polar_to_rect(r, theta, &x, &y);
  return { { x.val, y.val } };
}
};
