#pragma once

#include <cmath>

inline double
tanh_over(double x, double T)
{
  if (std::fabs(x) > 1e-5 * T) {
    return std::tanh(x / (2 * T)) / x;
  } else {
    return 1 / (2 * T) - x * x / (24 * std::pow(T, 3));
  }
}