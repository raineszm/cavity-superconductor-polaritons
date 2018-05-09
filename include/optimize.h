#pragma once
#include <vector>

double
wrap(const std::vector<double>& x, std::vector<double>&, void* data);

enum class OptResult
{
  SUCCESS = 1,
  STOPVAL_REACHED = 2,
  FTOL_REACHED = 3,
  XTOL_REACHED = 4,
  MAXEVAL_REACHED = 5,
  MAXTIME_REACHED = 6
};
