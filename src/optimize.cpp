#include <functional>

#include "optimize.h"

double
wrap(const std::vector<double>& x, std::vector<double>&, void* data)
{
  return (*static_cast<std::function<double(const std::vector<double>&)>*>(data))(x);
}
