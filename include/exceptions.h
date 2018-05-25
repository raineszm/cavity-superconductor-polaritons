#pragma once

#include <gsl/gsl_errno.h>

#include <stdexcept>
#include <string>

namespace gsl {
class GSLException : public std::runtime_error
{

public:
  GSLException(const char* reason, const char* file, int line, int gsl_errno)
    : std::runtime_error(gsl_strerror(gsl_errno) + std::string(": ") +
                         std::string(reason))
  {}
};

class RootException : public GSLException
{
public:
  RootException(const char* reason, const char* file, int line, int gsl_errno)
    : GSLException(reason, file, line, gsl_errno)
  {}
};

[[noreturn]] inline void
error_handler(const char* reason, const char* file, int line, int gsl_errno)
{
  if (gsl_errno == GSL_EBADFUNC || gsl_errno == GSL_EZERODIV ||
      gsl_errno == GSL_EINVAL) {
    throw RootException(reason, file, line, gsl_errno);
  } else {
    throw GSLException(reason, file, line, gsl_errno);
  }
}
}
