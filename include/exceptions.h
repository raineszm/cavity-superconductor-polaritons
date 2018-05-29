#pragma once

#include <gsl/gsl_errno.h>

#include <sstream>
#include <stdexcept>
#include <string>

namespace gsl {
class GSLException : public std::runtime_error
{
  static std::string format_exception(const char* reason,
                                      const char* file,
                                      int line,
                                      int gsl_errno)
  {
    std::stringstream str;
    str << gsl_strerror(gsl_errno) << std::endl;
    str << file << ":" << line << std::endl;
    str << std::string(reason);
    return str.str();
  }

public:
  GSLException(const char* reason, const char* file, int line, int gsl_errno)
    : std::runtime_error(
        GSLException::format_exception(reason, file, line, gsl_errno))
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
