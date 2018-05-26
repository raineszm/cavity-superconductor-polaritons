#pragma once

#include <gsl/gsl_math.h>
#include <tuple>

// Build gsl_function from lambda
template<typename F>
class gsl_function_pp : public gsl_function
{
  const F func;
  static double invoke(double x, void* params)
  {
    return static_cast<gsl_function_pp*>(params)->func(x);
  }

public:
  gsl_function_pp(const F& f)
    : func(f)
  {
    function = &gsl_function_pp::invoke; // inherited from gsl_function
    params = this;                       // inherited from gsl_function
  }
  operator gsl_function*() { return this; }
};

// Helper function for template construction
template<typename F>
gsl_function_pp<F>
make_gsl_function(const F& func)
{
  return gsl_function_pp<F>(func);
}

template<typename F, typename G>
class gsl_function_fdf_pp : public gsl_function_fdf
{
private:
  const F& _f;
  const G& _g;

public:
  gsl_function_fdf_pp(const F& ff, const G& gg)
    : _f(ff)
    , _g(gg)
  {
    f = &gsl_function_fdf_pp::invoke_f;
    df = &gsl_function_fdf_pp::invoke_df;
    fdf = &gsl_function_fdf_pp::invoke_fdf;
    params = this;
  }

  static double invoke_f(double x, void* params)
  {
    return static_cast<gsl_function_fdf_pp*>(params)->_f(x);
  }

  static double invoke_df(double x, void* params)
  {
    return std::get<1>(static_cast<gsl_function_fdf_pp*>(params)->_g(x));
  }

  static void invoke_fdf(double x, void* params, double* f, double* df)
  {
    std::tie(*f, *df) = static_cast<gsl_function_fdf_pp*>(params)->_g(x);
  }
};
