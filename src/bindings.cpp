/* #define EIGEN_STACK_ALLOCATION_LIMIT 0 */

/* #include <pybind11/eigen.h> */
#include <gsl/gsl_errno.h>
#include <pybind11/pybind11.h>

#include <cmath>

namespace py = pybind11;
using namespace py::literals;

void
bind_system(py::module&);
void
bind_state(py::module&);
void
bind_coupling(py::module&);
void
bind_bs(py::module&);
void
bind_polariton(py::module&);
void
bind_cavity(py::module&);

void
handler(const char* reason, const char* file, int line, int gsl_errno)
{
  std::string msg(reason);
  throw std::runtime_error("GSL: " + msg);
}

PYBIND11_MODULE(bardasis_schrieffer, m)
{
  // gsl_set_error_handler_off();
  gsl_set_error_handler(&handler);
  m.doc() =
    "Calculate the hybridization of photons and Bardasis-Schrieffer modes";

  bind_system(m);
  bind_state(m);
  bind_coupling(m);
  bind_bs(m);
  bind_polariton(m);
  bind_cavity(m);
}
