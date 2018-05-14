/* #define EIGEN_STACK_ALLOCATION_LIMIT 0 */

/* #include <pybind11/eigen.h> */
#include <pybind11/pybind11.h>

#include "cuba.h"
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

PYBIND11_MODULE(bardasis_schrieffer, m)
{
  // Disable forking in cuba
  // This causes hangs when interfacing with python
  cubacores(0, 0);

  m.doc() =
    "Calculate the hybridization of photons and Bardasis-Schrieffer modes";

  bind_system(m);
  bind_state(m);
  bind_coupling(m);
  bind_bs(m);
  bind_polariton(m);
  bind_cavity(m);
}
