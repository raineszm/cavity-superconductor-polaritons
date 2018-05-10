/* #define EIGEN_STACK_ALLOCATION_LIMIT 0 */

/* #include <pybind11/eigen.h> */
#include <pybind11/pybind11.h>

#include "cuba.h"

#include "optimize.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_system(py::module&);
void
bind_meanfield(py::module&);
void
bind_coupling(py::module&);
void
bind_bs(py::module&);
void
bind_polariton(py::module&);

PYBIND11_MODULE(bardasis_schrieffer, m)
{
  // Disable forking in cuba
  // This causes hangs when interfacing with python
  cubacores(0, 0);

  m.doc() = "A short description of the project";

  bind_system(m);
  bind_meanfield(m);
  bind_coupling(m);
  bind_bs(m);
  bind_polariton(m);

  py::enum_<OptResult>(m, "OptResult")
    .value("SUCCESS", OptResult::SUCCESS)
    .value("STOPVAL_REACHED", OptResult::STOPVAL_REACHED)
    .value("FTOL_REACHED", OptResult::FTOL_REACHED)
    .value("XTOL_REACHED", OptResult::XTOL_REACHED)
    .value("MAXEVAL_REACHED", OptResult::MAXEVAL_REACHED)
    .value("MAXTIME_REACHED", OptResult::MAXTIME_REACHED);
}
