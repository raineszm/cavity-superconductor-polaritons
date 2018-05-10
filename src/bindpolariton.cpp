#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "bs.h"
#include "cavity.h"
#include "coupling.h"
#include "polariton.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_polariton(py::module& m)
{
  py::class_<Polariton>(m, "Polariton")
    .def(py::init<const BS&, const Cavity&, const Coupling&>(),
         "bs"_a,
         "cav"_a,
         "coupling"_a)
    // Attributes
    .def("action", &Polariton::action)
    .def("find_mode", &Polariton::find_mode);
}
