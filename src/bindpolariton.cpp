#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
    .def(py::init<const BS&, const Cavity&, const Coupling&, double>(),
         "bs"_a,
         "cav"_a,
         "coupling"_a,
         "big"_a = 1)
    // Attributes
    .def_readonly("bs", &Polariton::bs)
    .def_readonly("cav", &Polariton::cav)
    .def_readonly("coupling", &Polariton::coupling)
    .def("action", &Polariton::action)
    .def("d_action", &Polariton::d_action)
    .def("det_and_d", &Polariton::det_and_d)
    .def("eigval", &Polariton::eigval)
    .def("find_modes", &Polariton::find_modes);
}
