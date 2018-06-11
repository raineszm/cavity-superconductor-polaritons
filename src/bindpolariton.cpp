#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "coupling.h"
#include "polariton.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_polariton(py::module& m)
{
  py::class_<Polariton>(m, "Polariton")
    .def(py::init<const Coupling&, double, double>(),
         "coupling"_a,
         "paraX"_a = 1.,
         "dipoleX"_a = 1.)
    // Attributes
    .def_readonly("coupling", &Polariton::coupling)
    .def_property_readonly("sys", &Polariton::sys)
    .def_property_readonly("state", &Polariton::state)
    .def_property_readonly("bs", &Polariton::bs)
    .def_property_readonly("cav", &Polariton::cav)
    .def("action", &Polariton::action)
    .def("det_action", &Polariton::det_action)
    .def("d_action", &Polariton::d_action)
    .def("d_det_action", &Polariton::d_det_action)
    .def("det_and_d", &Polariton::det_and_d)
    .def("eigval", &Polariton::eigval)
    .def("_extrema", &Polariton::_extrema)
    .def("find_modes", &Polariton::find_modes);
}
