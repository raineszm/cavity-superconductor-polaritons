#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "cavity.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_cavity(py::module& m)
{
  py::class_<Cavity>(m, "Cavity")
    .def(py::init<double>(), "omega0"_a)
    // Attributes
    .def_readwrite("omega0", &Cavity::omega0)
    .def("action", &Cavity::action)
    .def("d_action", &Cavity::d_action)
    .def("omega", &Cavity::omega);
  m.attr("GPAR") = GPAR;
  m.attr("ALPHA") = ALPHA;
  m.attr("C") = C;
}
