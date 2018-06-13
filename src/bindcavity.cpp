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
    .def("inf_gf", &Cavity::inf_gf)
    .def("d_inf_gf", &Cavity::d_inf_gf)
    .def("omega", &Cavity::omega);
  m.attr("GPAR") = GPAR;
  m.attr("ALPHA") = ALPHA;
  m.attr("C") = C;
}
