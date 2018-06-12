#include <pybind11/pybind11.h>

#include "bs.h"
#include "state.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_bs(py::module& m)
{
  py::class_<BS>(m, "BS")
    .def(py::init<double, const State&>(), "mass"_a, "state"_a)
    // Attributes
    .def_readonly("state", &BS::state)
    .def_readwrite("mass", &BS::mass)
    .def("gf_int", &BS::gf_int)
    .def("d_gf_int", &BS::d_gf_int)
    .def("gf", &BS::gf)
    .def("d_gf", &BS::d_gf)
    .def_property_readonly("root", &BS::root);
}
