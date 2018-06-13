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
    .def("inf_gf_int", &BS::inf_gf_int)
    .def("d_inf_gf_int", &BS::d_inf_gf_int)
    .def("inf_gf", &BS::inf_gf)
    .def("d_inf_gf", &BS::d_inf_gf)
    .def_property_readonly("root", &BS::root);
}
