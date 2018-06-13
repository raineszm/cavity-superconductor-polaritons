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
    .def("inv_gf_int", &BS::inv_gf_int)
    .def("d_inv_gf_int", &BS::d_inv_gf_int)
    .def("inv_gf", &BS::inv_gf)
    .def("d_inv_gf", &BS::d_inv_gf)
    .def_property_readonly("root", &BS::root);
}
