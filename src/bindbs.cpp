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
    .def("action_int", &BS::action_int)
    .def("d_action_int", &BS::d_action_int)
    .def("action", &BS::action)
    .def("d_action", &BS::d_action)
    .def_property_readonly("root", &BS::root);
}
