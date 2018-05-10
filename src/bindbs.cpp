#include <pybind11/pybind11.h>

#include "bs.h"
#include "state.h"
#include "system.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_bs(py::module& m)
{
  py::class_<BS>(m, "BS")
    .def(py::init<double, const System&, const State&>(),
         "mass"_a,
         "sys"_a,
         "state"_a)
    // Attributes
    .def_readwrite("mass", &BS::mass)
    .def("action_int", &BS::action_int)
    .def("action", &BS::action);
}
