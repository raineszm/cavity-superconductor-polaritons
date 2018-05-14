#include <functional>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include "state.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_state(py::module& m)
{
  py::class_<State>(m, "State")
    .def(py::init<const System&, double, double>(), "sys"_a, "T"_a, "D"_a)
    // Attributes
    .def_readonly("delta", &State::delta, "The order parameter")
    .def_readonly("T", &State::T, "The temperature")
    .def_readonly("sys", &State::sys, "The associated system")
    // Functions
    .def_static(
      "solve", &State::solve, "solve the gap equations at temperature T");
}
