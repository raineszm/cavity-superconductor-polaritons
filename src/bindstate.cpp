#include <functional>
#include <tuple>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
      "solve", &State::solve, "solve the gap equations at temperature T")
    .def(py::pickle(
      [](const State& s) { return py::make_tuple(s.sys.pickle(), s.pickle()); },
      [](py::tuple t) {
        auto sys = System::unpickle(t[0].cast<System::pickle_type>());
        return State::unpickle(sys, t[1].cast<State::pickle_type>());
      }));
}
