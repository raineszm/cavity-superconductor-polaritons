#include <functional>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include "coupling.h"
#include "state.h"
#include "system.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_coupling(py::module& m)
{
  py::class_<Coupling>(m, "Coupling")
    .def(py::init<const State&>(), "state"_a)
    // Attributes
    .def_readonly("state", &Coupling::state)
    // Functions
    .def("ImDA", &Coupling::ImDA)
    .def("ImDA_int", &Coupling::ImDA_int)
    .def("photon_se", &Coupling::photon_se)
    .def("photon_se_int", &Coupling::photon_se_int);
}
