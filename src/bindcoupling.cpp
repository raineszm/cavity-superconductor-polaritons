#include <functional>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include "coupling.h"
#include "meanfield.h"
#include "system.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_coupling(py::module& m)
{
  py::class_<Coupling>(m, "Coupling")
    .def(py::init<const System&, const State&>(), "sys"_a, "state"_a)
    // Attributes
    .def_readonly("sys", &Coupling::sys)
    .def_readonly("state", &Coupling::state)
    // Functions
    .def("ImDA", &Coupling::ImDA)
    .def("photon_se", &Coupling::photon_se);
}
