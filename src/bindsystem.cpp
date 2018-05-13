#include <pybind11/pybind11.h>

#include "system.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_system(py::module& m)
{
  py::class_<System>(m, "System")
    .def(py::init<double, double, double>(), "m"_a, "mu"_a, "vs"_a)
    // Attributes
    .def_readwrite("m", &System::m)
    .def_readwrite("mu", &System::mu)
    .def_readwrite("vs", &System::vs)
    .def("xi", &System::xi)
    .def("drift", &System::drift)
    .def_property_readonly("kf", &System::kf)
    .def_property_readonly("vf", &System::vf)
    // Pickle
    .def(py::pickle(
      [](const System& sys) { return py::make_tuple(sys.m, sys.mu, sys.vs); },
      [](py::tuple t) {
        return new System(
          t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());
      }));
}
