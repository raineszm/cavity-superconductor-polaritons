#include <functional>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
    .def("d_ImDA", &Coupling::d_ImDA)
    .def("ImDA_int", &Coupling::ImDA_int)
    .def("pi0_elems", &Coupling::pi0_elems)
    .def("photon_se",
         py::overload_cast<double, double, double>(&Coupling::photon_se,
                                                   py::const_),
         "Photon self energy matrix")
    .def("photon_se",
         py::overload_cast<double, double, double, int, int>(
           &Coupling::photon_se, py::const_),
         "Photon self energy matrix components")
    .def("_photon_se", &Coupling::_photon_se)
    .def("d_photon_se", &Coupling::d_photon_se)
    .def("photon_se_int", &Coupling::photon_se_int);
}
