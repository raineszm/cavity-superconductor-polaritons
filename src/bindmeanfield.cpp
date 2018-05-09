#include <functional>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>

#include "meanfield.h"
#include "state.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_meanfield(py::module& m)
{
  py::class_<MeanField>(m, "MeanField")
    .def(py::init<double, double, const System&>(),
         "g"_a,
         "T"_a,
         "sys"_a)
    // Attributes
    .def_readwrite("delta", &MeanField::delta, "The order parameter")
    .def_readwrite("g", &MeanField::g, "The coupling strength")
    .def_readwrite("T", &MeanField::T, "The coupling temperature")
    // Functions
    .def_static("Tc", &MeanField::Tc, "the critical temperature")
    .def("F", &MeanField::F, "the free energy density")
    .def("xi", &MeanField::xi, "the normal particle energy")
    .def("delta_ratio",
         &MeanField::delta_ratio,
         "The integral appearing in the self-consistency equation")
    .def("solve",
         &MeanField::solve,
         "solve the gap equations for coupling g and temperature T");

  py::class_<State>(m, "State")
    .def_readwrite("T", &State::T)
    .def_readwrite("delta", &State::delta);

  py::class_<Solver>(m, "Solver")
    .def_property_readonly("Tc", &Solver::Tc)
    .def("state", &Solver::state);
  py::class_<MeanFieldSolver, Solver>(m, "MeanFieldSolver")
    .def(py::init<double, const System&>(), "g"_a, "sys"_a);
}
