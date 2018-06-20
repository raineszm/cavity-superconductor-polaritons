#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "coupling.h"
#include "polariton.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_polariton(py::module& m)
{
  py::class_<Polariton> polariton(m, "Polariton");
  polariton
    .def(py::init<const Coupling&, double, double>(),
         "coupling"_a,
         "paraX"_a = 1.,
         "dipoleX"_a = 1.)
    // Attributes
    .def_readonly("coupling", &Polariton::coupling)
    .def_readonly("paraX", &Polariton::paraX)
    .def_readonly("dipoleX", &Polariton::dipoleX)
    .def_property_readonly("sys", &Polariton::sys)
    .def_property_readonly("state", &Polariton::state)
    .def_property_readonly("bs", &Polariton::bs)
    .def_property_readonly("cav", &Polariton::cav)
    .def("inv_gf", &Polariton::inv_gf)
    .def("det_inv_gf", &Polariton::det_inv_gf)
    .def("d_inv_gf", &Polariton::d_inv_gf)
    .def("d_det_inv_gf", &Polariton::d_det_inv_gf)
    .def("det_and_d", &Polariton::det_and_d)
    .def("eigval", &Polariton::eigval)
    .def("_extrema",
         &Polariton::_extrema,
         "find extrema",
         "q"_a,
         "theta_q"_a,
         "xl"_a,
         "xu"_a,
         "ftol"_a = 1e-10)
    .def("find_modes",
         &Polariton::find_modes,
         "find the modes",
         "q"_a,
         "theta_q"_a,
         "ftol"_a = 1e-10,
         "double_root_tol"_a = 1e-17);

  py::class_<ModePolariton>(m, "ModePolariton", polariton)
    .def(py::init<const Coupling&, double, double>(),
         "coupling"_a,
         "paraX"_a = 1.,
         "dipoleX"_a = 1.)
    .def("hamiltonian", &ModePolariton::hamiltonian);
}
