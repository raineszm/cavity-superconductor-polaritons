#include <pybind11/pybind11.h>

#include "system.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_system(py::module& m)
{
  py::class_<System>(m, "System")
    .def(py::init<double, double, double, double, double>(),
         "m"_a,
         "mu"_a,
         "Tc"_a,
         "vs"_a,
         "theta_v"_a = 0.)
    // Attributes
    .def_readonly("m", &System::m)
    .def_readonly("mu", &System::mu)
    .def_readonly("Tc", &System::Tc)
    .def_readonly("vs", &System::vs)
    .def_readonly("theta_v", &System::theta_v)
    .def("gap_eq", &System::gap_eq, "the gap equation")
    .def("gap_eq_int",
         &System::gap_eq_int,
         "The integrand appearing in the self-consistency equation")
    .def("xi", &System::xi)
    .def("xi_k", &System::xi_k)
    .def("drift", &System::drift)
    .def("drift_theta", &System::drift_theta)
    .def_property_readonly("kf", &System::kf)
    .def_property_readonly("vf", &System::vf)
    .def_property_readonly("dos", &System::dos)
    // Pickle
    .def(py::pickle(
      [](const System& sys) {
        return py::make_tuple(sys.m, sys.mu, sys.Tc, sys.vs, sys.theta_v);
      },
      [](py::tuple t) {
        return new System(t[0].cast<double>(),
                          t[1].cast<double>(),
                          t[2].cast<double>(),
                          t[3].cast<double>(),
                          t[4].cast<double>());
      }));
}
