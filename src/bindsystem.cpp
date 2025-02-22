#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
         "theta_s"_a = 0.)
    // Attributes
    .def_readonly("m", &System::m)
    .def_readonly("mu", &System::mu)
    .def_readonly("Tc", &System::Tc)
    .def_readonly("vs", &System::vs)
    .def_readonly("theta_s", &System::theta_s)
    .def("gap_eq", &System::gap_eq, "the gap equation")
    .def("gap_eq_int",
         &System::gap_eq_int,
         "The integrand appearing in the self-consistency equation")
    .def("xi", &System::xi)
    .def("xi_k", &System::xi_k)
    .def("doppler", &System::doppler)
    .def_property_readonly("kf", &System::kf)
    .def_property_readonly("vf", &System::vf)
    .def_property_readonly("dos", &System::dos)
    .def_property_readonly("n", &System::n)
    // Pickle
    .def(
      py::pickle([](const System& sys) { return py::make_tuple(sys.pickle()); },
                 [](py::tuple t) {
                   return System::unpickle(t[0].cast<System::pickle_type>());
                 }));
}
