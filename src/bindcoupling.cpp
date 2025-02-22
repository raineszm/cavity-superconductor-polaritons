#include <functional>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cavity.h"
#include "coupling.h"
#include "state.h"
#include "system.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_coupling(py::module& m)
{
  py::class_<Coupling>(m, "Coupling")
    .def(py::init<const BS&, const Cavity&>(), "bs"_a, "cav"_a)
    // Attributes
    .def_readonly("bs", &Coupling::bs)
    .def_readonly("cav", &Coupling::cav)
    .def_property_readonly("state", &Coupling::state)
    // Functions
    .def("ImDA", &Coupling::ImDA)
    .def("d_ImDA", &Coupling::d_ImDA)
    .def("ImDA_int", &Coupling::ImDA_int)
    .def("mode_coupling", &Coupling::mode_coupling)
    .def("d_mode_coupling", &Coupling::d_mode_coupling)
    .def("coupling_I", &Coupling::coupling_I)
    .def("pi0_elems", &Coupling::pi0_elems)
    .def("photon_se", &Coupling::photon_se, "Photon self energy matrix")
    .def("_photon_se_or_deriv", &Coupling::_photon_se_or_deriv)
    .def("d_photon_se", &Coupling::d_photon_se)
    .def("photon_se_int", &Coupling::photon_se_int)
    .def("photon_se_mode", &Coupling::photon_se_mode)
    .def("d_photon_se_mode", &Coupling::d_photon_se_mode)
    .def("Z", &Coupling::Z)
    .def("wf_renorm", &Coupling::wf_renorm)
    .def(py::pickle(
      [](const Coupling& c) {
        auto sysp = c.bs.state.sys.pickle();
        auto statep = c.bs.state.pickle();
        auto bsp = c.bs.pickle();
        auto cavp = c.cav.pickle();
        return py::make_tuple(sysp, statep, bsp, cavp);
      },
      [](py::tuple t) {
        auto sys = System::unpickle(t[0].cast<System::pickle_type>());
        auto state = State::unpickle(sys, t[1].cast<State::pickle_type>());
        auto bs = BS::unpickle(state, t[2].cast<BS::pickle_type>());
        auto cav = Cavity::unpickle(t[3].cast<Cavity::pickle_type>());
        return Coupling(bs, cav);
      }));
}
