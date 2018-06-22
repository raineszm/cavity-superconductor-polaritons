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
         "double_root_tol"_a = 1e-17)
    .def(py::pickle(
      [](const Polariton& p) {
        auto sysp = p.coupling.bs.state.sys.pickle();
        auto statep = p.coupling.bs.state.pickle();
        auto bsp = p.coupling.bs.pickle();
        auto cavp = p.coupling.cav.pickle();
        return py::make_tuple(sysp, statep, bsp, cavp, p.pickle());
      },
      [](py::tuple t) {
        auto sys = System::unpickle(t[0].cast<System::pickle_type>());
        auto state = State::unpickle(sys, t[1].cast<State::pickle_type>());
        auto bs = BS::unpickle(state, t[2].cast<BS::pickle_type>());
        auto cav = Cavity::unpickle(t[3].cast<Cavity::pickle_type>());
        auto c = Coupling(bs, cav);
        return Polariton::unpickle(c, t[4].cast<Polariton::pickle_type>());
      }));
  ;

  py::class_<ModePolariton>(m, "ModePolariton", polariton)
    .def(py::init<const Coupling&, double, double>(),
         "coupling"_a,
         "paraX"_a = 1.,
         "dipoleX"_a = 1.)
    .def("hamiltonian", &ModePolariton::hamiltonian)
    .def("bands", &ModePolariton::bands)
    .def(py::pickle(
      [](const ModePolariton& p) {
        auto sysp = p.coupling.bs.state.sys.pickle();
        auto statep = p.coupling.bs.state.pickle();
        auto bsp = p.coupling.bs.pickle();
        auto cavp = p.coupling.cav.pickle();
        return py::make_tuple(sysp, statep, bsp, cavp, p.pickle());
      },
      [](py::tuple t) {
        auto sys = System::unpickle(t[0].cast<System::pickle_type>());
        auto state = State::unpickle(sys, t[1].cast<State::pickle_type>());
        auto bs = BS::unpickle(state, t[2].cast<BS::pickle_type>());
        auto cav = Cavity::unpickle(t[3].cast<Cavity::pickle_type>());
        auto c = Coupling(bs, cav);
        return ModePolariton::unpickle(c, t[4].cast<Polariton::pickle_type>());
      }));
}
