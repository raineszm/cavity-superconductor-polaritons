#include <pybind11/pybind11.h>

#include "bs.h"
#include "state.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_bs(py::module& m)
{
  py::class_<BS>(m, "BS")
    .def(py::init<double, const State&>(), "mass"_a, "state"_a)
    // Attributes
    .def_readonly("state", &BS::state)
    .def_readwrite("mass", &BS::mass)
    .def("inv_gf_int", &BS::inv_gf_int)
    .def("d_inv_gf_int", &BS::d_inv_gf_int)
    .def("inv_gf", &BS::inv_gf)
    .def("d_inv_gf", &BS::d_inv_gf)
    .def("I", &BS::I)
    .def_property_readonly("root", &BS::root)
    .def_property_readonly("M", &BS::M)
    .def(py::pickle(
      [](const BS& bs) {
        auto sysp = bs.state.sys.pickle();
        auto statep = bs.state.pickle();
        return py::make_tuple(sysp, statep, bs.pickle());
      },
      [](py::tuple t) {
        auto sys = System::unpickle(t[0].cast<System::pickle_type>());
        auto state = State::unpickle(sys, t[1].cast<State::pickle_type>());
        return BS::unpickle(state, t[2].cast<BS::pickle_type>());
      }));
}
