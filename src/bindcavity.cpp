#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "cavity.h"

namespace py = pybind11;
using namespace py::literals;

void
bind_cavity(py::module& m)
{
  py::class_<Cavity>(m, "Cavity")
    .def(py::init<double>(), "omega0"_a)
    // Attributes
    .def_readwrite("omega0", &Cavity::omega0)
    .def("inv_gf", &Cavity::inv_gf)
    .def("d_inv_gf", &Cavity::d_inv_gf)
    .def("omega", &Cavity::omega)
    .def("polarizations", &Cavity::polarizations)
    .def(py::pickle([](const Cavity& c) { return py::make_tuple(c.pickle()); },
                    [](py::tuple t) {
                      return Cavity::unpickle(t[0].cast<Cavity::pickle_type>());
                    }));
  m.attr("GPAR") = GPAR;
  m.attr("ALPHA") = ALPHA;
  m.attr("C") = C;
}
