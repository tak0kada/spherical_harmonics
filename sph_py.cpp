#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/detail/common.h>
#include "sph.hpp"

namespace py = pybind11;

PYBIND11_MODULE(cppsph, m)
{
    m.doc() = "C++ Spherical Harmonic Function";
    m.def("sph_harm", &sph::sph_harm,
          py::arg("l"), py::arg("m"), py::arg("theta"), py::arg("phi"));
}
