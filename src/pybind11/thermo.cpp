// clang-format off
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// cland-format on
#include <SHARPlib/thermo.h>

namespace py = pybind11;
// clang-format on
PYBIND11_MODULE(thermo, m) {
    m.doc() =
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Thermodynamic Routines";
    m.def("lcl_temperature", &sharp::lcl_temperature,
          "Compute the Lifted Condensation Level (LCL) temperature (K) using "
          "the air temperature (K) and dewpoint temperature(K)");
    m.def("lcl_temperature", py::vectorize(sharp::lcl_temperature),
          "Compute the Lifted Condensation Level (LCL) temperature (K) using "
          "the air temperature (K) and dewpoint temperature(K)");
}
