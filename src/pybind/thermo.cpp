// clang-format off
#include <pybind11/pybind11.h>

// cland-format on
#include <SHARPlib/thermo.h>


PYBIND11_MODULE(_thermo, m) {
    m.doc() = "Sounding and Hodograph Analysis and Research Program Library (SHARPlib) :: Thermodynamic Routines";
    m.def("lcl_temperature", &sharp::lcl_temperature);
}
