// clang-format off
#include <nanobind/nanobind.h>

// cland-format on
#include <SHARPlib/thermo.h>

namespace nb = nanobind;
// clang-format on
NB_MODULE(thermo, m) {
    m.doc() =
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Thermodynamic Routines";

    nb::enum_<sharp::adiabat>(m, "adiabat")
        .value("pseudo_liq", sharp::adiabat::pseudo_liq,
               "Use a pseudoadiabatic parcel with liquid only processes")
        .value("adiab_liq", sharp::adiabat::adiab_liq,
               "Use an adiabatic parcel with liquid only processes")
        .value("pseudo_ice", sharp::adiabat::pseudo_ice,
               "Use a pseudoadiabatic parcel with liquid and ice processes")
        .value("adiab_ice", sharp::adiabat::adiab_ice,
               "Use an adiabatic parcel with liquid and ice processes");

    m.def("wobf", &sharp::wobf, nb::arg("temperature"));

    m.def("lcl_temperature", &sharp::lcl_temperature,
          "Compute the Lifted Condensation Level (LCL) temperature (K) using "
          "the air temperature (K) and dewpoint temperature(K)");
}
