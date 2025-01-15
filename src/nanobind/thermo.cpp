// clang-format off
#include <nanobind/nanobind.h>

// cland-format on
#include <SHARPlib/thermo.h>
#include "sharplib_types.h"

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

    m.def("lcl_temperature", &sharp::lcl_temperature, nb::arg("temperature"),
          nb::arg("dewpoint"),
          R"pbdoc(
            Compute the Lifted Condensation Level (LCL) temperature (Kelvin) 
            given an air temperature (Kelvin) and dewpoint temperature (Kelvin).

            The LCL temperature is computed as in Bolton (1980) eq 15, and is
            considered to be within a 10th of a degree of the more exact 
            iterative formula.

            Parameters:
                temperature: The air temperature (Kelvin)
                dewpoint: The dewpoint temperature (Kelvin)

            Returns:
                The LCL temperature (Kelvin)

            )pbdoc");

    m.def(
        "lcl_temperature",
        [](const_prof_arr_t tmpk_arr, const_prof_arr_t dwpk_arr) {
            if ((tmpk_arr.size() != dwpk_arr.size())) {
                throw nb::buffer_error(
                    "tmpk_arr and dwpk_arr must have the same size!");
            }
            // Get a high performance view
            // into these arrays for access
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();

            float *lcl_tmpk_arr = new float[tmpk_arr.size()];
            for (size_t k = 0; k < tmpk_arr.size(); ++k) {
                lcl_tmpk_arr[k] = sharp::lcl_temperature(tmpk(k), dwpk(k));
            }

            nb::capsule owner(lcl_tmpk_arr,
                              [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                lcl_tmpk_arr, {tmpk_arr.size()}, owner);
        },
        nb::arg("tmpk_arr"), nb::arg("dwpk_arr"),
        R"pbdoc(
        Compute the Lifted Condensation Level (LCL) temperature (Kelvin) 
        given an array of air temperatures (Kelvin) and dewpoint 
        temperatures (Kelvin).

        The LCL temperatures are computed as in Bolton (1980) eq 15, and is
        considered to be within a 10th of a degree of the more exact 
        iterative formula.

        Parameters:
            tmpk_arr: The 1D air temperature array (Kelvin)
            dwpk_arr: The 1D dewpoint temperature array (Kelvin)

        Returns:
            The 1D array of LCL temperatures (Kelvin)
    
        )pbdoc");
}
