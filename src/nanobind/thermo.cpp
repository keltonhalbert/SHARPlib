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

    m.def("wobf", &sharp::wobf, nb::arg("temperature"),
          R"pbdoc(
Computes the difference between the wet-bulb potential temperatures
(Kelvin) for saturated and dry air, given the temperature (Kelvin).

The Wobus Function (wobf) is defined as the difference between the 
wet-bulb potential temperature for saturated air (WBPTS) and the 
wet-bulb potential temperature for dry air (WBPTD) given an iput
air temperature. This function is used in the computation of moist 
adiabatic ascent as part of the Wobus parcel lifter (lifter_wobus).

WOBF(T) := WBPTS - WBPTD

Although WBPTS and WBPTD are functions of both pressure and 
temperature, it is assumed their difference is a function of 
temperature only. The difference is also proportional to the 
heat imparted to a parcel of air.

This function uses a polynomial approximation to the wobus function, 
fitted to value in Table 78 of PP.319-322 of the Smithsonian Meteorological
Table by Rolan List (6th Revised Edition). Herman Wobus, a mathematician 
for the Navy Weather Research Facility in Norfolk, VA computed these 
coefficients a very long time ago, as he was retired as of the time of 
the documentation found on this routine written in 1981.

It was shown by Robert Davies-Jones (2007) that the Wobus function has
a slight dependence on pressure, which results in errors of up to 1.2
Kelvin in the temperature of a lifted parcel. 

Parameters:
    temperature: The air temperature (Kelvin)

Returns:
    The Wobus function temperature (Kelvin)

            )pbdoc");

    m.def(
        "wobf",
        [](const_prof_arr_t tmpk_arr) {
            // Get a high performance view
            // into the array for access
            auto tmpk = tmpk_arr.view();

            float *wobf_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                wobf_arr[k] = sharp::wobf(tmpk(k));
            }

            nb::capsule owner(wobf_arr,
                              [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wobf_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("tmpk_arr"),
        R"pbdoc(
Computes the difference between the wet-bulb potential temperatures
(Kelvin) for saturated and dry air, given the temperature (Kelvin).

The Wobus Function (wobf) is defined as the difference between the 
wet-bulb potential temperature for saturated air (WBPTS) and the 
wet-bulb potential temperature for dry air (WBPTD) given an iput
air temperature. This function is used in the computation of moist 
adiabatic ascent as part of the Wobus parcel lifter (lifter_wobus).

WOBF(T) := WBPTS - WBPTD

Although WBPTS and WBPTD are functions of both pressure and 
temperature, it is assumed their difference is a function of 
temperature only. The difference is also proportional to the 
heat imparted to a parcel of air.

This function uses a polynomial approximation to the wobus function, 
fitted to value in Table 78 of PP.319-322 of the Smithsonian Meteorological
Table by Rolan List (6th Revised Edition). Herman Wobus, a mathematician 
for the Navy Weather Research Facility in Norfolk, VA computed these 
coefficients a very long time ago, as he was retired as of the time of 
the documentation found on this routine written in 1981.

It was shown by Robert Davies-Jones (2007) that the Wobus function has
a slight dependence on pressure, which results in errors of up to 1.2
Kelvin in the temperature of a lifted parcel. 

Parameters:
    tmpk_arr: The 1D array of air temperatures (Kelvin)

Returns:
    The 1D array of Wobus function temperatures (Kelvin)

            )pbdoc");

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
            // Get a high performance view
            // into these arrays for access
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();
            if ((tmpk.shape(0) != dwpk.shape(0))) {
                throw nb::buffer_error(
                    "tmpk_arr and dwpk_arr must have the same size!");
            }

            float *lcl_tmpk_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                lcl_tmpk_arr[k] = sharp::lcl_temperature(tmpk(k), dwpk(k));
            }

            nb::capsule owner(lcl_tmpk_arr,
                              [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                lcl_tmpk_arr, {tmpk.shape(0)}, owner);
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

    m.def("vapor_pressure", &sharp::vapor_pressure, nb::arg("pressure"),
          nb::arg("temperature"),
          R"pbdoc(
Compute the vapor pressure with respect to liquid water.

The vapor pressure (or saturation vapor pressure) is computed with 
respect to liquid water when the dewpoint temperature (or air temperature)
is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
is used as a minimum floor value for extremely cold temperatures at low
pressures, and is consistent with how vapor pressure is treated in CM1. 

This function uses the formulation by Bolton (1980), and is
accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

Parameters:
    pressure: The air pressure (Pa)
    temperature: The air temperature (K) or dewpoint temperature (K)

Returns: 
    The vapor pressure (Pa) given dewpoint temperature, or the 
    saturation vapor pressure given the air temperature. 
    )pbdoc");

    m.def(
        "vapor_pressure",
        [](const_prof_arr_t pres_arr, const_prof_arr_t tmpk_arr) {
            auto tmpk = tmpk_arr.view();
            auto pres = pres_arr.view();
            if ((tmpk.shape(0) != pres.shape(0))) {
                throw nb::buffer_error(
                    "tmpk_arr and pres_arr must have the same size!");
            }

            float *vappres_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                vappres_arr[k] = sharp::vapor_pressure(pres(k), tmpk(k));
            }

            nb::capsule owner(vappres_arr,
                              [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                vappres_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("pres_arr"), nb::arg("tmpk_arr"),
        R"pbdoc(
Compute the vapor pressure with respect to liquid water.

The vapor pressure (or saturation vapor pressure) is computed with 
respect to liquid water when the dewpoint temperature (or air temperature)
is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
is used as a minimum floor value for extremely cold temperatures at low
pressures, and is consistent with how vapor pressure is treated in CM1. 

This function uses the formulation by Bolton (1980), and is
accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

Parameters:
    pres_arr: The 1D array of air pressure (Pa)
    tmpk_arr: The 1D array of air temperature (K) or dewpoint temperature (K)

Returns: 
    The 1D array of vapor pressure (Pa) given dewpoint temperature, or the 
    saturation vapor pressure given the air temperature. 
    )pbdoc");

    m.def("vapor_pressure_ice", &sharp::vapor_pressure_ice, nb::arg("pressure"),
          nb::arg("temperature"),
          R"pbdoc(
Compute the vapor pressure with respect to ice.

The vapor pressure (or saturation vapor pressure) is computed with 
respect to ice when the dewpoint temperature (or air temperature)
is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
is used as a minimum floor value for extremely cold temperatures at low
pressures, and is consistent with how vapor pressure is treated in CM1. 

This function uses the formulation by Bolton (1980), and is
accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

Parameters:
    pressure: The air pressure (Pa)
    temperature: The air temperature (K) or dewpoint temperature (K)

Returns: 
    The vapor pressure (Pa) given dewpoint temperature, or the 
    saturation vapor pressure given the air temperature. 
    )pbdoc");

    m.def(
        "vapor_pressure_ice",
        [](const_prof_arr_t pres_arr, const_prof_arr_t tmpk_arr) {
            auto tmpk = tmpk_arr.view();
            auto pres = pres_arr.view();
            if ((tmpk.shape(0) != pres.shape(0))) {
                throw nb::buffer_error(
                    "tmpk_arr and pres_arr must have the same size!");
            }

            float *vappres_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                vappres_arr[k] = sharp::vapor_pressure_ice(pres(k), tmpk(k));
            }

            nb::capsule owner(vappres_arr,
                              [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                vappres_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("pres_arr"), nb::arg("tmpk_arr"),
        R"pbdoc(
Compute the vapor pressure with respect to ice.

The vapor pressure (or saturation vapor pressure) is computed with 
respect to ice when the dewpoint temperature (or air temperature)
is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
is used as a minimum floor value for extremely cold temperatures at low
pressures, and is consistent with how vapor pressure is treated in CM1. 

This function uses the formulation by Bolton (1980), and is
accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

Parameters:
    pres_arr: The 1D array of air pressure (Pa)
    tmpk_arr: The 1D array of air temperature (K) or dewpoint temperature (K)

Returns: 
    The 1D array of vapor pressure (Pa) given dewpoint temperature, or the 
    saturation vapor pressure given the air temperature. 
    )pbdoc");
}
