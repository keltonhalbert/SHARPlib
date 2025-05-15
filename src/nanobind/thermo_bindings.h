#ifndef SHARPLIB_THERMO_BINDINGS
#define SHARPLIB_THERMO_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

#include <SHARPlib/thermo.h>
#include <SHARPlib/parcel.h>
#include "sharplib_types.h"

namespace nb = nanobind;


inline void make_thermo_bindings (nb::module_ m) {
    nb::module_ m_therm = m.def_submodule("thermo", 
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Thermodynamic Routines");

    nb::enum_<sharp::adiabat>(m_therm, "adiabat")
        .value("pseudo_liq", sharp::adiabat::pseudo_liq,
            "Use a pseudoadiabatic parcel with liquid only processes")
        .value("adiab_liq", sharp::adiabat::adiab_liq,
            "Use an adiabatic parcel with liquid only processes")
        .value("pseudo_ice", sharp::adiabat::pseudo_ice,
            "Use a pseudoadiabatic parcel with liquid and ice processes")
        .value("adiab_ice", sharp::adiabat::adiab_ice,
            "Use an adiabatic parcel with liquid and ice processes");

    m_therm.def("wobf", &sharp::wobf, nb::arg("temperature"),
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

    m_therm.def(
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

    m_therm.def("lcl_temperature", &sharp::lcl_temperature, nb::arg("temperature"),
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

    m_therm.def(
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

    m_therm.def("vapor_pressure", &sharp::vapor_pressure, nb::arg("pressure"),
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

    m_therm.def(
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

    m_therm.def("vapor_pressure_ice", &sharp::vapor_pressure_ice, nb::arg("pressure"),
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

    m_therm.def(
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

    m_therm.def("temperature_at_mixratio", &sharp::temperature_at_mixratio,
        nb::arg("wv_mixratio"), nb::arg("pressure"),
        R"pbdoc(
Computes the temperature (K) of air at the given water vapor mixing ratio
(kg/kg) and air pressure (Pa). Can be used to compute the dewpoint temperature 
from mixing ratio. 

The routine is implemented as in Bolton (1980) eq 11, and is considered to be 
accurate to 0.03 K for -35C <= T <= 35C. 

Parameters:
    wv_mixratio: The water vapor mixing ratio (kg/kg)
    pressure: The air pressure (Pa)

Returns:
    The temperature (K) of an air parcel at a given mixing ratio and pressure.
    )pbdoc");

    m_therm.def(
        "temperature_at_mixratio",
        [](const_prof_arr_t mixr_arr, const_prof_arr_t pres_arr) {
            auto mixr = mixr_arr.view();
            auto pres = pres_arr.view();
            if ((mixr.shape(0) != pres.shape(0))) {
                throw nb::buffer_error(
                    "mixr_arr and pres_arr must have the same size!");
            }

            float *dwpk_arr = new float[mixr.shape(0)];
            for (size_t k = 0; k < mixr.shape(0); ++k) {
                dwpk_arr[k] = sharp::temperature_at_mixratio(mixr(k), pres(k));
            }

            nb::capsule owner(dwpk_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                dwpk_arr, {mixr.shape(0)}, owner);
        },
        nb::arg("mixr_arr"), nb::arg("pres_arr"),
        R"pbdoc(
Computes the temperature (K) of air at the given water vapor mixing ratio
(kg/kg) and air pressure (Pa). Can be used to compute the dewpoint temperature 
from mixing ratio. 

The routine is implemented as in Bolton (1980) eq 11, and is considered to be 
accurate to 0.03 K for -35C <= T <= 35C. 

Parameters:
    mixr_arr: The 1D array of water vapor mixing ratio (kg/kg)
    pres_arr: The 1D array of air pressure (Pa)

Returns:
    The 1D array of temperature (K) of an air parcel at a given mixing ratio and pressure.

    )pbdoc");

    m_therm.def("theta_level", &sharp::theta_level, nb::arg("potential_temperature"),
        nb::arg("temperature"),
        R"pbdoc(
Computes the pressure level (Pa) of a parcel given the potential temperature (K) and air 
temperature (K).

Params:
    potential_temperature: The potential temperature, or theta (K)
    temperature: The air temperature (K)

Returns:
    The pressure level (Pa) corresponding to the potential temperature and air temperature.
    )pbdoc");

    m_therm.def(
        "theta_level",
        [](const_prof_arr_t theta_arr, const_prof_arr_t tmpk_arr) {
            auto theta = theta_arr.view();
            auto tmpk = tmpk_arr.view();
            if ((theta.shape(0) != tmpk.shape(0))) {
                throw nb::buffer_error(
                    "theta_arr and tmpk_arr must have the same size!");
            }
            float *pres_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                pres_arr[k] = sharp::theta_level(theta(k), tmpk(k));
            }

            nb::capsule owner(pres_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                pres_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("theta_arr"), nb::arg("tmpk_arr"),
        R"pbdoc(
Computes the pressure level (Pa) of a parcel given the potential temperature (K) and air 
temperature (K).

Params:
    theta_arr: The 1D array of potential temperature, or theta (K)
    tmpk_arr: The 1D array of air temperature (K)

Returns:
    The 1D pressure level (Pa) corresponding to the potential temperature and air temperature.

    )pbdoc");

    m_therm.def("theta", &sharp::theta, nb::arg("pressure"), nb::arg("temperature"),
        nb::arg("ref_pressure") = sharp::THETA_REF_PRESSURE,
        R"pbdoc(
Computes the potential temperature (K), or theta, given the air pressure (Pa), air temperature (K),
and a reference pressure (default value is 100000 Pa).

Parameters:
    pressure: The air pressure (Pa)
    temperature: The air temperature (K)
    ref_pressure (optional): The reference pressure (default 100000.0 Pa)

Returns:
    The potential temperature (K), or theta
    )pbdoc");

    m_therm.def(
        "theta",
        [](const_prof_arr_t pres_arr, const_prof_arr_t tmpk_arr,
        const float ref_pres) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            if ((pres.shape(0) != tmpk.shape(0))) {
                throw nb::buffer_error(
                    "pres_arr and tmpk_arr must have the same size!");
            }
            float *theta_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                theta_arr[k] = sharp::theta(pres(k), tmpk(k), ref_pres);
            }

            nb::capsule owner(theta_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                theta_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("pres_arr"), nb::arg("tmpk_arr"),
        nb::arg("ref_pressure") = sharp::THETA_REF_PRESSURE,
        R"pbdoc(
Computes the potential temperature (K), or theta, given the air pressure (Pa), air temperature (K),
and a reference pressure (default value is 100000 Pa).

Parameters:
    pressure: The 1D array of air pressure (Pa)
    temperature: The 1D array of air temperature (K)
    ref_pressure (optional): The reference pressure (default 100000.0 Pa)

Returns:
    The 1D array of potential temperature (K), or theta
    )pbdoc");

    m_therm.def("mixratio", static_cast<float (*)(float)>(&sharp::mixratio),
        nb::arg("q"),
        R"pbdoc(
Compute the water vapor mixing ratio (kg/kg) from specific humidity (kg/kg).

Parameters:
    q: The specific humidity (kg/kg)

Returns:
    The water vapor mixing ratio (kg/kg)
    )pbdoc");

    m_therm.def("mixratio", static_cast<float (*)(float, float)>(&sharp::mixratio),
        nb::arg("pressure"), nb::arg("temperature"),
        R"pbdoc(
Compute the water vapor mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
If given the air temperature, this is the saturation mixing ratio. If given the dewpoint 
temperature, tis is the mixing ratio. 

Parameters:
    pressure: The air pressure (Pa)
    temperature: The air temperature (K)

Returns:
    The water vapor mixing ratio (kg/kg)
    )pbdoc");

    m_therm.def(
        "mixratio",
        [](const_prof_arr_t spfh_arr) {
            auto spfh = spfh_arr.view();

            float *mixr_arr = new float[spfh.shape(0)];
            for (size_t k = 0; k < spfh.shape(0); ++k) {
                mixr_arr[k] = sharp::mixratio(spfh(k));
            }

            nb::capsule owner(mixr_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                mixr_arr, {spfh.shape(0)}, owner);
        },
        nb::arg("spfh_arr"),
        R"pbdoc(
Compute the water vapor mixing ratio (kg/kg) from the specific humidity (kg/kg).

Parameters:
    spfh_arr: The 1D array of specific humidity (kg/kg)

Returns:
    The 1D array of water vapor mixing ratio (kg/kg)

    )pbdoc");

    m_therm.def(
        "mixratio",
        [](const_prof_arr_t pres_arr, const_prof_arr_t tmpk_arr) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            if ((pres.shape(0) != tmpk.shape(0))) {
                throw nb::buffer_error(
                    "pres_arr and tmpk_arr must have the same size!");
            }
            float *mixr_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                mixr_arr[k] = sharp::mixratio(pres(k), tmpk(k));
            }

            nb::capsule owner(mixr_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                mixr_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("pres_arr"), nb::arg("tmpk_arr"),
        R"pbdoc(
Compute the water vapor mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
If given the air temperature, this is the saturation mixing ratio. If given the dewpoint 
temperature, this is the mixing ratio. 

Parameters:
    pres_arr: The 1D array of air pressure (Pa)
    tmpk_arr: The 1D array of air temperature (K)

Returns:
    The 1D array of water vapor mixing ratio (kg/kg)

    )pbdoc");

    m_therm.def("mixratio_ice", &sharp::mixratio_ice, nb::arg("pressure"),
        nb::arg("temperature"),
        R"pbdoc(
Compute the ice water mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
If given the air temperatuer, this is the saturation mixing ratio. If given the dewpoint 
temperature, this is the mixing ratio. 

Parameters:
    pressure: The air pressure (Pa)
    temperature: The air temperature (K)
Returns:
    The ice water mixing ratio (kg/kg)
    )pbdoc");

    m_therm.def(
        "mixratio_ice",
        [](const_prof_arr_t pres_arr, const_prof_arr_t tmpk_arr) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            if ((pres.shape(0) != tmpk.shape(0))) {
                throw nb::buffer_error(
                    "pres_arr and tmpk_arr must have the same size!");
            }
            float *mixr_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                mixr_arr[k] = sharp::mixratio_ice(pres(k), tmpk(k));
            }

            nb::capsule owner(mixr_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                mixr_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("pres_arr"), nb::arg("tmpk_arr"),
        R"pbdoc(
Compute the ice water mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
If given the air temperatuer, this is the saturation mixing ratio. If given the dewpoint 
temperature, this is the mixing ratio. 

Parameters:
    pres_arr: The 1D array of air pressure (Pa)
    tmpl_arr: The 1D array of air temperature (K)
Returns:
    The ice water mixing ratio (kg/kg)

    )pbdoc");

    m_therm.def("specific_humidity", &sharp::specific_humidity, nb::arg("rv"),
        R"pbdoc(
Compute the specific humidity (kg/kg) from a mixing ratio (kg/kg).

Parameters:
    rv: The water vapor mixing ratio (kg/kg)

Returns:
    The specific humidity (kg/kg)
    )pbdoc");

    m_therm.def(
        "specific_humidity",
        [](const_prof_arr_t mixr_arr) {
            auto mixr = mixr_arr.view();
            float *spfh_arr = new float[mixr.shape(0)];
            for (size_t k = 0; k < mixr.shape(0); ++k) {
                spfh_arr[k] = sharp::specific_humidity(mixr(k));
            }

            nb::capsule owner(spfh_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                spfh_arr, {mixr.shape(0)}, owner);
        },
        nb::arg("mixr_arr"),
        R"pbdoc(
Compute the specific humidity (kg/kg) from a mixing ratio (kg/kg).

Parameters:
    mixr_arr: The 1D array of water vapor mixing ratios (kg/kg)
    )pbdoc");

    m_therm.def(
        "relative_humidity",
        &sharp::relative_humidity, 
        nb::arg("pressure"), nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the relative humidity of water vapor with respect to liquid water.
NOTE: The pressure variable is only used as a sanity check when computing 
vapor pressure at extremely low pressures and temperatures. If you do not 
want or need this behavior, you can pass the THETA_REF_PRESSURE constant
in place of an air pressure.

Parameters:
    pressure: the ambient air pressure (Pa)
    temperature: the ambient air temperature (K)
    dewpoint: the ambient dewpoint temperature (K)

Returns:
    Relative Humidity (fraction, unitless)
    )pbdoc"
    );
    
    m_therm.def(
        "relative_humidity",
        [](
            const_prof_arr_t pres_arr,
            const_prof_arr_t tmpk_arr,
            const_prof_arr_t dwpk_arr
        ) {

            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();

            float* relh = new float[pres.shape(0)];

            for (std::ptrdiff_t k = 0; k < pres.shape(0); ++k) {
                relh[k] = sharp::relative_humidity(pres(k), tmpk(k), dwpk(k));
            }

            nb::capsule owner(relh,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                relh, {tmpk.shape(0)}, owner);
        },
        nb::arg("pressure"), nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the relative humidity of water vapor with respect to liquid water.
NOTE: The pressure variable is only used as a sanity check when computing 
vapor pressure at extremely low pressures and temperatures. If you do not
want or need this behavior, you canpass in an array of THETA_REF_PRESSURE
in place of air pressure.

Paramters:
    pressure: 1D NumPy array of ambient air pressure Pa)
    temperature: 1D NumPy array of ambient air temperature (K)
    dewpoint: 1D NumPy array of ambient dewpoint temperature (K)

Returns:
    1D NumPy array of relative humidity values (fraction, unitless)
    )pbdoc"

    );

    m_therm.def("virtual_temperature", &sharp::virtual_temperature,
        nb::arg("temperature"), nb::arg("rv"), nb::arg("rl") = 0.0f,
        nb::arg("ri") = 0.0f,
        R"pbdoc(
Returns the virtual temperature in Kelvin given the dry-bulb 
temperature (Kelvin), the water vapor mixing ratio (kg/kg), the 
liquid water mixing ratio (kg/kg), and the ice water mixing ratios
(kg/kg). The liquid and ice water mixing ratios have default values 
of zero, if unspecified. 

Parameters:
    temperature: The dry-bulb temperature (K)
    rv: The water vapor mixing ratio (kg/kg)
    rl: The liquid water mixing ratio (kg/kg) 
    ri: The ice water mixing ratio (kg/kg)

Returns:
    The virtual temperature (K)
    )pbdoc");

    m_therm.def(
        "virtual_temperature",
        [](const_prof_arr_t tmpk_arr, const_prof_arr_t rv_arr,
        const_prof_arr_t rl_arr, const_prof_arr_t ri_arr) {
            auto tmpk = tmpk_arr.view();
            auto rv = rv_arr.view();
            float *vtmp_arr;

            // tmpk and rv are always defined, check their
            // sizes and ensure they're equal
            if (tmpk.shape(0) != rv.shape(0)) {
                throw nb::buffer_error(
                    "tmpk_arr and rv_arr must have the same size!");
            }

            // check if rl_arr and ri_arr are not none
            if (rl_arr.is_valid() && ri_arr.is_valid()) {
                auto rl = rl_arr.view();
                auto ri = ri_arr.view();
                // ensure they're the same shape/size
                if ((tmpk.shape(0) != rl.shape(0)) ||
                    (tmpk.shape(0) != ri.shape(0))) {
                    throw nb::buffer_error(
                        "tmpk_arr and rv_arr must have the same size!");
                }
                vtmp_arr = new float[tmpk.shape(0)];

                for (size_t k = 0; k < tmpk.shape(0); ++k) {
                    vtmp_arr[k] = sharp::virtual_temperature(tmpk(k), rv(k),
                                                            rl(k), ri(k));
                }
            } else if (rl_arr.is_valid()) {
                auto rl = rl_arr.view();
                if (tmpk.shape(0) != rl.shape(0)) {
                    throw nb::buffer_error(
                        "tmpk_arr, rv_arr, and rl_arr must have the same "
                        "size!");
                }
                vtmp_arr = new float[tmpk.shape(0)];
                for (size_t k = 0; k < tmpk.shape(0); ++k) {
                    vtmp_arr[k] =
                        sharp::virtual_temperature(tmpk(k), rv(k), rl(k));
                }
            } else if (ri_arr.is_valid()) {
                auto ri = ri_arr.view();
                if (tmpk.shape(0) != ri.shape(0)) {
                    throw nb::buffer_error(
                        "tmpk_arr, rv_arr, and ri_arr must have the same "
                        "size!");
                }
                vtmp_arr = new float[tmpk.shape(0)];
                for (size_t k = 0; k < tmpk.shape(0); ++k) {
                    vtmp_arr[k] =
                        sharp::virtual_temperature(tmpk(k), rv(k), 0.0, ri(k));
                }
            } else {
                vtmp_arr = new float[tmpk.shape(0)];
                for (size_t k = 0; k < tmpk.shape(0); ++k) {
                    vtmp_arr[k] = sharp::virtual_temperature(tmpk(k), rv(k));
                }
            }

            nb::capsule owner(vtmp_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                vtmp_arr, {tmpk.shape(0)}, owner);
        },
        nb::arg("tmpk_arr"), nb::arg("rv_arr"), nb::arg("rl_arr") = nb::none(),
        nb::arg("ri_arr") = nb::none(),
        R"pbdoc(
Returns the virtual temperature in Kelvin given the dry-bulb 
temperature (Kelvin), the water vapor mixing ratio (kg/kg), the 
liquid water mixing ratio (kg/kg), and the ice water mixing ratios
(kg/kg). The liquid and ice water mixing ratios have default values 
of zero, if unspecified. 

Parameters:
    tmpk_arr: The 1D array of dry-bulb temperature (K)
    rv_arr: The 1D array of water vapor mixing ratio (kg/kg)
    rl_arr: (optional) The 1D array of liquid water mixing ratio (kg/kg) 
    ri_arr: (optional) The 1D array of ice water mixing ratio (kg/kg)

Returns:
    The 1D array of virtual temperature (K)

    )pbdoc");

    // Skipping saturated_lift, because I don't intend
    // for it to be called from Python directly.

    m_therm.def("wetlift", &sharp::wetlift, nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("lifted_pressure"),
        R"pbdoc(
Compute the temperature of a parcel lifted moist adiabatically to a new level. 

With a given parcel defined by a pressure (Pa) and temperature (K), lift it 
moist adiabatically to a new pressure level (Pa) and return the temperature of 
the parcel at that level. 

This function relies on the Wobus Function (thermo.wobf), and it was shown by 
Robert Davies-Jones (2007) that the WObus function has a slight dependence on 
pressure, which results in errors of up to 1.2 K in the temperature of a lifted 
parcel. 

Parameters:
    pressure: The air pressure (Pa)
    temperature: The saturated air temperature (K)
    lifted_pressure: The new pressure level to lift to (Pa)

Returns:
    The new temperature (K) when lifted moist adiabatically to the new pressure level
        )pbdoc");

    m_therm.def(
        "wetlift",
        [](const_prof_arr_t pres_arr, const_prof_arr_t tmpk_arr,
        const_prof_arr_t lifted_pres_arr) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto lifted_pres = lifted_pres_arr.view();

            if ((pres.shape(0) != tmpk.shape(0)) ||
                (pres.shape(0) != lifted_pres.shape(0))) {
                throw nb::buffer_error(
                    "pres_arr, tmpk_arr, and lifted_pres_arr must have the "
                    "same sizes!");
            }
            float *tmpk_out_arr = new float[tmpk.shape(0)];
            for (size_t k = 0; k < tmpk.shape(0); ++k) {
                tmpk_out_arr[k] =
                    sharp::wetlift(pres(k), tmpk(k), lifted_pres(k));
            }

            nb::capsule owner(tmpk_out_arr,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                tmpk_out_arr, {tmpk.shape(0)}, owner);
        },
        R"pbdoc(
Compute the temperature of a parcel lifted moist adiabatically to a new level. 

With a given parcel defined by a pressure (Pa) and temperature (K), lift it 
moist adiabatically to a new pressure level (Pa) and return the temperature of 
the parcel at that level. 

This function relies on the Wobus Function (thermo.wobf), and it was shown by 
Robert Davies-Jones (2007) that the WObus function has a slight dependence on 
pressure, which results in errors of up to 1.2 K in the temperature of a lifted 
parcel. 

Parameters:
    pressure: The 1D array of air pressures (Pa)
    temperature: The 1D array of saturated air temperatures (K)
    lifted_pressure: The 1D array of new pressure levels to lift to (Pa)

Returns:
    The 1D array of new temperatures (K) when lifted moist adiabatically to the new pressure levels

        )pbdoc");

    m_therm.def(
        "drylift",
        [](float pressure, float temperature, float dewpoint) {
            float pressure_at_lcl, temperature_at_lcl;
            
            sharp::drylift(pressure, temperature, dewpoint, pressure_at_lcl, temperature_at_lcl);

            return nb::make_tuple(pressure_at_lcl, temperature_at_lcl);
        },
        nb::arg("pressure"), nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Given a parcel's initial pressure (Pa), temperature (K), and 
dewpoint temperature (K), lift the parcel dry adiabaitically to its
Lifted Condensation Level and return the resulting LCL pressure (Pa)
and temperature (K).

The LCL temperature is computed using an approximation. See the 
lcl_temperature function documentation for more information.

Parameters:
    pressure: Parcel starting pressure (Pa)
    temperature: Parcel starting temperature (K)
    dewpoint: Parcel starting dewpoint temperature (K)

Returns:
    A tuple of (lcl_pres, lcl_temperature)
    )pbdoc"
    );

    m_therm.def(
        "wetbulb",
        [](
            sharp::lifter_wobus& lifter,
            float pressure, 
            float temperature,
            float dewpoint
        ) {

            return sharp::wetbulb(lifter, pressure, temperature, dewpoint);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, and
dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
temperature and pressure of the LCL, an approximation is used. See the 
lcl_temperature function for further detail. 

After the parcel has reached the LCL, the lifter passed to the function 
lowers the parcel to its initial pressure level along a moist adiabat or 
pseudoadiabat. 

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: The ambient pressure (Pa)
    temperature: The ambient temperature (K)
    dewpoint: The ambient dewpoint temperature (K)

Returns:
    The wetbulb temperature (K)
    )pbdoc"
    );

    m_therm.def(
        "wetbulb",
        [](
            sharp::lifter_cm1& lifter,
            float pressure, 
            float temperature,
            float dewpoint
        ) {

            return sharp::wetbulb(lifter, pressure, temperature, dewpoint);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, and
dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
temperature and pressure of the LCL, an approximation is used. See the 
lcl_temperature function for further detail. 

After the parcel has reached the LCL, the lifter passed to the function 
lowers the parcel to its initial pressure level along a moist adiabat or 
pseudoadiabat. 

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: The ambient pressure (Pa)
    temperature: The ambient temperature (K)
    dewpoint: The ambient dewpoint temperature (K)

Returns:
    The wetbulb temperature (K)
    )pbdoc"
    );

    m_therm.def(
        "wetbulb",
        [](
            sharp::lifter_wobus& lifter,
            const_prof_arr_t pres_arr,
            const_prof_arr_t tmpk_arr,
            const_prof_arr_t dwpk_arr
        ) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();
            
            float* wb_out = new float[pres.shape(0)];

            for (size_t k = 0; k <= pres.shape(0); ++k) {
                wb_out[k] = sharp::wetbulb(lifter, pres(k), tmpk(k), dwpk(k));
            }
            
            nb::capsule owner(wb_out,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wb_out, {tmpk.shape(0)}, owner);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, and
dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
temperature and pressure of the LCL, an approximation is used. See the 
lcl_temperature function for further detail. 

After the parcel has reached the LCL, the lifter passed to the function 
lowers the parcel to its initial pressure level along a moist adiabat or 
pseudoadiabat. 

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: 1D NumPy array of ambient pressures (Pa)
    temperature: 1D NumPy array of ambient temperatureds (K)
    dewpoint: 1D NumPy array of ambient dewpoint temperatures (K)

Returns:
    1D NumPy array of wetbulb temperatures (K)
    )pbdoc"
    );

    m_therm.def(
        "wetbulb",
        [](
            sharp::lifter_cm1& lifter,
            const_prof_arr_t pres_arr,
            const_prof_arr_t tmpk_arr,
            const_prof_arr_t dwpk_arr
        ) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();
            
            float* wb_out = new float[pres.shape(0)];

            for (size_t k = 0; k <= pres.shape(0); ++k) {
                wb_out[k] = sharp::wetbulb(lifter, pres(k), tmpk(k), dwpk(k));
            }
            
            nb::capsule owner(wb_out,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wb_out, {tmpk.shape(0)}, owner);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, and
dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
temperature and pressure of the LCL, an approximation is used. See the 
lcl_temperature function for further detail. 

After the parcel has reached the LCL, the lifter passed to the function 
lowers the parcel to its initial pressure level along a moist adiabat or 
pseudoadiabat. 

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: 1D NumPy array of ambient pressures (Pa)
    temperature: 1D NumPy array of ambient temperatureds (K)
    dewpoint: 1D NumPy array of ambient dewpoint temperatures (K)

Returns:
    1D NumPy array of wetbulb temperatures (K)
    )pbdoc"
    );

    m_therm.def(
        "theta_wetbulb",
        [](
            sharp::lifter_wobus& lifter,
            float pressure, 
            float temperature,
            float dewpoint
        ) {
            return sharp::theta_wetbulb(lifter, pressure, temperature, dewpoint);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, 
and dewpoint temperature to its Lifted Condensation Level (LCL). 
To compute the temperature and pressure of the LCL, an approximation 
is used. See the lcl_temperature fuction for further detail. 

After the parcel has reached the LCL, the lifted passed to the function 
lowers the parcel to the standard parcel reference pressure level 
(1000 hPa) along a moist adiabat.

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: the ambient air pressure (Pa)
    temperature: the ambient air temperature (K)
    dewpoint: the ambient dewpoint temperature

Returns:
    The wet-bulb potential temperature (K)
    )pbdoc"
    );

    m_therm.def(
        "theta_wetbulb",
        [](
            sharp::lifter_cm1& lifter,
            float pressure, 
            float temperature,
            float dewpoint
        ) {
            return sharp::theta_wetbulb(lifter, pressure, temperature, dewpoint);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, 
and dewpoint temperature to its Lifted Condensation Level (LCL). 
To compute the temperature and pressure of the LCL, an approximation 
is used. See the lcl_temperature fuction for further detail. 

After the parcel has reached the LCL, the lifted passed to the function 
lowers the parcel to the standard parcel reference pressure level 
(1000 hPa) along a moist adiabat.

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: the ambient air pressure (Pa)
    temperature: the ambient air temperature (K)
    dewpoint: the ambient dewpoint temperature

Returns:
    The wet-bulb potential temperature (K)
    )pbdoc"
    );
    
    m_therm.def(
        "theta_wetbulb",
        [](
            sharp::lifter_wobus& lifter,
            const_prof_arr_t pres_arr,
            const_prof_arr_t tmpk_arr,
            const_prof_arr_t dwpk_arr
        ) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();

            float* theta_wb = new float[pres.shape(0)];
            
            for (size_t k = 0; k < pres.shape(0); ++k) {
                theta_wb[k] = sharp::theta_wetbulb(lifter, pres(k), tmpk(k), dwpk(k));
            }
            
            nb::capsule owner(theta_wb,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                theta_wb, {tmpk.shape(0)}, owner);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, 
and dewpoint temperature to its Lifted Condensation Level (LCL). 
To compute the temperature and pressure of the LCL, an approximation 
is used. See the lcl_temperature fuction for further detail. 

After the parcel has reached the LCL, the lifted passed to the function 
lowers the parcel to the standard parcel reference pressure level 
(1000 hPa) along a moist adiabat.

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: 1D NumPy array of ambient air pressure (Pa)
    temperature: 1D NumPy array of ambient air temperature (K)
    dewpoint: 1D NumPy array of ambient dewpoint temperature

Returns:
    1D NumPy array of wet-bulb potential temperatures (K)
    )pbdoc"
    );
    
    m_therm.def(
        "theta_wetbulb",
        [](
            sharp::lifter_cm1& lifter,
            const_prof_arr_t pres_arr,
            const_prof_arr_t tmpk_arr,
            const_prof_arr_t dwpk_arr
        ) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();

            float* theta_wb = new float[pres.shape(0)];
            
            for (size_t k = 0; k < pres.shape(0); ++k) {
                theta_wb[k] = sharp::theta_wetbulb(lifter, pres(k), tmpk(k), dwpk(k));
            }
            
            nb::capsule owner(theta_wb,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                theta_wb, {tmpk.shape(0)}, owner);
        },
        nb::arg("lifter"), nb::arg("pressure"),
        nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
temperature (K), and dewpoint temperature (K).

First, it lifts a parcel with the given pressure, temperature, 
and dewpoint temperature to its Lifted Condensation Level (LCL). 
To compute the temperature and pressure of the LCL, an approximation 
is used. See the lcl_temperature fuction for further detail. 

After the parcel has reached the LCL, the lifted passed to the function 
lowers the parcel to the standard parcel reference pressure level 
(1000 hPa) along a moist adiabat.

Parameters:
    lifter: a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure: 1D NumPy array of ambient air pressure (Pa)
    temperature: 1D NumPy array of ambient air temperature (K)
    dewpoint: 1D NumPy array of ambient dewpoint temperature

Returns:
    1D NumPy array of wet-bulb potential temperatures (K)
    )pbdoc"
    );

    m_therm.def(
        "thetae",
        &sharp::thetae,
        nb::arg("pressure"), nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the equivalent potential temperature (K) given the 
air pressure (Pa), air temperature (K), and dewpoint temperature (K).
The equivalent potential temperature is computed as in Bolton 1980. 
It was found in Davies-Jones 2009 to be the most accurate non-iterative 
formulation of theta-e.

 Bolton 1980: https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2
 Davies-Jones 2009: https://doi.org/10.1175/2009MWR2774.1

Parameters:
    pressure: the ambient air pressure (Pa)
    temperature: the ambient air temperature (K)
    dewpoint: the dewpoint temperature (K)

Returns:
    The equivalent potential temperature (K)
    )pbdoc"
    );
    
    m_therm.def(
        "thetae",
        [](
            const_prof_arr_t pres_arr,
            const_prof_arr_t tmpk_arr,
            const_prof_arr_t dwpk_arr
        ) {
            auto pres = pres_arr.view();
            auto tmpk = tmpk_arr.view();
            auto dwpk = dwpk_arr.view();

            float* thetae = new float[pres.shape(0)];

            for (size_t k = 0; k < pres.shape(0); ++k) {
                thetae[k] = sharp::thetae(pres(k), tmpk(k), dwpk(k));
            }
            
            nb::capsule owner(thetae,
                            [](void *p) noexcept { delete[] (float *)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                thetae, {tmpk.shape(0)}, owner);
        },
        nb::arg("pressure"), nb::arg("temperature"), nb::arg("dewpoint"),
        R"pbdoc(
Compute the equivalent potential temperature (K) given the 
air pressure (Pa), air temperature (K), and dewpoint temperature (K).
The equivalent potential temperature is computed as in Bolton 1980. 
It was found in Davies-Jones 2009 to be the most accurate non-iterative 
formulation of theta-e.

 Bolton 1980: https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2
 Davies-Jones 2009: https://doi.org/10.1175/2009MWR2774.1

Parameters:
    pressure: 1D NumPy array of ambient air pressure (Pa)
    temperature: 1D NumPy array of ambient air temperature (K)
    dewpoint: 1D NumPy array of dewpoint temperature (K)

Returns:
    1D NumPy array of equivalent potential temperature (K)
    )pbdoc"
    );

    m_therm.def(
        "lapse_rate",
        [](
            sharp::HeightLayer layer,
            const_prof_arr_t hght_arr,
            const_prof_arr_t tmpk_arr
        ) {
            return sharp::lapse_rate(layer, hght_arr.data(), tmpk_arr.data(), tmpk_arr.size());
        },
        nb::arg("layer_agl"), nb::arg("height"),
        nb::arg("temperature"),
        R"pbdoc(
Computes the lapse rate over a given HeightLayer (meters AGL).
This routine handles converting the height AGL to MSL by adding
the surface height value to the HeightLayer.

Parameters:
    layer_agl: a HeightLayer (meters AGL) to compute the lapse rate over
    height: 1D NumPy array of height values (meters)
    temperature: 1D NumPy array of temperature values (K)

Returns:
    The temperature lapse rate (K)
    )pbdoc"
    );

    m_therm.def(
        "lapse_rate",
        [](
            sharp::PressureLayer layer,
            const_prof_arr_t pres_arr,
            const_prof_arr_t hght_arr,
            const_prof_arr_t tmpk_arr
        ) {
            return sharp::lapse_rate(layer, pres_arr.data(), hght_arr.data(), tmpk_arr.data(), tmpk_arr.size());
        },
        nb::arg("layer"), nb::arg("pressure"), 
        nb::arg("height"), nb::arg("temperature"),
        R"pbdoc(
Computes the lapse rate over a given PressureLayer (Pa).

Parameters:
    layer: a PressureLayer (Pa) to compute the lapse rate over
    pressure: 1D NumPy array of pressure values (Pa)
    height: 1D NumPy array of height values (meters)
    temperature: 1D NumPy array of temperature values (K)

Returns:
    The temperature lapse rate (K)
    )pbdoc"
    );

    m_therm.def(
        "lapse_rate_max",
        [](
            sharp::HeightLayer layer,
            float depth, 
            const_prof_arr_t hght_arr, 
            const_prof_arr_t tmpk_arr
        ) {
            
            sharp::HeightLayer out_lyr = {0, 0};
            float max_lr = sharp::lapse_rate_max(
                layer, 
                depth, 
                hght_arr.data(), 
                tmpk_arr.data(), 
                tmpk_arr.size(), 
                &out_lyr
            );

            return nb::make_tuple(max_lr, out_lyr);
        },
        nb::arg("layer"), nb::arg("depth"),
        nb::arg("height"), nb::arg("temperature"),
        R"pbdoc(
Given a layer of the atmosphere (e.g. 2 - 6 km), find the maximum
lapse rate over the provided depth (e.g. 2 km) within that given layer. 
Returns the maximum lapse rate, as well as the layer it was found in. 

Parameters:
    layer: a HeightLayer to search for the maximum lapse lapse_rate (meters AGL) 
    depth: the depth used to compute a lapse rate within a layer (meters)
    height: 1D NumPy array of height values (meters)
    temperature: 1D NumPy array of temperature values (K)

Returns:
    A tuple containing the maximum lapse rate and the layer it was found in
    )pbdoc"
    );

    m_therm.def(
        "lapse_rate_max",
        [](
            sharp::PressureLayer layer,
            float depth, 
            const_prof_arr_t pres_arr,
            const_prof_arr_t hght_arr, 
            const_prof_arr_t tmpk_arr
        ) {
            
            sharp::PressureLayer out_lyr = {0, 0};
            float max_lr = sharp::lapse_rate_max(
                layer, 
                depth, 
                pres_arr.data(),
                hght_arr.data(), 
                tmpk_arr.data(), 
                tmpk_arr.size(), 
                &out_lyr
            );

            return nb::make_tuple(max_lr, out_lyr);
        },
        nb::arg("layer"), nb::arg("depth"),
        nb::arg("pressure"), nb::arg("height"), nb::arg("temperature"),
        R"pbdoc(
Given a layer of the atmosphere (e.g. 800 hPa - 500 hPa), find the maximum
lapse rate over the provided depth (e.g. 100 hPa) within that given layer. 
Returns the maximum lapse rate, as well as the layer it was found in. 

Parameters:
    layer: a PressureLayer to search for the maximum lapse lapse_rate (Pa) 
    depth: the depth used to compute a lapse rate within a layer (Pa)
    pressure: 1D NumPy array of pressure values (Pa)
    height: 1D NumPy array of height values (meters)
    temperature: 1D NumPy array of temperature values (K)

Returns:
    A tuple containing the maximum lapse rate and the layer it was found in
    )pbdoc"
    );

    m_therm.def(
        "buoyancy",
        [](float pcl_t, float env_t) {
            return sharp::buoyancy(pcl_t, env_t);
        },
        nb::arg("parcel_temperature"),
        nb::arg("environment_temperature"),
        R"pbdoc(
Compute buoyancy given parcel & environment temperatures.

Parameters:
    parcel_temperature: parcel temperature (K)
    environment_temperature: environment temperature (K)

Returns:
    Buoyancy (m/s^2)
    )pbdoc"
    );

    m_therm.def(
        "buoyancy",
        [](const_prof_arr_t pcl_tmpk, const_prof_arr_t env_tmpk) {

            float* buoy_arr = new float[pcl_tmpk.size()];

            sharp::buoyancy(pcl_tmpk.data(), env_tmpk.data(), buoy_arr, pcl_tmpk.size());

        nb::capsule owner(buoy_arr,
                        [](void *p) noexcept { delete[] (float *)p; });

        return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
            buoy_arr, {pcl_tmpk.size()}, owner);
        },
        nb::arg("parcel_temperature"), nb::arg("environment_temperature"),
        R"pbdoc(
Compute buoyancy given arrays of parcel & environment temperatures.

Parameters:
    parcel_temperature: 1D NumPy array of parcel temperatures (K)
    environment_temperature: 1D NumPy array of environment temperatures (K)

Returns:
    1D NumPy array of buoyancy values (m/s^2)
    )pbdoc"
    );
}

#endif
