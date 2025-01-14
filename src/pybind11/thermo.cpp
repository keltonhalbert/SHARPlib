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

    py::enum_<sharp::adiabat>(m, "adiabat")
        .value("pseudo_liq", sharp::adiabat::pseudo_liq,
               "Use a pseudoadiabatic parcel with liquid only processes")
        .value("adiab_liq", sharp::adiabat::adiab_liq,
               "Use an adiabatic parcel with liquid only processes")
        .value("pseudo_ice", sharp::adiabat::pseudo_ice,
               "Use a pseudoadiabatic parcel with liquid and ice processes")
        .value("adiab_ice", sharp::adiabat::adiab_ice,
               "Use an adiabatic parcel with liquid and ice processes");

    m.def("wobf", &sharp::wobf, py::arg("temperature"),
          R"pbdoc(
           Computes the difference between the wet-bulb potential temperatures for sarturated and           dry air, given a temperature.

        The Wobus Function (wobf) is defined as the difference between the wet-bulb potential 
        temperature for saturated air (WBPTS) and the wet-bulb potential temperature for dry air         (WBPTD) given an input temperature in Kelvin.

        WOBF(T) := WBPTS - WBPTD

        Although WBPTS and WBPTD are functions of both pressure and temperature, it is assumed 
        their difference is a function of temperature only. The difference is also proportional 
        to the heat imparted to an air parcel. 

        This function uses a polynomial approximation to the Wobus function, fitted to the 
        values in Table 78 of PP.319-322 of the Smithsonian Meteorologial Tables by Roland 
        List (6th Revised Edition). Herman Wobus, a mathematician for the Navy Weather 
        Research Facility in Norfolk, VA, computed these coefficients a very long time ago, 
        as he was retired at the time of the documentation found on this routine written in 1981.

        It has been shown by Robert Davies-Jones (2007) that the Wobus function has a slight 
        dependence on pressure, which results in errors of up to 1.2 Kelvin in the temperature 
        of a lifted parcel. 

        Parameters:
            temperature: Air temperature in Kelvin

        Returns:
            The value of the Wobus function in Kelvin 

          )pbdoc");

    m.def("lcl_temperature", &sharp::lcl_temperature,
          "Compute the Lifted Condensation Level (LCL) temperature (K) using "
          "the air temperature (K) and dewpoint temperature(K)");
    m.def("lcl_temperature", py::vectorize(sharp::lcl_temperature),
          "Compute the Lifted Condensation Level (LCL) temperature (K) using "
          "the air temperature (K) and dewpoint temperature(K)");
}
