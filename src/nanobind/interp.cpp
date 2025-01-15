// clang-format off
#include <nanobind/nanobind.h>

// cland-format on
#include <SHARPlib/interp.h>
#include "sharplib_types.h"

namespace nb = nanobind;
// clang-format on
NB_MODULE(interp, m) {
    m.doc() =
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Interpolation Routines";

    m.def(
        "interp_height",
        [](const float hght_val, const_prof_arr_t hght_arr,
           const_prof_arr_t data_arr) -> float {
            if ((hght_arr.size() != data_arr.size())) {
                throw nb::buffer_error(
                    "hght_arr and data_arr must have the same size!");
            }
            return sharp::interp_height(hght_val, hght_arr.data(),
                                        data_arr.data(), hght_arr.size());
        },
        nb::arg("hght_val"), nb::arg("hght_arr"), nb::arg("data_arr"),
        R"pbdoc(
Interpolate a value from an array in height coordinates (meters).
The coordinate array (hght_arr) is assumed to be sorted (ascending).

Parameters:
    hght_val: The coordinate height value to interpolate to (meters)
    hght_arr: 1D numpy array of height values to interpolate from (meters)
    data_arr: 1D numpy array of data values to interpolate from 
Returns:
    An interpolated data value from data_arr corresponding to hght_val

        )pbdoc");

    m.def(
        "interp_pressure",
        [](const float pres_val, const_prof_arr_t pres_arr,
           const_prof_arr_t data_arr) -> float {
            if ((pres_arr.size() != data_arr.size())) {
                throw nb::buffer_error(
                    "hght_arr and data_arr must have the same size!");
            }
            return sharp::interp_pressure(pres_val, pres_arr.data(),
                                          data_arr.data(), pres_arr.size());
        },
        nb::arg("pres_val"), nb::arg("pres_arr"), nb::arg("data_arr"),
        R"pbdoc(
Interpolate a value from an array in pressure coordinates (Pa).
All pressure interpolation happens in log10 space.
The coordinate array (pres_arr) is assumed to be sorted (descending).

Parameters:
    pres_val: The coordinate pressure value to interpolate to (Pa)
    pres_arr: 1D numpy array of pressure values to interpolate from (Pa)
    data_arr: 1D numpy array of data values to interpolate from 
Returns:
    An interpolated data value from data_arr corresponding to pres_val 

        )pbdoc");
}
