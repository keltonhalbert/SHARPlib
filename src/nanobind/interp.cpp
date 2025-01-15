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
        [](const float height_val, const_prof_arr_t hght_arr,
           const_prof_arr_t data_arr) -> float {
            return sharp::interp_height(height_val, hght_arr.data(),
                                        data_arr.data(), hght_arr.size());
        },
        nb::arg("height_val"), nb::arg("hght_arr"), nb::arg("data_arr"),
        R"pbdoc(
            Interpolate a value from an array in height coordinates (meters).

            Parameters:
                height_val: The coordinate height value to interpolate to (meters)
                hght_arr: 1D numpy array of height values to interpolate from (meters)
                data_arr: 1D numpy array of data values to interpolate from 
            Returns:
                An interpolated data value from data_arr corresponding to height_val

        )pbdoc");
}
