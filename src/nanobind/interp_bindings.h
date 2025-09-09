#ifndef SHARPLIB_INTERP_BINDINGS
#define SHARPLIB_INTERP_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>

// clang-format on
#include <SHARPlib/interp.h>

#include "sharplib_types.h"

namespace nb = nanobind;

inline void make_interp_bindings(nb::module_ m) {
    nb::module_ m_interp = m.def_submodule(
        "interp",
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Interpolation Routines");

    m_interp.def(
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

Parameters
----------
hght_val : float 
    The coordinate height value to interpolate to (meters)
hght_arr : numpy.ndarray[dtype=float32]
    1D numpy array of height values to interpolate from (meters)
data_arr : numpy.ndarray[dtype=float32]
    float 1D numpy array of data values to interpolate from 

Returns
-------
float 
    An interpolated data value from data_arr corresponding to hght_val
        )pbdoc");

    m_interp.def(
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

Parameters
----------
pres_val : float 
    The coordinate pressure value to interpolate to (Pa)
pres_arr : nump.ndarray[dtype=float32] 
    1D numpy array of pressure values to interpolate from (Pa)
data_arr : numpy.ndarray[dtype=float32]
    1D numpy array of data values to interpolate from 

Returns
-------
float 
    An interpolated data value from data_arr corresponding to pres_val 
        )pbdoc");

    m_interp.def(
        "find_first_pressure",
        [](float data_val, const_prof_arr_t pres_arr,
           const_prof_arr_t data_arr) {
            return sharp::find_first_pressure(data_val, pres_arr.data(),
                                              data_arr.data(), pres_arr.size());
        },
        nb::arg("data_val"), nb::arg("pressure"), nb::arg("data_array"),
        R"pbdoc(
Conducts a bottom-up search for the first occurrence of a given value,
and interpolates in order to get the pressure level it occurs at.

Parameters
----------
data_val : float 
    The value to search for 
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values (Pa)
data_array : numpy.ndarray[dtype=float32]
    1D NumPy array of values to search over

Returns
-------
float
    The pressure level of first occurrence (Pa)
    )pbdoc");

    m_interp.def(
        "find_first_height",
        [](float data_val, const_prof_arr_t hght_arr,
           const_prof_arr_t data_arr) {
            return sharp::find_first_height(data_val, hght_arr.data(),
                                            data_arr.data(), hght_arr.size());
        },
        nb::arg("data_val"), nb::arg("height"), nb::arg("data_array"),
        R"pbdoc(
Conducts a bottom-up search for the first occurrence of a given value,
and interpolates in order to get the height level it occurs at.

Parameters
----------
data_val : float 
    The value to search for 
height : numpy.ndarray[dtype=float32]
    1D NumPy array of height values (meters)
data_array : numpy.ndarray[dtype=float32]
    1D NumPy array of values to search over

Returns
-------
float
    The height level of first occurrence (meters)
    )pbdoc");
}

#endif
