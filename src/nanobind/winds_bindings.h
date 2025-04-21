#ifndef SHARPLIB_WINDS_BINDINGS
#define SHARPLIB_WINDS_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>

// clang-format on 
#include <SHARPlib/layer.h>
#include <SHARPlib/winds.h>

#include "sharplib_types.h"

// clang-format on
inline void make_winds_bindings(nb::module_ m) {
    nb::module_ m_wind =
        m.def_submodule("winds",
                        "Sounding and Hodograph Analysis and Research Program "
                        "Library (SHARPlib) :: Kinematic Routines");

    nb::class_<sharp::WindVector>(m_wind, "WindVector")
        .def(nb::init<>())
        .def_rw("speed", &sharp::WindVector::speed, "Wind Speed (m/s)")
        .def_rw("direction", &sharp::WindVector::direction,
                "Wind Direction (degrees from North)");

    nb::class_<sharp::WindComponents>(m_wind, "WindComponents")
        .def(nb::init<>())
        .def_rw("u", &sharp::WindComponents::u, "U wind component (m/s)")
        .def_rw("v", &sharp::WindComponents::v, "V wind component (m/s)");

    m_wind.def(
        "helicity",
        [](sharp::HeightLayer& layer, sharp::WindComponents storm_motion,
           const_prof_arr_t height, const_prof_arr_t u_wind,
           const_prof_arr_t v_wind) {
            return sharp::helicity(layer, storm_motion, height.data(),
                                   u_wind.data(), v_wind.data(), height.size());
        },
        nb::arg("layer"), nb::arg("storm_motion"), nb::arg("height"),
        nb::arg("u_wind"), nb::arg("v_wind"),
        R"pbdoc(
Computes the Storm Relative Helicity (SRH) over a given layer using 
storm motion vector components stored in a WindComponents object.

This integration occurs over the given arrays, using interpolation 
for the top and bottom of the integration layer, and native data levels
in between. 

Parameters:
    layer: A HeightLayer to integrate over (meters AGL)
    storm_motion: A WindComponents object with the storm motion in m/s 
    height: 1D NumPy array of height values (meters)
    u_wind: 1D NumPy array of environment U-wind components (m/s)
    v_wind: 1D NumPy array of environment V-wind components (m/s)

Returns:
    Storm Relative Helicity (m^2/s^2)
    )pbdoc");

    m_wind.def(
        "helicity",
        [](sharp::PressureLayer& layer, sharp::WindComponents storm_motion,
           const_prof_arr_t pressure, const_prof_arr_t u_wind,
           const_prof_arr_t v_wind) {
            return sharp::helicity(layer, storm_motion, pressure.data(),
                                   u_wind.data(), v_wind.data(),
                                   pressure.size());
        },
        nb::arg("layer"), nb::arg("storm_motion"), nb::arg("pressure"),
        nb::arg("u_wind"), nb::arg("v_wind"),
        R"pbdoc(
Computes the Storm Relative Helicity (SRH) over a given layer using 
storm motion vector components stored in a WindComponents object.

This integration occurs over the given arrays, using interpolation 
for the top and bottom of the integration layer, and native data levels
in between. 

Parameters:
    layer: A PressureLayer to integrate over (Pa)
    storm_motion: A WindComponents object with the storm motion in m/s 
    pressure: 1D NumPy array of pressure values (Pa)
    u_wind: 1D NumPy array of environment U-wind components (m/s)
    v_wind: 1D NumPy array of environment V-wind components (m/s)

Returns:
    Storm Relative Helicity (m^2/s^2)
    )pbdoc");

    m_wind.def(
        "wind_shear",
        [](sharp::HeightLayer layer, const_prof_arr_t height,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind) {
            return sharp::wind_shear(layer, height.data(), u_wind.data(),
                                     v_wind.data(), height.size());
        },
        nb::arg("layer"), nb::arg("height"), nb::arg("u_wind"),
        nb::arg("v_wind"),
        R"pbdoc(
Computes the U and V components of the wind shear over a 
layer given the vertical sounding arrays of height, u_wind,
and v_wind.

Parameters:
    layer: HeightLayer for which to compute wind shear 
    height: 1D NumPy array of height values (meters)
    u_wind: 1D NumPy array of U-wind components (m/s)
    v_wind: 1D NumPy array of V-wind components (m/s)

Returns:
    WindComponents of U and V wind shear components (m/s)
    )pbdoc");

    m_wind.def(
        "wind_shear",
        [](sharp::PressureLayer layer, const_prof_arr_t pressure,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind) {
            return sharp::wind_shear(layer, pressure.data(), u_wind.data(),
                                     v_wind.data(), pressure.size());
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("u_wind"),
        nb::arg("v_wind"),
        R"pbdoc(
Computes the U and V components of the wind shear over a 
layer given the vertical sounding arrays of height, u_wind,
and v_wind.

Parameters:
    layer: PressureLayer for which to compute wind shear 
    pressure: 1D NumPy array of pressure values (Pa)
    u_wind: 1D NumPy array of U-wind components (m/s)
    v_wind: 1D NumPy array of V-wind components (m/s)

Returns:
    WindComponents of U and V wind shear components (m/s)
    )pbdoc");

    m_wind.def(
        "mean_wind",
        [](sharp::PressureLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t u_wind, const prof_arr_t v_wind,
           const bool weighted) {
            return sharp::mean_wind(layer, pressure.data(), u_wind.data(),
                                    v_wind.data(), pressure.size(), weighted);
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("u_wind"),
        nb::arg("v_wind"), nb::arg("weighted") = false,
        R"pbdoc(
Computes the mean wind over the given PressureLayer and input profile
arrays of pressure, U-wind, and V-wind components. 

Parameters:
    layer: PressureLayer over which to compute mean
    pressure: 1D NumPy array of pressure coordinate values (Pa)
    u_wind: 1D NumPy array of U-wind component values (m/s)
    v_wind: 1D NumPy array of V-wind component values (m/s)
    weighted: Boolean flag to compute pressure-weighted mean wind (default: False)

Returns:
    WindComponents of U anf V mean wind components (m/s)
    )pbdoc");
}

#endif
