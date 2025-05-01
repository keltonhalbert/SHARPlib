#ifndef SHARPLIB_WINDS_BINDINGS
#define SHARPLIB_WINDS_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

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
        .def(nb::init<float, float>(), nb::arg("speed"), nb::arg("direction"))
        .def_rw("speed", &sharp::WindVector::speed, "Wind Speed (m/s)")
        .def_rw("direction", &sharp::WindVector::direction,
                "Wind Direction (degrees from North)");

    nb::class_<sharp::WindComponents>(m_wind, "WindComponents")
        .def(nb::init<>())
        .def(nb::init<float, float>(), nb::arg("u_comp"), nb::arg("v_comp"))
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

    m_wind.def("vector_magnitude", &sharp::vector_magnitude, nb::arg("u_comp"),
               nb::arg("v_comp"),
               R"pbdoc(
Given the zonal (U) and meridional (V) wind components of a vector,
compute and return the magnitude (m/s) of the vector.

Parameters:
    u_comp: U-wind component (m/s)
    v_comp: V-wind component (m/s)

Returns: 
    Wind speed (m/s)

    )pbdoc");

    m_wind.def(
        "vector_magnitude",
        [](float u_comp, float v_comp) {
            return sharp::vector_magnitude(u_comp, v_comp);
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector, 
compute and return the magnitude (m/s) of the vector. 

Parameters:
    u_comp: U-wind component (m/s)
    v_comp: V-wind component (m/s)

Returns:
    Wind speed (m/s)
    )pbdoc");

    m_wind.def(
        "vector_magnitude",
        [](const_prof_arr_t u_comp_arr, const_prof_arr_t v_comp_arr) {
            auto u_comp = u_comp_arr.view();
            auto v_comp = v_comp_arr.view();
            auto len = u_comp.shape(0);

            float* wspd_arr = new float[u_comp.shape(0)];

            for (int idx = 0; idx < len; ++idx) {
                wspd_arr[idx] =
                    sharp::vector_magnitude(u_comp(idx), v_comp(idx));
            }

            nb::capsule owner(wspd_arr,
                              [](void* p) noexcept { delete[] (float*)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wspd_arr, {u_comp.shape(0)}, owner);
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector, 
compute and return the magnitude (m/s) of the vector. 

Parameters:
    u_comp: 1D NumPy array of U-wind component (m/s)
    v_comp: 1D NumPy array of V-wind component (m/s)

Returns:
    1D NumPy array of wind speed (m/s)
    )pbdoc");

    m_wind.def("vector_angle", &sharp::vector_angle, nb::arg("u_comp"),
               nb::arg("v_comp"),
               R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector,
compute and return the angle (from North) of the vector. 

Parameters:
    u_comp: The U-wind component
    v_comp: The V-wind component

Returns:
    Wind direvtion (degrees from North)
    )pbdoc");

    m_wind.def(
        "vector_angle",
        [](const_prof_arr_t u_comp_arr, const_prof_arr_t v_comp_arr) {
            auto u_comp = u_comp_arr.view();
            auto v_comp = v_comp_arr.view();
            auto len = u_comp.shape(0);

            float* wdir_arr = new float[u_comp.shape(0)];

            for (int idx = 0; idx < len; ++idx) {
                wdir_arr[idx] = sharp::vector_angle(u_comp(idx), v_comp(idx));
            }

            nb::capsule owner(wdir_arr,
                              [](void* p) noexcept { delete[] (float*)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wdir_arr, {u_comp.shape(0)}, owner);
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector, 
compute and return the angle (from North) of the vector. 

Parameters:
    u_comp: 1D NumPy array of U-wind component (m/s)
    v_comp: 1D NumPy array of V-wind component (m/s)

Returns:
    1D NumPy array of wind direction (degrees from North)
    )pbdoc");

    m_wind.def("u_component", &sharp::u_component, nb::arg("wind_speed"),
               nb::arg("wind_direction"),
               R"pbdoc(
Computes the zonal (U) wind component from a wind vector.

Parameters: 
    wind_speed: The vector speed (m/s)
    wind_direction: The vector direction (degrees from North)

Returns:
    The U-wind component (m/s)
    )pbdoc");

    m_wind.def(
        "u_component",
        [](const_prof_arr_t wind_speed_arr,
           const_prof_arr_t wind_direction_arr) {
            auto wind_speed = wind_speed_arr.view();
            auto wind_direction = wind_direction_arr.view();
            float* uwin_arr = new float[wind_speed.shape(0)];

            for (size_t k = 0; k < wind_speed.shape(0); ++k) {
                uwin_arr[k] =
                    sharp::u_component(wind_speed(k), wind_direction(k));
            }

            nb::capsule owner(uwin_arr,
                              [](void* p) noexcept { delete[] (float*)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                uwin_arr, {wind_speed.shape(0)}, owner);
        },
        nb::arg("wind_speed"), nb::arg("wind_direction"),
        R"pbdoc(
Computes the zonal (U) wind component from a wind vector.

Parameters:
    wind_speed: 1D NumPy array of wind speeds (m/s)
    wind_direction: 1D NumPy array of wind direction (degrees)

Returns:
    1D NumPy array of U-wind component values (m/s)
    )pbdoc");

    m_wind.def("v_component", &sharp::v_component, nb::arg("wind_speed"),
               nb::arg("wind_direction"),
               R"pbdoc(
Computes the meridional (V) wind component from a wind vector.

Parameters:
    wind_speed: Vector speed (m/s)
    wind_direction: Vector angle (degrees from North)

Returns:
    The V-wind component
    )pbdoc");

    m_wind.def(
        "v_component",
        [](const_prof_arr_t wind_speed_arr,
           const_prof_arr_t wind_direction_arr) {
            auto wind_speed = wind_speed_arr.view();
            auto wind_direction = wind_direction_arr.view();
            float* vwin_arr = new float[wind_speed.shape(0)];

            for (size_t k = 0; k < wind_speed.shape(0); ++k) {
                vwin_arr[k] =
                    sharp::v_component(wind_speed(k), wind_direction(k));
            }

            nb::capsule owner(vwin_arr,
                              [](void* p) noexcept { delete[] (float*)p; });

            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                vwin_arr, {wind_speed.shape(0)}, owner);
        },
        nb::arg("wind_speed"), nb::arg("wind_direction"),
        R"pbdoc(
Computes the meridional (V) wind component from a wind vector.

Parameters:
    wind_speed: 1D NumPy array of wind speeds (m/s)
    wind_direction: 1D NumPy array of wind direction (degrees)

Returns:
    1D NumPy array of V-wind component values (m/s)
    )pbdoc");

    m_wind.def(
        "components_to_vector",
        [](float u_comp, float v_comp) {
            return sharp::components_to_vector(u_comp, v_comp);
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector, 
compute and return the wind speed (m/s) and direction 
(degrees from North) as a WindVector type.
        
Parameters:
    u_comp: The U-wind component (m/s)
    v_comp: The V-wind component (m/s)

Returns:
    WindVector containing wind speed (m/s) and direction (degrees from North)
    )pbdoc");

    m_wind.def(
        "components_to_vector",
        [](sharp::WindComponents& cmp) {
            return sharp::components_to_vector(cmp);
        },
        nb::arg("wind_comp"),
        R"pbdoc(
Given the components of a vector via a WindComponents object (m/s),
compute and return the wind speed (m/s) and wind direction (degrees from North)
as a WindVector type.
        
Parameters:
    wind_comp: WindComponents (m/s)

Returns:
    WindVector containing wind speed (m/s) and direction (degrees from North)
    )pbdoc");

    m_wind.def(
        "components_to_vector",
        [](const_prof_arr_t u_comp_arr, const_prof_arr_t v_comp_arr) {
            auto u_comp = u_comp_arr.view();
            auto v_comp = v_comp_arr.view();

            float* wspd = new float[u_comp.shape(0)];
            float* wdir = new float[v_comp.shape(0)];

            for (size_t k = 0; k <= u_comp.shape(0); ++k) {
                wspd[k] = sharp::vector_magnitude(u_comp(k), v_comp(k));
                wdir[k] = sharp::vector_angle(u_comp(k), v_comp(k));
            }

            nb::capsule owner_wspd(
                wspd, [](void* p) noexcept { delete[] (float*)p; });
            nb::capsule owner_wdir(
                wdir, [](void* p) noexcept { delete[] (float*)p; });

            auto wspd_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wspd, {u_comp.shape(0)}, owner_wspd);
            auto wdir_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                wdir, {u_comp.shape(0)}, owner_wdir);

            return nb::make_tuple(wspd_arr, wdir_arr);
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given 1D NumPy arrays of zonal (U) and meridional (V) wind components,
compute the wind speed (m/s) and direction (degrees from North).

Parameters:
    u_comp: 1D NumPy array of U-wind component values (m/s)
    v_comp: 1D NumPy array of V-wind component values (m/s)

Returns:
    wspd: 1D NumPy array of wind speeds (m/s)
    wdir: 1D NumPy array of wind directiond (degrees from North)
    )pbdoc");

    m_wind.def(
        "vector_to_components",
        [](float wind_speed, float wind_direction) {
            return sharp::vector_to_components(wind_speed, wind_direction);
        },
        nb::arg("wind_speed"), nb::arg("wind_direction"),
        R"pbdoc(
Given the wind speed (m/s) and wind direction (degrees from North),
compute and return the zonal and meridional vector components as WindComponents.

Parameters:
    wind_speed: The magnitude/speed of the vector (m/s)
    wind_direction: The direction of the vector (degrees from North)

Returns:
    The U-wind and V-wind components (m/s) as WindComponents
    )pbdoc");

    m_wind.def(
        "vector_to_components",
        [](sharp::WindVector& vec) { return sharp::vector_to_components(vec); },
        nb::arg("wind_vector"),
        R"pbdoc(
Given the wind speed (m/s) and direction (degrees from North),
compute and return the zonal (U) and meridional (V) vector 
components as WindComponents.

Parameters:
    wind_vector: WindComponents containing wind wpeed (m/s) and direction (degrees from North)

Returns:
    U and V wind components (m/s) in a WindComponents object
    )pbdoc");

    m_wind.def(
        "vector_to_components",
        [](const_prof_arr_t wspd_arr, const_prof_arr_t wdir_arr) {
            auto wspd = wspd_arr.view();
            auto wdir = wdir_arr.view();

            float* u_comp = new float[wspd.shape(0)];
            float* v_comp = new float[wspd.shape(0)];

            for (size_t k = 0; k <= wspd.shape(0); ++k) {
                u_comp[k] = sharp::u_component(wspd(k), wdir(k));
                v_comp[k] = sharp::v_component(wspd(k), wdir(k));
            }

            nb::capsule owner_uwin(
                u_comp, [](void* p) noexcept { delete[] (float*)p; });
            nb::capsule owner_vwin(
                v_comp, [](void* p) noexcept { delete[] (float*)p; });

            auto uwin_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                u_comp, {wspd.shape(0)}, owner_uwin);
            auto vwin_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                v_comp, {wdir.shape(0)}, owner_vwin);

            return nb::make_tuple(uwin_arr, vwin_arr);
        },
        nb::arg("wspd"), nb::arg("wdir"),
        R"pbdoc(
Given 1D NumPy arrays of the wind speed (m/s) and direction (degrees from North),
compute and return 1D NumPy arays of the zonal (U) and meridional (V) vector 
componenWindComponents.

Parameters:
    wspd: 1D NumPy array of wind speeds (m/s)
    wdir: 1D NumPy array of wind directions (degrees from North)

Returns:
    uwin: 1D NumPy array of U-wind components
    vwin: 1D NumPy array of V-wind components
    )pbdoc");
}

#endif
