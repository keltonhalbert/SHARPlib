#ifndef SHARPLIB_WINDS_BINDINGS
#define SHARPLIB_WINDS_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <SHARPlib/winds.h>
#include <fmt/core.h>

#include <cstddef>

#include "binding_utils.h"
#include "sharplib_types.h"

namespace nb = nanobind;

template <typename Layer>
void bind_helicity(nb::module_& mod, const char* doc_template,
                   const char* layer_name) {
    const char *coord_name, *coord_unit, *layer_unit;
    if constexpr (Layer::coord == sharp::LayerCoordinate::pressure) {
        coord_name = "pressure";
        coord_unit = "Pa";
        layer_unit = "Pa";
    } else {
        coord_name = "height";
        coord_unit = "meters";
        layer_unit = "meters AGL";
    }

    std::string doc = fmt::format(doc_template, layer_name, coord_name,
                                  coord_unit, layer_unit);
    mod.def(
        "helicity",
        [](Layer layer, sharp::WindComponents storm_motion,
           const_prof_arr_t coord, const_prof_arr_t u_wind,
           const_prof_arr_t v_wind) {
            check_equal_sizes(coord, u_wind, v_wind);
            const std::size_t NZ = coord.size();
            return sharp::helicity(layer, storm_motion, coord.data(),
                                   u_wind.data(), v_wind.data(), NZ);
        },
        nb::arg("layer"), nb::arg("storm_motion"), nb::arg(coord_name),
        nb::arg("u_wind"), nb::arg("v_wind"), doc.c_str());
}

template <typename Layer>
void bind_wind_shear(nb::module_& mod, const char* doc_template,
                     const char* layer_name) {
    const char *coord_name, *coord_unit, *layer_unit;
    if constexpr (Layer::coord == sharp::LayerCoordinate::pressure) {
        coord_name = "pressure";
        coord_unit = "Pa";
        layer_unit = "Pa";
    } else {
        coord_name = "height";
        coord_unit = "meters";
        layer_unit = "meters AGL";
    }

    std::string doc = fmt::format(doc_template, layer_name, coord_name,
                                  coord_unit, layer_unit);
    mod.def(
        "wind_shear",
        [](Layer layer, const_prof_arr_t height, const_prof_arr_t u_wind,
           const_prof_arr_t v_wind) {
            check_equal_sizes(height, u_wind, v_wind);
            return sharp::wind_shear(layer, height.data(), u_wind.data(),
                                     v_wind.data(), height.size());
        },
        nb::arg("layer"), nb::arg(coord_name), nb::arg("u_wind"),
        nb::arg("v_wind"), doc.c_str());
}

template <typename Layer>
void bind_max_wind(nb::module_& mod, const char* doc_template,
                   const char* layer_name) {
    const char *coord_name, *coord_unit, *layer_unit;
    if constexpr (Layer::coord == sharp::LayerCoordinate::pressure) {
        coord_name = "pressure";
        coord_unit = "Pa";
        layer_unit = "Pa";
    } else {
        coord_name = "height";
        coord_unit = "meters";
        layer_unit = "meters AGL";
    }

    std::string doc = fmt::format(doc_template, layer_name, coord_name,
                                  coord_unit, layer_unit);
    mod.def(
        "max_wind",
        [](Layer layer, const_prof_arr_t coord, const_prof_arr_t u_wind,
           const_prof_arr_t v_wind) {
            check_equal_sizes(coord, u_wind, v_wind);
            const std::size_t NZ = u_wind.size();

            return sharp::max_wind(layer, coord.data(), u_wind.data(),
                                   v_wind.data(), NZ);
        },
        nb::arg("layer"), nb::arg(coord_name), nb::arg("u_wind"),
        nb::arg("v_wind"), doc.c_str());
}

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

    const char* helicity_template_doc =
        R"pbdoc(
Computes the Storm Relative Helicity (SRH) over a given layer using 
storm motion vector components stored in a WindComponents object.

This integration occurs over the given arrays, using interpolation 
for the top and bottom of the integration layer, and native data levels
in between. 

Parameters
----------
layer : {0}
    A layer to integrate over ({3})
storm_motion : nwsspc.sharp.calc.winds.WindComponents 
    A WindComponents object with the storm motion in m/s 
{1} : numpy.ndarray[dtype=float32] 
    1D NumPy array of {1} values ({2})
u_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of environment U-wind components (m/s)
v_wind : numpy.ndarray[dtype=float32]
    1D NumPy array of environment V-wind components (m/s)

Returns
-------
float
    Storm Relative Helicity (m^2/s^2)
    )pbdoc";
    bind_helicity<sharp::HeightLayer>(m_wind, helicity_template_doc,
                                      "nwsspc.sharp.calc.layer.HeightLayer");
    bind_helicity<sharp::PressureLayer>(
        m_wind, helicity_template_doc, "nwsspc.sharp.calc.layer.PressureLayer");

    const char* wind_shear_template_doc =
        R"pbdoc(
Computes the U and V components of the wind shear over a 
layer given the vertical sounding arrays of height, u_wind,
and v_wind.

Parameters
----------
layer : {0} 
    Layer for which to compute wind shear ({3})
{1} : numpy.ndarray[dtype=float32] 
    1D NumPy array of {1} values ({2})
u_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of U-wind components (m/s)
v_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of V-wind components (m/s)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
    WindComponents of U and V wind shear components (m/s)
    )pbdoc";
    bind_wind_shear<sharp::HeightLayer>(m_wind, wind_shear_template_doc,
                                        "nwsspc.sharp.calc.layer.HeightLayer");
    bind_wind_shear<sharp::PressureLayer>(
        m_wind, wind_shear_template_doc,
        "nwsspc.sharp.calc.layer.PressureLayer");

    m_wind.def(
        "mean_wind",
        [](sharp::PressureLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t u_wind, const_prof_arr_t v_wind,
           const bool weighted) {
            check_equal_sizes(pressure, u_wind, v_wind);
            return sharp::mean_wind(layer, pressure.data(), u_wind.data(),
                                    v_wind.data(), pressure.size(), weighted);
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("u_wind"),
        nb::arg("v_wind"), nb::arg("weighted") = false,
        R"pbdoc(
Computes the mean wind over the given PressureLayer and input profile
arrays of pressure, U-wind, and V-wind components. 

Parameters
----------
layer : nwsspc.sharp.calc.layer.PressureLayer 
    PressureLayer over which to compute mean
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure coordinate values (Pa)
u_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of U-wind component values (m/s)
v_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of V-wind component values (m/s)
weighted : bool 
    Boolean flag to compute pressure-weighted mean wind (default: False)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
    WindComponents of U and V mean wind components (m/s)
    )pbdoc");

    const char* max_wind_template_doc =
        R"pbdoc(
Finds the maximum wind speed over a given layer,
returning the vector components.

Parameters
----------
layer : {0}
    Layer over which to compute the max ({3})
{1} : numpy.ndarray[dtype=float32]
    1D NumPy array of {1} coordinate values ({2})
u_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of U-wind component values (m/s)
v_wind : numpy.ndarray[dtype=float32] 
    1D NumPy array of V-wind component values (m/s)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
)pbdoc";

    bind_max_wind<sharp::HeightLayer>(m_wind, max_wind_template_doc,
                                      "nwsspc.sharp.calc.layer.HeightLayer");
    bind_max_wind<sharp::PressureLayer>(
        m_wind, max_wind_template_doc, "nwsspc.sharp.calc.layer.PressureLayer");

    m_wind.def("vector_magnitude", &sharp::vector_magnitude, nb::arg("u_comp"),
               nb::arg("v_comp"),
               R"pbdoc(
Given the zonal (U) and meridional (V) wind components of a vector,
compute and return the magnitude (m/s) of the vector.

Parameters
----------
u_comp : float
    U-wind component (m/s)
v_comp : float
    V-wind component (m/s)

Returns
-------
float
    Wind speed (m/s)

    )pbdoc");

    m_wind.def(
        "vector_magnitude",
        [](const_prof_arr_t u_comp_arr, const_prof_arr_t v_comp_arr) {
            check_equal_sizes(u_comp_arr, v_comp_arr);
            const std::size_t NZ = u_comp_arr.size();
            auto u_comp = u_comp_arr.view();
            auto v_comp = v_comp_arr.view();

            return make_output_array(NZ, [&](float* out) {
                for (std::size_t idx = 0; idx < NZ; ++idx) {
                    out[idx] =
                        sharp::vector_magnitude(u_comp(idx), v_comp(idx));
                }
            });
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector, 
compute and return the magnitude (m/s) of the vector. 

Parameters
----------
u_comp : numpy.ndarray[dtype=float32] 
    1D NumPy array of U-wind component (m/s)
v_comp : numpy.ndarray[dtype=float32] 
    1D NumPy array of V-wind component (m/s)

Returns
-------
numpy.ndarray[dtype=float32]
    1D NumPy array of wind speed (m/s)
    )pbdoc");

    m_wind.def("vector_angle", &sharp::vector_angle, nb::arg("u_comp"),
               nb::arg("v_comp"),
               R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector,
compute and return the angle (from North) of the vector. 

Parameters
----------
u_comp : float 
    The U-wind component
v_comp : float 
    The V-wind component

Returns
-------
float
    Wind direction (degrees from North)
    )pbdoc");

    m_wind.def(
        "vector_angle",
        [](const_prof_arr_t u_comp_arr, const_prof_arr_t v_comp_arr) {
            check_equal_sizes(u_comp_arr, v_comp_arr);
            const std::size_t NZ = u_comp_arr.size();
            auto u_comp = u_comp_arr.view();
            auto v_comp = v_comp_arr.view();

            return make_output_array(NZ, [&](float* out) {
                for (std::size_t idx = 0; idx < NZ; ++idx) {
                    out[idx] = sharp::vector_angle(u_comp(idx), v_comp(idx));
                }
            });
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given the zonal (U) and meridional (V) components of a vector, 
compute and return the angle (from North) of the vector. 

Parameters
----------
u_comp : numpy.ndarray[dtype=float32] 
    1D NumPy array of U-wind component (m/s)
v_comp : numpy.ndarray[dtype=float32] 
    1D NumPy array of V-wind component (m/s)

Returns
-------
numpy.ndarray[dtype=float32]
    1D NumPy array of wind direction (degrees from North)
    )pbdoc");

    m_wind.def("u_component", &sharp::u_component, nb::arg("wind_speed"),
               nb::arg("wind_direction"),
               R"pbdoc(
Computes the zonal (U) wind component from a wind vector.

Parameters
----------
wind_speed : float 
    The vector speed (m/s)
wind_direction : float 
    The vector direction (degrees from North)

Returns
-------
float
    The U-wind component (m/s)
    )pbdoc");

    m_wind.def(
        "u_component",
        [](const_prof_arr_t wind_speed_arr,
           const_prof_arr_t wind_direction_arr) {
            check_equal_sizes(wind_speed_arr, wind_direction_arr);
            const std::size_t NZ = wind_speed_arr.size();
            auto wind_speed = wind_speed_arr.view();
            auto wind_direction = wind_direction_arr.view();

            return make_output_array(NZ, [&](float* out) {
                for (size_t k = 0; k < NZ; ++k) {
                    out[k] =
                        sharp::u_component(wind_speed(k), wind_direction(k));
                }
            });
        },
        nb::arg("wind_speed"), nb::arg("wind_direction"),
        R"pbdoc(
Computes the zonal (U) wind component from a wind vector.

Parameters
----------
wind_speed : numpy.ndarray[dtype=float32] 
    1D NumPy array of wind speeds (m/s)
wind_direction : numpy.ndarray[dtype=float32] 
    1D NumPy array of wind direction (degrees)

Returns
-------
numpy.ndarray[dtype=float32]
    1D NumPy array of U-wind component values (m/s)
    )pbdoc");

    m_wind.def("v_component", &sharp::v_component, nb::arg("wind_speed"),
               nb::arg("wind_direction"),
               R"pbdoc(
Computes the meridional (V) wind component from a wind vector.

Parameters
----------
wind_speed : float 
    Vector speed (m/s)
wind_direction : float 
    Vector angle (degrees from North)

Returns
-------
float
    The V-wind component
    )pbdoc");

    m_wind.def(
        "v_component",
        [](const_prof_arr_t wind_speed_arr,
           const_prof_arr_t wind_direction_arr) {
            check_equal_sizes(wind_speed_arr, wind_direction_arr);
            const std::size_t NZ = wind_speed_arr.size();
            auto wind_speed = wind_speed_arr.view();
            auto wind_direction = wind_direction_arr.view();

            return make_output_array(NZ, [&](float* out) {
                for (size_t k = 0; k < NZ; ++k) {
                    out[k] =
                        sharp::v_component(wind_speed(k), wind_direction(k));
                }
            });
        },
        nb::arg("wind_speed"), nb::arg("wind_direction"),
        R"pbdoc(
Computes the meridional (V) wind component from a wind vector.

Parameters
----------
wind_speed : numpy.ndarray[dtype=float32] 
    1D NumPy array of wind speeds (m/s)
wind_direction : numpy.ndarray[dtype=float32] 
    1D NumPy array of wind direction (degrees)

Returns
-------
numpy.ndarray[dtype=float32]
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
        
Parameters
----------
u_comp : float 
    The U-wind component (m/s)
v_comp : float 
    The V-wind component (m/s)

Returns
-------
nwsspc.sharp.calc.winds.WindVector
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
        
Parameters
----------
wind_comp : nwsspc.sharp.calc.winds.WindComponents 
    WindComponents (m/s)

Returns
-------
nwsspc.sharp.calc.winds.WindVector
    WindVector containing wind speed (m/s) and direction (degrees from North)
    )pbdoc");

    m_wind.def(
        "components_to_vector",
        [](const_prof_arr_t u_comp_arr,
           const_prof_arr_t v_comp_arr) -> std::tuple<out_arr_t, out_arr_t> {
            check_equal_sizes(u_comp_arr, v_comp_arr);
            const std::size_t NZ = u_comp_arr.size();
            auto u_comp = u_comp_arr.view();
            auto v_comp = v_comp_arr.view();
            auto wspd_buf = std::make_unique<float[]>(NZ);
            auto wdir_buf = std::make_unique<float[]>(NZ);

            for (size_t k = 0; k < NZ; ++k) {
                wspd_buf[k] = sharp::vector_magnitude(u_comp(k), v_comp(k));
                wdir_buf[k] = sharp::vector_angle(u_comp(k), v_comp(k));
            }

            float* wspd = wspd_buf.release();
            float* wdir = wdir_buf.release();
            nb::capsule owner_wspd(
                wspd, [](void* p) noexcept { delete[] (float*)p; });
            nb::capsule owner_wdir(
                wdir, [](void* p) noexcept { delete[] (float*)p; });

            auto wspd_arr = out_arr_t(wspd, {NZ}, owner_wspd);
            auto wdir_arr = out_arr_t(wdir, {NZ}, owner_wdir);

            return std::make_tuple(wspd_arr, wdir_arr);
        },
        nb::arg("u_comp"), nb::arg("v_comp"),
        R"pbdoc(
Given 1D NumPy arrays of zonal (U) and meridional (V) wind components,
compute the wind speed (m/s) and direction (degrees from North).

Parameters
----------
u_comp : numpy.ndarray[dtype=float32] 
    1D NumPy array of U-wind component values (m/s)
v_comp : numpy.ndarray[dtype=float32]
    1D NumPy array of V-wind component values (m/s)

Returns
-------
tuple[numpy.ndarray[dtype=float32], numpy.ndarray[dtype=float32]]
    wspd: 1D NumPy array of wind speeds (m/s)

    wdir: 1D NumPy array of wind directions (degrees from North)
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

Parameters
----------
wind_speed : float 
    The magnitude/speed of the vector (m/s)
wind_direction : float 
    The direction of the vector (degrees from North)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
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

Parameters
----------
wind_vector : nwsspc.sharp.calc.winds.WindVector 
    WindVector containing wind speed (m/s) and direction (degrees from North)

Returns
-------
nwsspc.sharp.calc.winds.WindComponents
    U and V wind components (m/s) in a WindComponents object
    )pbdoc");

    m_wind.def(
        "vector_to_components",
        [](const_prof_arr_t wspd_arr,
           const_prof_arr_t wdir_arr) -> std::tuple<out_arr_t, out_arr_t> {
            check_equal_sizes(wspd_arr, wdir_arr);
            const std::size_t NZ = wspd_arr.size();
            auto wspd = wspd_arr.view();
            auto wdir = wdir_arr.view();
            auto u_comp = std::make_unique<float[]>(NZ);
            auto v_comp = std::make_unique<float[]>(NZ);

            for (std::size_t k = 0; k < NZ; ++k) {
                u_comp[k] = sharp::u_component(wspd(k), wdir(k));
                v_comp[k] = sharp::v_component(wspd(k), wdir(k));
            }

            float* uwin = u_comp.release();
            float* vwin = v_comp.release();
            nb::capsule owner_uwin(
                uwin, [](void* p) noexcept { delete[] (float*)p; });
            nb::capsule owner_vwin(
                vwin, [](void* p) noexcept { delete[] (float*)p; });

            auto uwin_arr = out_arr_t(uwin, {NZ}, owner_uwin);
            auto vwin_arr = out_arr_t(vwin, {NZ}, owner_vwin);

            return std::make_tuple(uwin_arr, vwin_arr);
        },
        nb::arg("wspd"), nb::arg("wdir"),
        R"pbdoc(
Given 1D NumPy arrays of the wind speed (m/s) and direction (degrees from North),
compute and return 1D NumPy arrays of the zonal (U) and meridional (V) vector 
components.

Parameters
----------
wspd : numpy.ndarray[dtype=float32] 
    1D NumPy array of wind speeds (m/s)
wdir : numpy.ndarray[dtype=float32]
    1D NumPy array of wind directions (degrees from North)

Returns
-------
tuple[numpy.ndarray[dtype=float32], numpy.ndarray[dtype=float32]]
    uwin: 1D NumPy array of U-wind components

    vwin: 1D NumPy array of V-wind components
    )pbdoc");
}

#endif
