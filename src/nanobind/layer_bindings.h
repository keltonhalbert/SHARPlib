#ifndef SHARPLIB_LAYER_BINDINGS
#define SHARPLIB_LAYER_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/string.h>
#include <fmt/core.h>

// clang-format on
#include <SHARPlib/layer.h>

#include "sharplib_types.h"

namespace nb = nanobind;

inline void make_layer_bindings(nb::module_ m) {
    nb::module_ m_layer = m.def_submodule(
        "layer",
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Atmospheric Layer Routines");

    // Bind the constructors, named fields, and default arguments
    nb::class_<sharp::HeightLayer>(m_layer, "HeightLayer")
        .def(nb::init<float, float, float>(), nb::arg("bottom"), nb::arg("top"),
             nb::arg("delta") = 100.0f)
        .def_rw("bottom", &sharp::HeightLayer::bottom,
                "The bottom of the HeightLayer (meters)")
        .def_rw("top", &sharp::HeightLayer::top,
                "The top of the HeightLayer (meters)")
        .def_rw(
            "delta", &sharp::HeightLayer::delta,
            "The HeightLayer delta (increment) to use if iterating (meters).")
        .def("__repr__", [](const sharp::HeightLayer& layer) -> std::string {
            return fmt::format("HeightLayer: [{0} meters, {1} meters]",
                               layer.bottom, layer.top);
        });

    nb::class_<sharp::PressureLayer>(m_layer, "PressureLayer")
        .def(nb::init<float, float, float>(), nb::arg("bottom"), nb::arg("top"),
             nb::arg("delta") = -1000.0f)
        .def_rw("bottom", &sharp::PressureLayer::bottom,
                "The bottom of the PressureLayer (Pa)")
        .def_rw("top", &sharp::PressureLayer::top,
                "The top of the PressureLayer (Pa)")
        .def_rw("delta", &sharp::PressureLayer::delta,
                "The PressureLayer delta (increment) to use if iterating (Pa)")
        .def("__repr__", [](const sharp::PressureLayer& layer) -> std::string {
            return fmt::format("PressureLayer: [{0} Pa, {1} Pa]", layer.bottom,
                               layer.top);
        });

    nb::class_<sharp::LayerIndex>(m_layer, "LayerIndex")
        .def(nb::init<std::ptrdiff_t, std::ptrdiff_t>(), nb::arg("kbot"),
             nb::arg("ktop"))
        .def_rw("kbot", &sharp::LayerIndex::kbot,
                "The bottom index of a layer on a coordinate array "
                "(pressure or height).")
        .def_rw("ktop", &sharp::LayerIndex::ktop,
                "The top index of a layer on a coordinate array "
                "(pressure or height).")
        .def("__repr__", [](const sharp::LayerIndex& layer) -> std::string {
            return fmt::format("LayerIndex: [{0}, {1}]", layer.kbot,
                               layer.ktop);
        });

    m_layer.def(
        "get_layer_index",
        [](sharp::PressureLayer& layer,
           const_prof_arr_t pressure_array) -> sharp::LayerIndex {
            return sharp::get_layer_index(layer, pressure_array.data(),
                                          pressure_array.size());
        },
        nb::arg("layer"), nb::arg("pressure"),
        R"pbdoc(
Finds the interior array indices encapsulated by the given PressureLayer. 
More specifically, the returned LayerIndex excludes the exact top and bottom 
values (in mathematical notation, [bottom, top]). This behavior is due to the 
fact many algorithms use interpolation to get the exact values of arbitrary
top/bottom locations within a profile. 

If layer.bottom or layer.top are out of bounds, this function 
will truncate the layer to the coordinate range of data provided 
by coord[] in an attempt to gracefully continue and produce a result. 
This will modify the value of the given PressureLayer.

Parameters:
    layer: PressureLayer
    pressure: 1D NumPy array of pressures

Returns:
    A LayerIndex with {kbot, ktop}.
            )pbdoc");

    m_layer.def(
        "get_layer_index",
        [](sharp::HeightLayer& layer,
           const_prof_arr_t height_array) -> sharp::LayerIndex {
            return sharp::get_layer_index(layer, height_array.data(),
                                          height_array.size());
        },
        nb::arg("layer"), nb::arg("height"),
        R"pbdoc(
Finds the interior array indices encapsulated by the given HeightLayer. 
More specifically, the returned LayerIndex excludes the exact top and bottom 
values (in mathematical notation, [bottom, top]). This behavior is due to the 
fact many algorithms use interpolation to get the exact values of arbitrary
top/bottom locations within a profile. 

If layer.bottom or layer.top are out of bounds, this function 
will truncate the layer to the coordinate range of data provided 
by coord[] in an attempt to gracefully continue and produce a result. 
This will modify the value of the given HeightLayer.

Parameters:
    layer: HeightLayer
    height: 1D NumPy array of heights

Returns:
    A LayerIndex with {kbot, ktop}.
        )pbdoc");

    m_layer.def(
        "height_layer_to_pressure",
        [](sharp::HeightLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t height, const bool isAGL) -> sharp::PressureLayer {
            return sharp::height_layer_to_pressure(
                layer, pressure.data(), height.data(), height.size(), isAGL);
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("height"),
        nb::arg("isAGL") = false,
        R"pbdoc(
Converts a HeightLayer to a PressureLayer via interpolation, with the
optional argument to convert the HeightLayer from meters AGL to MSL by adding
the station height to the HeightLayer. If for some strange reason 
you provide a HeightLayer that is out of the bounds of height[], then
the bottom and top of the output layer will be set to MISSING.

Parameters:
    layer: HeightLayer
    pressure: 1D NumPy array of pressure
    height: 1D NumPy array of heights
    isAGL: Whether or not the station height needs to be added for interpolation (default: False) 

    )pbdoc");

    m_layer.def(
        "pressure_layer_to_height",
        [](sharp::PressureLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t height, const bool toAGL) -> sharp::HeightLayer {
            return sharp::pressure_layer_to_height(
                layer, pressure.data(), height.data(), height.size(), toAGL);
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("height"),
        nb::arg("toAGL") = false,
        R"pbdoc(
Converts a PressureLayer to a HeightLayer via interpolation, with the 
optional argument to convert the HeightLayer to meters AGL by subtracting off 
the station height from the returned HeightLayer. If for some strange reason
you provide a PressureLayer that is out of the bounds of pressure[], then 
the bottom and top of the output layer will be set to MISSING. 

Parameters:
    layer: PressureLayer
    pressure: 1D NumPy array of pressure 
    height: 1D NumPy array of heights
    toAGL: Whether or not to subtract the station height from the HeightLayer (default: False)

    )pbdoc");

    m_layer.def(
        "layer_min",
        [](sharp::HeightLayer& layer, const_prof_arr_t height,
           const_prof_arr_t data) {
            float lvl_of_min;
            float min = sharp::layer_min(layer, height.data(), data.data(),
                                         data.size(), &lvl_of_min);

            return std::make_tuple(min, lvl_of_min);
        },
        nb::arg("layer"), nb::arg("height"), nb::arg("data"),
        R"pbdoc(
Returns the minimum value of the data array within the given HeightLayer. The 
function bounds checks the layer by calling get_layer_index. 

Parameters:
    layer: HeightLayer 
    height: 1D NumPy array of heights 
    data: 1D array of data 

Returns:
    tuple: (min_value, level_of_min)

    )pbdoc");

    m_layer.def(
        "layer_min",
        [](sharp::PressureLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t data) {
            float lvl_of_min;
            float min = sharp::layer_min(layer, pressure.data(), data.data(),
                                         data.size(), &lvl_of_min);

            return std::make_tuple(min, lvl_of_min);
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("data"),
        R"pbdoc(
Returns the minimum value of the data array within the given PressureLayer. The 
function bounds checks the layer by calling get_layer_index. 

Parameters:
    layer: PressureLayer 
    pressure: 1D NumPy array of pressure 
    data: 1D array of data 

Returns:
    tuple: (min_value, level_of_min)
    )pbdoc");

    m_layer.def(
        "layer_max",
        [](sharp::HeightLayer& layer, const_prof_arr_t height,
           const_prof_arr_t data) {
            float lvl_of_max;
            float max = sharp::layer_max(layer, height.data(), data.data(),
                                         data.size(), &lvl_of_max);

            return std::make_tuple(max, lvl_of_max);
        },
        nb::arg("layer"), nb::arg("height"), nb::arg("data"),
        R"pbdoc(
Returns the maximum value observed within the given data array over
the given HeightLayer. The function bounds checks the layer by calling 
get_layer_index. 

Parameters:
    layer: HeightLayer 
    height: 1D NumPy array of height values 
    data: 1D NumPy array of data values

Returns:
    tuple: (max_value, level_of_max)

    )pbdoc");

    m_layer.def(
        "layer_max",
        [](sharp::PressureLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t data) {
            float lvl_of_max;
            float max = sharp::layer_max(layer, pressure.data(), data.data(),
                                         data.size(), &lvl_of_max);

            return std::make_tuple(max, lvl_of_max);
        },
        nb::arg("layer"), nb::arg("pressure"), nb::arg("data"),
        R"pbdoc(
Returns the maximum value observed within the given data array over the given 
PressureLayer. The function bounds checks the layer by calling get_layer_index. 

Parameters:
    layer: PressureLayer 
    pressure: 1D NumPy array of pressure values 
    data: 1D NumPy array of data values

Returns:
    tuple: (max_value, level_of_max)

    )pbdoc");

    m_layer.def(
        "layer_mean",
        [](sharp::PressureLayer layer, const_prof_arr_t pressure,
           const_prof_arr_t data) {
            return sharp::layer_mean(layer, pressure.data(), data.data(),
                                     data.size());
        },
        nb::arg("PressureLayer"), nb::arg("pressure"), nb::arg("data"),
        R"pbdoc(
Computes the pressure-weighted mean value of a field over 
a given PressureLayer. 

Parameters:
    layer: PressureLayer 
    pressure: 1D NumPy array of pressure values 
    data: 1D NumPy array of data values

    )pbdoc");

    m_layer.def(
        "layer_mean",
        [](sharp::HeightLayer layer, const_prof_arr_t height,
           const_prof_arr_t pressure, const_prof_arr_t data, const bool isAGL) {
            return sharp::layer_mean(layer, height.data(), pressure.data(),
                                     data.data(), data.size(), isAGL);
        },
        nb::arg("HeightLayer"), nb::arg("height"), nb::arg("pressure"),
        nb::arg("data"), nb::arg("isAGL") = false,
        R"pbdoc(
Computes the pressure-weighted mean value of a field over 
a given HeightLayer. 

Parameters:
    layer: HeightLayer 
    height: 1D NumPy array of height values
    pressure: 1D NumPy array of pressure values 
    data: 1D NumPy array of data values
    isAGL: Whether or not the surface station height should be added to the HeightLayer (default: False)

    )pbdoc");

    m_layer.def(
        "integrate_layer_trapz",
        [](sharp::HeightLayer layer, const_prof_arr_t data,
           const_prof_arr_t height, const int integ_sign, const bool weighted) {
            return sharp::integrate_layer_trapz(layer, data.data(),
                                                height.data(), data.size(),
                                                integ_sign, weighted);
        },
        nb::arg("layer"), nb::arg("data"), nb::arg("height"),
        nb::arg("integ_sign") = 0, nb::arg("weighted") = false,
        R"pbdoc(
Returns a trapezoidal integration of the given data array over 
the given HeightLayer. There is an additional argument that 
determines whether this is a weighted average or not. The sign 
of the integration may be passed as well, i.e. integrating only 
positive or negative area, by passing a 1, 0, or -1 to integ_sign. 

Parameters:
    layer: HeightLayer 
    data: 1D NumPy array of data values 
    height: 1D NumPy array of height values 
    integ_sign: The sign of the area to integrate (-1: negative; 1: positive; 0: both; default: 0)
    weighted: Boolean determining whether or not the integration is weighted by the coordinate array 

    )pbdoc");

    m_layer.def(
        "integrate_layer_trapz",
        [](sharp::PressureLayer layer, const_prof_arr_t data,
           const_prof_arr_t pressure, const int integ_sign,
           const bool weighted) {
            return sharp::integrate_layer_trapz(layer, data.data(),
                                                pressure.data(), data.size(),
                                                integ_sign, weighted);
        },
        nb::arg("layer"), nb::arg("data"), nb::arg("pressure"),
        nb::arg("integ_sign") = 0, nb::arg("weighted") = false,
        R"pbdoc(
Returns the trapezoidal integration of the given data array over 
the given PressureLayer. There is an additional argument that 
determines whether this is a weighted average or not. The sign 
of the integration may be passed as well, i.e. integrating only 
positive or negative area, by passing a 1, 0, or -1 to integ_sign. 

Parameters:
    layer: PressureLayer
    data: 1D NumPy array of data values
    pressure: 1D NumPy array of pressure values 
    integ_sign: The sign of the area to integrate (-1: negative; 1: positive; 0: both; default: 0)
    weighted: Boolean determining whether or not the integration is weighted by the coordinate array

    )pbdoc");
}

#endif
