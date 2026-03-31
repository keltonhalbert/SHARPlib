#ifndef SHARPLIB_LAYER_BINDINGS
#define SHARPLIB_LAYER_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/string.h>

// clang-format on
#include <SHARPlib/layer.h>
#include <fmt/core.h>

#include "binding_utils.h"
#include "sharplib_types.h"

namespace nb = nanobind;

template <typename Layer>
void bind_get_layer_index(nb::module_& mod, const char* doc_template,
                          const char* layer_type, const char* coord_name) {
    std::string doc = fmt::format(doc_template, layer_type, coord_name,
                                  coord_name, layer_type);
    mod.def(
        "get_layer_index",
        [](Layer& layer, const_prof_arr_t coord) -> sharp::LayerIndex {
            return sharp::get_layer_index(layer, coord.data(), coord.size());
        },
        nb::arg("layer"), nb::arg(coord_name), doc.c_str());
}

template <typename Layer>
void bind_layer_min(nb::module_& mod, const char* doc_template,
                    const char* layer_type, const char* coord_name) {
    std::string doc = fmt::format(doc_template, layer_type, coord_name);
    mod.def(
        "layer_min",
        [](Layer& layer, const_prof_arr_t coord, const_prof_arr_t data) {
            check_equal_sizes(coord, data);
            float lvl_of_min;
            float min = sharp::layer_min(layer, coord.data(), data.data(),
                                         data.size(), &lvl_of_min);
            return std::make_tuple(min, lvl_of_min);
        },
        nb::arg("layer"), nb::arg(coord_name), nb::arg("data"), doc.c_str());
}

template <typename Layer>
void bind_layer_max(nb::module_& mod, const char* doc_template,
                    const char* layer_type, const char* coord_name) {
    std::string doc = fmt::format(doc_template, layer_type, coord_name);
    mod.def(
        "layer_max",
        [](Layer& layer, const_prof_arr_t coord, const_prof_arr_t data) {
            check_equal_sizes(coord, data);
            float lvl_of_max;
            float max = sharp::layer_max(layer, coord.data(), data.data(),
                                         data.size(), &lvl_of_max);
            return std::make_tuple(max, lvl_of_max);
        },
        nb::arg("layer"), nb::arg(coord_name), nb::arg("data"), doc.c_str());
}

template <typename Layer>
void bind_layer_mean(nb::module_& mod, const char* doc) {
    if constexpr (Layer::coord == sharp::LayerCoordinate::height) {
        mod.def(
            "layer_mean",
            [](Layer layer, const_prof_arr_t height, const_prof_arr_t pressure,
               const_prof_arr_t data, const bool isAGL) {
                check_equal_sizes(height, pressure, data);
                return sharp::layer_mean(layer, height.data(), pressure.data(),
                                         data.data(), data.size(), isAGL);
            },
            nb::arg("layer"), nb::arg("height"), nb::arg("pressure"),
            nb::arg("data"), nb::arg("isAGL") = false, doc);
    } else {
        mod.def(
            "layer_mean",
            [](Layer layer, const_prof_arr_t pressure, const_prof_arr_t data) {
                check_equal_sizes(pressure, data);
                return sharp::layer_mean(layer, pressure.data(), data.data(),
                                         data.size());
            },
            nb::arg("layer"), nb::arg("pressure"), nb::arg("data"), doc);
    }
}

template <typename Layer>
void bind_integrate_layer_trapz(nb::module_& mod, const char* doc_template,
                                const char* layer_type,
                                const char* coord_name) {
    std::string doc = fmt::format(doc_template, layer_type, coord_name);
    mod.def(
        "integrate_layer_trapz",
        [](Layer layer, const_prof_arr_t data, const_prof_arr_t coord,
           const int integ_sign, const bool weighted) {
            check_equal_sizes(coord, data);
            return sharp::integrate_layer_trapz(layer, data.data(),
                                                coord.data(), data.size(),
                                                integ_sign, weighted);
        },
        nb::arg("layer"), nb::arg("data"), nb::arg(coord_name),
        nb::arg("integ_sign") = 0, nb::arg("weighted") = false, doc.c_str());
}

inline void make_layer_bindings(nb::module_ m) {
    nb::module_ m_layer = m.def_submodule(
        "layer",
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib) :: Atmospheric Layer Routines");

    nb::class_<sharp::HeightLayer>(m_layer, "HeightLayer")
        .def(nb::init<>())
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
        .def(nb::init<>())
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

    const char* get_layer_index_template_doc =
        R"pbdoc(
Finds the interior array indices encapsulated by the given {0}. 
More specifically, the returned LayerIndex excludes the exact top and bottom 
values (in mathematical notation, [bottom, top]). This behavior is due to the 
fact many algorithms use interpolation to get the exact values of arbitrary
top/bottom locations within a profile. 

If layer.bottom or layer.top are out of bounds, this function 
will truncate the layer to the coordinate range of data provided 
by coord[] in an attempt to gracefully continue and produce a result. 
This will modify the value of the given {3}.

Parameters
----------
layer : nwsspc.sharp.calc.layer.{0}
{1} : numpy.ndarray[dtype=float32]
    1D NumPy array of {2} values

Returns
-------
nwsspc.sharp.calc.layer.LayerIndex
    A LayerIndex with {{kbot, ktop}}.
)pbdoc";

    bind_get_layer_index<sharp::PressureLayer>(
        m_layer, get_layer_index_template_doc, "PressureLayer", "pressure");
    bind_get_layer_index<sharp::HeightLayer>(
        m_layer, get_layer_index_template_doc, "HeightLayer", "height");

    m_layer.def(
        "height_layer_to_pressure",
        [](sharp::HeightLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t height, const bool isAGL) -> sharp::PressureLayer {
            check_equal_sizes(pressure, height);
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

Parameters
----------
layer : nwsspc.sharp.calc.layer.HeightLayer
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure
height : numpy.ndarray[dtype=float32]
    1D NumPy array of heights
isAGL : bool 
    Whether or not the station height needs to be added for interpolation (default: False) 

Returns
-------
nwsspc.sharp.calc.layer.PressureLayer
    The HeightLayer converted to a PressureLayer
)pbdoc");

    m_layer.def(
        "pressure_layer_to_height",
        [](sharp::PressureLayer& layer, const_prof_arr_t pressure,
           const_prof_arr_t height, const bool toAGL) -> sharp::HeightLayer {
            check_equal_sizes(pressure, height);
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

Parameters
----------
layer : nwsspc.sharp.calc.layer.PressureLayer
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure 
height : numpy.ndarray[dtype=float32]
    1D NumPy array of heights
toAGL : bool
    Whether or not to subtract the station height from the HeightLayer (default: False)

Returns 
-------
nwsspc.sharp.calc.layer.HeightLayer
    The PressureLayer converted to a HeightLayer
)pbdoc");

    const char* layer_min_template_doc =
        R"pbdoc(
Returns the minimum value of the data array within the given {0}. The 
function bounds checks the layer by calling get_layer_index. 

Parameters
----------
layer : nwsspc.sharp.calc.layer.{0} 
{1} : numpy.ndarray[dtype=float32]
    1D NumPy array of {1} values 
data : numpy.ndarray[dtype=float32]
    1D NumPy array of data values

Returns
-------
tuple[float, float]
    (min_value, level_of_min)
)pbdoc";

    bind_layer_min<sharp::HeightLayer>(m_layer, layer_min_template_doc,
                                       "HeightLayer", "height");
    bind_layer_min<sharp::PressureLayer>(m_layer, layer_min_template_doc,
                                         "PressureLayer", "pressure");

    const char* layer_max_template_doc =
        R"pbdoc(
Returns the maximum value of the data array within the given {0}. The 
function bounds checks the layer by calling get_layer_index. 

Parameters
----------
layer : nwsspc.sharp.calc.layer.{0} 
{1} : numpy.ndarray[dtype=float32] 
    1D NumPy array of {1} values 
data : numpy.ndarray[dtype=float32] 
    1D NumPy array of data values

Returns
-------
tuple[float, float]
    (max_value, level_of_max)
)pbdoc";

    bind_layer_max<sharp::HeightLayer>(m_layer, layer_max_template_doc,
                                       "HeightLayer", "height");
    bind_layer_max<sharp::PressureLayer>(m_layer, layer_max_template_doc,
                                         "PressureLayer", "pressure");

    const char* layer_mean_pres_doc =
        R"pbdoc(
Computes the pressure-weighted mean value of a field over 
a given PressureLayer. 

Parameters
----------
layer : nwsspc.sharp.calc.layer.PressureLayer 
pressure : numpy.ndarray[dtype=float32] 
    1D NumPy array of pressure values 
data : numpy.ndarray[dtype=float32] 
    1D NumPy array of data values

Returns
-------
float
    The layer mean
)pbdoc";

    const char* layer_mean_hght_doc =
        R"pbdoc(
Computes the pressure-weighted mean value of a field over 
a given HeightLayer. 

Parameters
----------
layer : nwsspc.sharp.calc.layer.HeightLayer 
height : numpy.ndarray[dtype=float32] 
    1D NumPy array of height values
pressure : numpy.ndarray[dtype=float32]
    1D NumPy array of pressure values 
data : numpy.ndarray[dtype=float32]
    1D NumPy array of data values
isAGL : bool
    Whether or not the surface station height should be added to the HeightLayer (default: False)

Returns
-------
float
    The layer mean
)pbdoc";

    bind_layer_mean<sharp::PressureLayer>(m_layer, layer_mean_pres_doc);
    bind_layer_mean<sharp::HeightLayer>(m_layer, layer_mean_hght_doc);

    const char* integrate_trapz_template_doc =
        R"pbdoc(
Returns the trapezoidal integration of the given data array over 
the given {0}. There is an additional argument that 
determines whether this is a weighted average or not. The sign 
of the integration may be passed as well, i.e. integrating only 
positive or negative area, by passing a 1, 0, or -1 to integ_sign. 

Parameters
----------
layer : nwsspc.sharp.calc.layer.{0}
data : numpy.ndarray[dtype=float32] 
    1D NumPy array of data values
{1} : numpy.ndarray[dtype=float32] 
    1D NumPy array of {1} values 
integ_sign : int 
    The sign of the area to integrate (-1: negative; 1: positive; 0: both; default: 0)
weighted : bool 
    Boolean determining whether or not the integration is weighted by the coordinate array

Returns 
-------
float
    The integrated value
)pbdoc";

    bind_integrate_layer_trapz<sharp::HeightLayer>(
        m_layer, integrate_trapz_template_doc, "HeightLayer", "height");
    bind_integrate_layer_trapz<sharp::PressureLayer>(
        m_layer, integrate_trapz_template_doc, "PressureLayer", "pressure");
}

#endif
