#ifndef __SHARPLIB_LAYER_BINDINGS
#define __SHARPLIB_LAYER_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>

#include <SHARPlib/layer.h>
#include "sharplib_types.h"

namespace nb = nanobind;

void make_layer_bindings(nb::module_ m) {
    nb::module_ m_layer = m.def_submodule("layer", 
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
            "The HeightLayer delta (increment) to use if iterating (meters).");

    nb::class_<sharp::PressureLayer>(m_layer, "PressureLayer")
        .def(nb::init<float, float, float>(), nb::arg("bottom"), nb::arg("top"),
            nb::arg("delta") = -1000.0f)
        .def_rw("bottom", &sharp::PressureLayer::bottom,
                "The bottom of the PressureLayer (Pa)")
        .def_rw("top", &sharp::PressureLayer::top,
                "The top of the PressureLayer (Pa)")
        .def_rw("delta", &sharp::PressureLayer::delta,
                "The PressureLayer delta (increment) to use if iterating (Pa)");

    nb::class_<sharp::LayerIndex>(m_layer, "LayerIndex")
        .def(nb::init<std::ptrdiff_t, std::ptrdiff_t>(), nb::arg("kbot"),
            nb::arg("ktop"))
        .def_rw("kbot", &sharp::LayerIndex::kbot,
                "The bottom index of a layer on a coordinate array "
                "(pressure or height).")
        .def_rw("ktop", &sharp::LayerIndex::ktop,
                "The top index of a layer on a coordinate array "
                "(pressure or height).");

    m_layer.def(
        "get_layer_index",
        [](sharp::PressureLayer& layer,
        const_prof_arr_t pressure_array) -> sharp::LayerIndex {
            return sharp::get_layer_index(layer, pressure_array.data(),
                                        pressure_array.size());
        },
        nb::arg("layer"), nb::arg("pressure"),
        R"pbdoc(
Finds the array indices corresponding to the given PressureLayer.

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
Finds the array indices corresponding to the given HeightLayer.

Parameters:
    layer: HeightLayer
    pressure: 1D NumPy array of heights

Returns:
    A LayerIndex with {kbot, ktop}.
        )pbdoc");
}

#endif