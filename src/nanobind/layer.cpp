// clang-format off
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

// cland-format on
#include <SHARPlib/layer.h>
#include <cstddef>

namespace nb = nanobind;

using profile_arr_t = nb::ndarray<float, nb::ndim<1>, nb::device::cpu, nb::c_contig, nb::ro>;
// clang-format on
NB_MODULE(layer, m) {
    m.doc() = "Utility functions for operating on layers of the atmosphere.";

    // Bind the constructors, named fields, and default arguments
    nb::class_<sharp::HeightLayer>(m, "HeightLayer")
        .def(nb::init<float, float, float>(), nb::arg("bottom"), nb::arg("top"),
             nb::arg("delta") = 100.0f)
        .def_rw("bottom", &sharp::HeightLayer::bottom,
                "The bottom of the HeightLayer (meters)")
        .def_rw("top", &sharp::HeightLayer::top,
                "The top of the HeightLayer (meters)")
        .def_rw(
            "delta", &sharp::HeightLayer::delta,
            "The HeightLayer delta (increment) to use if iterating (meters).");

    nb::class_<sharp::PressureLayer>(m, "PressureLayer")
        .def(nb::init<float, float, float>(), nb::arg("bottom"), nb::arg("top"),
             nb::arg("delta") = -1000.0f)
        .def_rw("bottom", &sharp::PressureLayer::bottom,
                "The bottom of the PressureLayer (Pa)")
        .def_rw("top", &sharp::PressureLayer::top,
                "The top of the PressureLayer (Pa)")
        .def_rw("delta", &sharp::PressureLayer::delta,
                "The PressureLayer delta (increment) to use if iterating (Pa)");

    nb::class_<sharp::LayerIndex>(m, "LayerIndex")
        .def(nb::init<std::ptrdiff_t, std::ptrdiff_t>(), nb::arg("kbot"),
             nb::arg("ktop"))
        .def_rw("kbot", &sharp::LayerIndex::kbot,
                "The bottom index of a layer on a coordinate array "
                "(pressure or height).")
        .def_rw("ktop", &sharp::LayerIndex::ktop,
                "The top index of a layer on a coordinate array "
                "(pressure or height).");

    m.def(
        "get_layer_index",
        [](sharp::PressureLayer& layer,
           profile_arr_t pressure_array) -> sharp::LayerIndex {
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

    m.def(
        "get_layer_index",
        [](sharp::HeightLayer& layer,
           profile_arr_t height_array) -> sharp::LayerIndex {
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
