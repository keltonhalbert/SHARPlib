// clang-format off
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

// cland-format on
#include <SHARPlib/layer.h>
#include <cstddef>

namespace py = pybind11;
// clang-format on
PYBIND11_MODULE(layer, m) {
    m.doc() = "Utility functions for operating on layers of the atmosphere.";

    // Bind the constructors, named fields, and default arguments
    py::class_<sharp::HeightLayer>(m, "HeightLayer")
        .def(py::init<float, float, float>(), py::arg("bottom"), py::arg("top"),
             py::arg("delta") = 100.0f)
        .def_readwrite("bottom", &sharp::HeightLayer::bottom,
                       "The bottom of the HeightLayer (meters)")
        .def_readwrite("top", &sharp::HeightLayer::top,
                       "The top of the HeightLayer (meters)")
        .def_readwrite(
            "delta", &sharp::HeightLayer::delta,
            "The HeightLayer delta (increment) to use if iterating (meters).");

    py::class_<sharp::PressureLayer>(m, "PressureLayer")
        .def(py::init<float, float, float>(), py::arg("bottom"), py::arg("top"),
             py::arg("delta") = -1000.0f)
        .def_readwrite("bottom", &sharp::PressureLayer::bottom,
                       "The bottom of the PressureLayer (Pa)")
        .def_readwrite("top", &sharp::PressureLayer::top,
                       "The top of the PressureLayer (Pa)")
        .def_readwrite(
            "delta", &sharp::PressureLayer::delta,
            "The PressureLayer delta (increment) to use if iterating (Pa)");

    py::class_<sharp::LayerIndex>(m, "LayerIndex")
        .def(py::init<std::ptrdiff_t, std::ptrdiff_t>(), py::arg("kbot"),
             py::arg("ktop"))
        .def_readwrite("kbot", &sharp::LayerIndex::kbot,
                       "The bottom index of a layer on a coordinate array "
                       "(pressure or height).")
        .def_readwrite("ktop", &sharp::LayerIndex::ktop,
                       "The top index of a layer on a coordinate array "
                       "(pressure or height).");

    // Set up the numpy dtypes so we can vectorize functions that
    // return a sharp::HeightLayer or sharp::PressureLayer
    PYBIND11_NUMPY_DTYPE(sharp::HeightLayer, bottom, top, delta);
    PYBIND11_NUMPY_DTYPE(sharp::PressureLayer, bottom, top, delta);
    PYBIND11_NUMPY_DTYPE(sharp::LayerIndex, kbot, ktop);

    m.def(
        "get_layer_index",
        [](sharp::PressureLayer& layer,
           py::array_t<float, py::array::c_style | py::array::forcecast>
               pressure_array) -> sharp::LayerIndex {
            // Request buffer info from the py::array
            py::buffer_info buf = pressure_array.request();

            // Check the shape of the array (1D expected)
            if (buf.ndim != 1) {
                throw std::runtime_error("pressure must be a 1D array");
            }

            // Call the original function
            return sharp::get_layer_index(layer, static_cast<float*>(buf.ptr),
                                          buf.size);
        },
        py::arg("layer"), py::arg("pressure"),
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
           py::array_t<float, py::array::c_style | py::array::forcecast>
               height_array) -> sharp::LayerIndex {
            // Request buffer info from the py::array
            py::buffer_info buf = height_array.request();

            // Check the shape of the array (1D expected)
            if (buf.ndim != 1) {
                throw std::runtime_error("height must be a 1D array");
            }

            // Call the original function
            return sharp::get_layer_index(layer, static_cast<float*>(buf.ptr),
                                          buf.size);
        },
        py::arg("layer"), py::arg("height"),
        R"pbdoc(
              Finds the array indices corresponding to the given HeightLayer.

              Parameters:
                  layer: HeightLayer
                  height: 1D NumPy array of heights

              Returns:
                  A LayerIndex with {kbot, ktop}.
          )pbdoc");
}
