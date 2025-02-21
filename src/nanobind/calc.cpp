// clang-format off
#include <nanobind/nanobind.h>

// cland-format on
#include "thermo_bindings.h"
#include "interp_bindings.h"
#include "layer_bindings.h"
#include "parcel_bindings.h"

// clang-format on
NB_MODULE(calc, m) {
    m.doc() = "Sounding and Hodograph Analysis and Research Program Library "
              "(SHARPlib)";

    make_thermo_bindings(m);
    make_interp_bindings(m);
    make_layer_bindings(m);
}