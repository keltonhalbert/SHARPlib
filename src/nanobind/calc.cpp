// clang-format off
#include <nanobind/nanobind.h>

// cland-format on
#include "constants_bindings.h"
#include "interp_bindings.h"
#include "layer_bindings.h"
#include "thermo_bindings.h"
#include "winds_bindings.h"
#include "parcel_bindings.h"
#include "params_bindings.h"

// clang-format on
NB_MODULE(_calc, m) {
    m.doc() =
        "Sounding and Hodograph Analysis and Research Program Library "
        "(SHARPlib)";

    make_constants_bindings(m);
    make_interp_bindings(m);
    make_layer_bindings(m);
    make_thermo_bindings(m);
    make_winds_bindings(m);
    make_parcel_bindings(m);
    make_params_bindings(m);
}
