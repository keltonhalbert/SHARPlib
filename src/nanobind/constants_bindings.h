#ifndef SHARPLIB_CONSTANTS_BINDINGS
#define SHARPLIB_CONSTANTS_BINDINGS

// clang-format off
#include <nanobind/nanobind.h>

// clang-format on
#include <SHARPlib/constants.h>

namespace nb = nanobind;

void make_constants_bindings(nb::module_ m) {
    nb::module_ m_const =
        m.def_submodule("constants",
                        "Sounding and Hodograph Analysis and Research Program "
                        "Library (SHARPlib) :: Constants");

    m_const.attr("MISSING") = sharp::MISSING;
}

#endif
