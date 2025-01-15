#ifndef COMMON_SHARPLIB_TYPES
#define COMMON_SHARPLIB_TYPES

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

// read only input arrays
using const_prof_arr_t =
    nb::ndarray<float, nb::ndim<1>, nb::device::cpu, nb::c_contig, nb::ro>;

// read/write arrays
using prof_arr_t =
    nb::ndarray<float, nb::ndim<1>, nb::device::cpu, nb::c_contig>;

#endif  // COMMON_SHARPLIB_TYPES
