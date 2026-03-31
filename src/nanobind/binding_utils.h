#ifndef SHARPLIB_BINDING_UTILS_H
#define SHARPLIB_BINDING_UTILS_H

#include <SHARPlib/parcel.h>

#include <memory>

#include "sharplib_types.h"

namespace nb = nanobind;

template <typename T>
struct type_tag {
    using type = T;
};

template <typename First, typename... Rest>
void check_equal_sizes(const First& first, const Rest&... rest) {
    if (((first.size() != rest.size()) || ...)) {
        throw nb::buffer_error("All input arrays must have the same size!");
    }
}

template <typename Func>
out_arr_t make_output_array(std::size_t NZ, Func&& compute) {
    auto buf = std::make_unique<float[]>(NZ);
    compute(buf.get());
    float* raw = buf.release();
    nb::capsule owner(raw, [](void* p) noexcept { delete[] (float*)p; });
    return out_arr_t(raw, {NZ}, owner);
}

template <typename F>
void for_each_base_lifter(F&& f) {
    f(type_tag<sharp::lifter_wobus>{}, "wobus", "Wobus",
      "nwsspc.sharp.calc.parcel.lifter_wobus",
      "nwsspc.sharp.calc.parcel.lut_data_wobus",
      "nwsspc.sharp.calc.parcel.lifter_lut_wobus");
    f(type_tag<sharp::lifter_cm1>{}, "cm1", "CM1",
      "nwsspc.sharp.calc.parcel.lifter_cm1",
      "nwsspc.sharp.calc.parcel.lut_data_cm1",
      "nwsspc.sharp.calc.parcel.lifter_lut_cm1");
}

template <typename F>
void for_each_lifter(F&& f) {
    for_each_base_lifter([&](auto tag, const char* suffix, const char*,
                             const char* lifter_fqn, const char*,
                             const char* lifter_lut_fqn) {
        using Lft = typename decltype(tag)::type;
        f(type_tag<Lft>{}, lifter_fqn);
        f(type_tag<sharp::lifter_lut<Lft>>{}, lifter_lut_fqn);
    });
}

#endif  // SHARPLIB_BINDING_UTILS_H
