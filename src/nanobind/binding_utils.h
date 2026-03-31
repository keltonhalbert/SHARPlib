#ifndef SHARPLIB_BINDING_UTILS_H
#define SHARPLIB_BINDING_UTILS_H

#include <memory>

#include "sharplib_types.h"

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

#endif  // SHARPLIB_BINDING_UTILS_H
