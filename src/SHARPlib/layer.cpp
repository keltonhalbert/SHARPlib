
/**
 * \file
 * \brief Utilities and structures useful for grouping objects,<!--
 * --> conducting QC, and things that otherwise don't have a home.
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-11-02
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#include <cmath>
#include <functional>
#include <stdexcept>

#define FMT_HEADER_ONLY
#include <fmt/core.h>

namespace sharp {

HeightLayer::HeightLayer(float bottom, float top, float delta) {
    if (bottom > top) {
        throw std::range_error(fmt::format(
            "RangeError: The top of the height layer must be > the bottom of "
            "the height layer. Got hbot: {0} and htop: {1}",
            bottom, top));
    }
    this->bottom = bottom;
    this->top = top;
    this->delta = delta;
}

PressureLayer::PressureLayer(float bottom, float top, float delta) {
    if (bottom < top) {
        throw std::range_error(fmt::format(
            "RangeError: The bottom of the pressure layer must be > the top of "
            "the pressure layer. Got pbot: {0} and ptop: {1}",
            bottom, top));
    }
    this->bottom = bottom;
    this->top = top;
    this->delta = delta;
}

LayerIndex get_layer_index(PressureLayer& layer, const float pressure[],
                           const int N) noexcept {
    static constexpr auto bottom_comp = std::greater<float>();
    static constexpr auto top_comp = std::less<float>();

    return get_layer_index(layer, pressure, N, bottom_comp, top_comp);
}

LayerIndex get_layer_index(HeightLayer& layer, const float height[],
                           const int N) noexcept {
    static constexpr auto bottom_comp = std::less<float>();
    static constexpr auto top_comp = std::greater<float>();

    return get_layer_index(layer, height, N, bottom_comp, top_comp);
}

PressureLayer height_layer_to_pressure(HeightLayer layer,
                                       const float pressure[],
                                       const float height[], const int N,
                                       const bool isAGL) noexcept {
    if (isAGL) {
        layer.bottom += height[0];
        layer.top += height[0];
    }

    if ((layer.bottom < height[0]) || (layer.top > height[N-1])) {
        return {MISSING, MISSING};
    }

    const float pbot = interp_height(layer.bottom, height, pressure, N);
    const float ptop = interp_height(layer.top, height, pressure, N);

    return {pbot, ptop};
}

HeightLayer pressure_layer_to_height(PressureLayer layer,
                                     const float pressure[],
                                     const float height[], const int N,
                                     const bool toAGL) noexcept {
    if ((layer.bottom > pressure[0]) || (layer.top < pressure[N-1])) {
        return {MISSING, MISSING};
    }

    float zbot = interp_pressure(layer.bottom, pressure, height, N);
    float ztop = interp_pressure(layer.top, pressure, height, N);

    if (toAGL) {
        zbot -= height[0];
        ztop -= height[0];
    }

    return {zbot, ztop};
}

float layer_mean(PressureLayer layer, const float pressure[],
                 const float data_arr[], const int N) noexcept {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    if (layer.bottom > pressure[0]) {
        layer.bottom = pressure[0];
    }

    if (layer.top < pressure[N - 1]) {
        layer.top = pressure[N - 1];
    }

    return integrate_layer_trapz(layer, data_arr, pressure, N, 0, true);
}

float layer_mean(HeightLayer layer, const float height[],
                 const float pressure[], const float data_arr[], const int N,
                 const bool isAGL) noexcept {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    if (isAGL) {
        layer.bottom += height[0];
        layer.top += height[0];
    }

    if (layer.bottom < height[0]) {
        layer.bottom = height[0];
    }

    if (layer.top > height[N - 1]) {
        layer.top = height[N - 1];
    }

    PressureLayer pres_layer =
        height_layer_to_pressure(layer, pressure, height, N, false);

    return integrate_layer_trapz(pres_layer, data_arr, pressure, N, 0, true);
}

}  // end namespace sharp

