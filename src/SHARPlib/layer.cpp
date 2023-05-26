
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
                           int N) noexcept {
    constexpr auto bottom_comp = std::greater<float>();
    constexpr auto top_comp = std::less<float>();

    return get_layer_index(layer, pressure, N, bottom_comp, top_comp);
}

LayerIndex get_layer_index(HeightLayer& layer, const float height[],
                           int N) noexcept {
    constexpr auto bottom_comp = std::less<float>();
    constexpr auto top_comp = std::greater<float>();

    return get_layer_index(layer, height, N, bottom_comp, top_comp);
}

PressureLayer height_layer_to_pressure(HeightLayer layer,
                                       const float pressure[],
                                       const float height[], int num_levs,
                                       bool isAGL) noexcept {
    if (isAGL) {
        layer.bottom += height[0];
        layer.top += height[0];
    }

    float pbot = interp_height(layer.bottom, height, pressure, num_levs);
    float ptop = interp_height(layer.top, height, pressure, num_levs);

    return {pbot, ptop};
}

HeightLayer pressure_layer_to_height(PressureLayer layer,
                                     const float pressure[],
                                     const float height[], int num_levs,
                                     bool toAGL) noexcept {
    float zbot = interp_pressure(layer.bottom, pressure, height, num_levs);
    float ztop = interp_pressure(layer.top, pressure, height, num_levs);

    if (toAGL) {
        zbot -= height[0];
        ztop -= height[0];
    }

    return {zbot, ztop};
}

float layer_mean(PressureLayer layer, const float pressure[],
                 const float data_arr[], int num_levs) noexcept {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    if (layer.bottom > pressure[0]) {
        layer.bottom = pressure[0];
    }

    if (layer.top < pressure[num_levs - 1]) {
        layer.top = pressure[num_levs - 1];
    }

    float mean =
        integrate_layer_trapz(layer, data_arr, pressure, num_levs, 0, true);
    return mean;
}

float layer_mean(HeightLayer layer, const float height[],
                 const float pressure[], const float data_arr[], int num_levs,
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

    if (layer.top > height[num_levs - 1]) {
        layer.top = height[num_levs - 1];
    }

    PressureLayer pres_layer =
        height_layer_to_pressure(layer, pressure, height, num_levs, false);

    float mean = integrate_layer_trapz(pres_layer, data_arr, pressure, num_levs,
                                       0, true);
    return mean;
}

}  // end namespace sharp

namespace sharp::exper {}  // end namespace sharp::exper
