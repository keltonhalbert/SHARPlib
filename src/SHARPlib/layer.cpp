
/**
 * \file
 * \brief Utilities and structures useful for grouping objects, conducting QC, and things that otherwise don't have a home.  
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
#include <functional>
#include <stdexcept>
#include <cmath>

#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#define FMT_HEADER_ONLY
#include <fmt/core.h>

namespace sharp {


HeightLayer::HeightLayer(float bottom, float top, float delta) {
    if (bottom > top) {
        throw std::range_error(
            fmt::format("RangeError: The top of the height layer must be > the bottom of the height layer. Got hbot: {0} and htop: {1}", bottom, top)
        );
    }
    this->bottom = bottom;
    this->top = top;
	this->delta = delta;
	this->coord = LayerCoordinate::height;
}


PressureLayer::PressureLayer(float bottom, float top, float delta) {
    if (bottom < top) {
        throw std::range_error(
            fmt::format("RangeError: The bottom of the pressure layer must be > the top of the pressure layer. Got pbot: {0} and ptop: {1}", bottom, top)
        );
    }
    this->bottom = bottom;
    this->top = top;
    this->delta = delta;
	this->coord = LayerCoordinate::pressure;
}

LayerIndex get_layer_index(PressureLayer& layer, const float* pressure, int N) noexcept {

	auto bottom_comp = std::greater_equal<float>();
	auto top_comp = std::less_equal<float>();

	LayerIndex idx = get_layer_index(layer, pressure, N, bottom_comp, top_comp);
	return {idx.kbot, idx.ktop};
}

LayerIndex get_layer_index(HeightLayer& layer, const float* height, int N) noexcept {

	auto bottom_comp = std::less_equal<float>();
	auto top_comp = std::greater_equal<float>();

	LayerIndex idx = get_layer_index(layer, height, N, bottom_comp, top_comp);
	return {idx.kbot, idx.ktop};
}


PressureLayer height_layer_to_pressure(HeightLayer layer, 
                const float* pressure, const float* height, 
                int num_levs, bool isAGL) noexcept {

    if (isAGL) {
        layer.bottom += height[0];
        layer.top += height[0];
    }

    float pbot = interp_height(layer.bottom, height, pressure, num_levs);
    float ptop = interp_height(layer.top, height, pressure, num_levs);
    
    return {pbot, ptop};
}

HeightLayer pressure_layer_to_height(PressureLayer layer, 
                const float* pressure, const float* height, 
                int num_levs, bool toAGL) noexcept {
    float zbot = interp_pressure(layer.bottom, pressure, height, num_levs);
    float ztop = interp_pressure(layer.top, pressure, height, num_levs);

    if (toAGL) {
        zbot -= height[0];
        ztop -= height[0];
    }

    return {zbot, ztop};
}

float layer_min(PressureLayer layer, const float* pressure,
                const float* data_arr, int N, 
                float* pres_of_min) noexcept {

    auto comp = std::less<float>();
    return layer_minmax(layer, pressure, data_arr, N, pres_of_min, comp);
}

float layer_min(HeightLayer layer, const float* height,
                const float* data_arr, int N, 
                float* hght_of_min) noexcept {

    auto comp = std::less<float>();
    return layer_minmax(layer, height, data_arr, N, hght_of_min, comp);
}

float layer_max(PressureLayer layer, const float* pressure,
                const float* data_arr, int N, 
                float* pres_of_max) noexcept {

    auto comp = std::greater<float>();
    return layer_minmax(layer, pressure, data_arr, N, pres_of_max, comp);
}

float layer_max(HeightLayer layer, const float* height,
                const float* data_arr, int N, 
                float* hght_of_max) noexcept {

    auto comp = std::greater<float>();
    return layer_minmax(layer, height, data_arr, N, hght_of_max, comp);
}

float layer_mean(PressureLayer layer,   const float* pressure,
                 const float* data_arr, int num_levs) noexcept {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search.
    LayerIndex layer_idx = get_layer_index(layer, pressure, num_levs);

    // start with interpolated bottom layer
    float val_bot = interp_pressure(layer.bottom, pressure, data_arr, num_levs);
    float val_top = MISSING;
    float pbot = layer.bottom;
    float ptop = 0.0;
    float avg_val = 0.0;
    float weight = 0.0;

    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((pressure[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        val_top = data_arr[k];
        ptop = pressure[k];
        avg_val += ((val_top + val_bot) / 2.0) * (pbot - ptop);
        weight += (pbot - ptop);
        
        val_bot = val_top;
        pbot = ptop;
    }
    
    val_top = interp_pressure(layer.top, pressure, data_arr, num_levs);
    ptop = layer.top;

    avg_val += ((val_top + val_bot) / 2.0) * (pbot - ptop);
    weight += (pbot - ptop);

    return avg_val / weight;
}


float layer_mean(HeightLayer layer_agl, const float* height, 
                 const float* pressure, const float* data_arr, 
                 int num_levs) noexcept {
#ifndef NO_QC
    if ((layer_agl.bottom == MISSING) || (layer_agl.top == MISSING)) {
        return MISSING;
    }
#endif

    PressureLayer pres_layer = height_layer_to_pressure(
                                layer_agl, pressure, height,
                                num_levs, true
                            );

    return layer_mean(pres_layer, pressure, data_arr, num_levs);
}


} // end namespace sharp


namespace sharp::exper {

} // end namespace sharp::exper



