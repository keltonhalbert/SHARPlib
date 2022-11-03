
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
#include <stdexcept>
#include <cmath>

#include "constants.h"
#include "interp.h"
#include "utils.h"

namespace sharp {


HeightLayer::HeightLayer(float bot, float top, float delta) {
    if (bot > top) {
        throw std::range_error("The top of the height layer must be > the bottom of the height layer.");
    }
    zbot = bot;
    ztop = top;
    dz = delta;
}


HeightLayer::HeightLayer(float bot, float top) {
    if (bot > top) {
        throw std::range_error("The top of the height layer must be > the bottom of the height layer.");
    }
    zbot = bot;
    ztop = top;
    dz = 100.0;
}


PressureLayer::PressureLayer(float bot, float top, float delta) {
    if (bot < top) {
        throw std::range_error("The top of the pressure layer must be < the bottom of the pressure layer.");
    }
    pbot = bot;
    ptop = top;
    dp = delta;
}


PressureLayer::PressureLayer(float bot, float top) {
    if (bot < top) {
        throw std::range_error("The top of the pressure layer must be < the bottom of the pressure layer.");
    }
    pbot = bot;
    ptop = top;
    dp = 1;
}

LayerIndex get_layer_index(PressureLayer& layer, 
                           const float* pressure, int num_levs) {
    // bounds check our search, modifying
    // our PressureLayer if need be.
    if ((layer.pbot > pressure[0])) {
        layer.pbot = pressure[0];
    }
    if ((layer.ptop < pressure[num_levs-1])) {
        layer.ptop = pressure[num_levs-1];
    }

    // find the lowest observation, excluding
    // the exact level for interpolation reasons
    int lower_idx = 0;
    while (pressure[lower_idx] >= layer.pbot) {
        lower_idx++;
    }

    // find the highest observation, excluding
    // the exact level for interpolation reasons
    int upper_idx = num_levs-1;
    while (pressure[upper_idx] <= layer.ptop) {
        upper_idx--;
    }

    return {lower_idx, upper_idx};
}


LayerIndex get_layer_index(HeightLayer& layer, 
                           const float* height, int num_levs) {
    // bounds check our search, modifying
    // our PressureLayer if need be. 
    if ((layer.zbot < height[0])) {
        layer.zbot = height[0];
    }
    if ((layer.ztop > height[num_levs-1])) {
        layer.ztop = height[num_levs-1];
    }

    // find the lowest observation, excluding
    // the exact level for interpolation reasons
    int lower_idx = 0;
    while (height[lower_idx] <= layer.zbot) {
        lower_idx++;
    }

    // find the highest observation, excluding
    // the exact level for interpolation reasons
    int upper_idx = num_levs-1;
    while (height[upper_idx] >= layer.ztop) {
        upper_idx--;
    }

    return {lower_idx, upper_idx};
}


float max_value(PressureLayer layer,   const float* pressure,
                const float* data_arr, int num_levs) {
    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. 
    LayerIndex layer_idx = get_layer_index(layer, pressure, num_levs);

    // start with the interpolated bottom of the layer
    float max_val = interp_pressure(layer.pbot, pressure, data_arr, num_levs);
    float pres_of_max = layer.pbot;

    int lower_idx = layer_idx.kbot;
    int upper_idx = layer_idx.ktop;
    float pval, dval;
    for (int k = lower_idx; k <= upper_idx; k++) {
#ifndef NO_QC
        if ((pressure[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        pval = pressure[k];
        dval = data_arr[k];

        if (dval > max_val) {
            max_val = dval;
            pres_of_max = pval;
        }
    }

    // now check the interpolated top of the layer
    pval = layer.ptop;
    dval = interp_pressure(layer.ptop, pressure, data_arr, num_levs);
    if (dval > max_val) {
        max_val = dval;
        pres_of_max = pval;
    }

    return max_val;
}


float max_value(HeightLayer layer,     const float* height,
                const float* data_arr, int num_levs) {
    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. 
    LayerIndex layer_idx = get_layer_index(layer, height, num_levs);

    // start with the interpolated bottom of the layer
    float max_val = interp_height(layer.zbot, height, data_arr, num_levs);
    float hght_of_max = layer.zbot;

    int lower_idx = layer_idx.kbot;
    int upper_idx = layer_idx.ktop;
    float zval, dval;
    for (int k = lower_idx; k <= upper_idx; k++) {
#ifndef NO_QC
        if ((height[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        zval = height[k];
        dval = data_arr[k];

        if (dval > max_val) {
            max_val = dval;
            hght_of_max = zval;
        }
    }

    // now check the interpolated top of the layer
    zval = layer.ztop;
    dval = interp_height(layer.ztop, height, data_arr, num_levs);
    if (dval > max_val) {
        max_val = dval;
        hght_of_max = zval;
    }

    return max_val;
}


float min_value(PressureLayer layer,   const float* pressure,
                const float* data_arr, int num_levs) {
    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. 
    LayerIndex layer_idx = get_layer_index(layer, pressure, num_levs);

    // start with the interpolated bottom of the layer
    float min_val = interp_pressure(layer.pbot, pressure, data_arr, num_levs);
    float pres_of_min = layer.pbot;

    int lower_idx = layer_idx.kbot;
    int upper_idx = layer_idx.ktop;
    float pval, dval;
    for (int k = lower_idx; k <= upper_idx; k++) {
#ifndef NO_QC
        if ((pressure[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        pval = pressure[k];
        dval = data_arr[k];

        if (dval < min_val) {
            min_val = dval;
            pres_of_min = pval;
        }
    }

    // now check the interpolated top of the layer
    pval = layer.ptop;
    dval = interp_pressure(layer.ptop, pressure, data_arr, num_levs);
    if (dval < min_val) {
        min_val = dval;
        pres_of_min = pval;
    }

    return min_val;
}


float min_value(HeightLayer layer,     const float* height,
                const float* data_arr, int num_levs) {
    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search.
    LayerIndex layer_idx = get_layer_index(layer, height, num_levs);

    // start with the interpolated bottom of the layer
    float min_val = interp_height(layer.zbot, height, data_arr, num_levs);
    float hght_of_min = layer.zbot;

    int lower_idx = layer_idx.kbot;
    int upper_idx = layer_idx.ktop;
    float zval, dval;
    for (int k = lower_idx; k <= upper_idx; k++) {
#ifndef NO_QC
        if ((height[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        zval = height[k];
        dval = data_arr[k];

        if (dval < min_val) {
            min_val = dval;
            hght_of_min = zval;
        }
    }

    // now check the interpolated top of the layer
    zval = layer.ztop;
    dval = interp_height(layer.ztop, height, data_arr, num_levs);
    if (dval < min_val) {
        min_val = dval;
        hght_of_min = zval;
    }

    return min_val;
}


} // end namespace sharp


namespace sharp::exper {

} // end namespace sharp::exper



