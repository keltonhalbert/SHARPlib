
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

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/utils.h>

namespace sharp {


HeightLayer::HeightLayer(float bot, float top, float delta) {
    if (bot > top) {
        throw std::range_error("RangeError: The top of the height layer must be > the bottom of the height layer.");
    }
    zbot = bot;
    ztop = top;
    dz = delta;
}


HeightLayer::HeightLayer(float bot, float top) {
    if (bot > top) {
        throw std::range_error("RangeError: The top of the height layer must be > the bottom of the height layer.");
    }
    zbot = bot;
    ztop = top;
    dz = 100.0;
}


PressureLayer::PressureLayer(float bot, float top, float delta) {
    if (bot < top) {
        throw std::range_error("RangeError: The top of the pressure layer must be < the bottom of the pressure layer.");
    }
    pbot = bot;
    ptop = top;
    dp = delta;
}


PressureLayer::PressureLayer(float bot, float top) {
    if (bot < top) {
        throw std::range_error("RangeError: The top of the pressure layer must be < the bottom of the pressure layer.");
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
                const float* data_arr, int num_levs,
                float* pres_of_max) {
#ifndef NO_QC
    if ((layer.pbot == MISSING) || (layer.ptop == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. 
    LayerIndex layer_idx = get_layer_index(layer, pressure, num_levs);

    // start with the interpolated bottom of the layer
    float max_val = interp_pressure(layer.pbot, pressure, data_arr, num_levs);
    float pr_max = layer.pbot;

    float pval, dval;
    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((pressure[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        pval = pressure[k];
        dval = data_arr[k];

        if (dval > max_val) {
            max_val = dval;
            pr_max = pval;
        }
    }

    // now check the interpolated top of the layer
    pval = layer.ptop;
    dval = interp_pressure(layer.ptop, pressure, data_arr, num_levs);
    if (dval > max_val) {
        max_val = dval;
        pr_max = pval;
    }

    // if it isn't nullptr
    if (pres_of_max) *pres_of_max = pr_max;
    return max_val;
}


float max_value(HeightLayer layer,     const float* height,
                const float* data_arr, int num_levs,
                float* hght_of_max) {
#ifndef NO_QC
    if ((layer.zbot == MISSING) || (layer.ztop == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. 
    LayerIndex layer_idx = get_layer_index(layer, height, num_levs);

    // start with the interpolated bottom of the layer
    float max_val = interp_height(layer.zbot, height, data_arr, num_levs);
    float ht_max = layer.zbot;

    float zval, dval;
    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((height[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        zval = height[k];
        dval = data_arr[k];

        if (dval > max_val) {
            max_val = dval;
            ht_max = zval;
        }
    }

    // now check the interpolated top of the layer
    zval = layer.ztop;
    dval = interp_height(layer.ztop, height, data_arr, num_levs);
    if (dval > max_val) {
        max_val = dval;
        ht_max = zval;
    }

    // if not nullptr
    if (hght_of_max) *hght_of_max = ht_max;
    return max_val;
}


float min_value(PressureLayer layer,   const float* pressure,
                const float* data_arr, int num_levs,
                float* pres_of_min) {
#ifndef NO_QC
    if ((layer.pbot == MISSING) || (layer.ptop == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. 
    LayerIndex layer_idx = get_layer_index(layer, pressure, num_levs);

    // start with the interpolated bottom of the layer
    float min_val = interp_pressure(layer.pbot, pressure, data_arr, num_levs);
    float pr_min = layer.pbot;

    float pval, dval;
    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((pressure[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        pval = pressure[k];
        dval = data_arr[k];

        if (dval < min_val) {
            min_val = dval;
            pr_min = pval;
        }
    }

    // now check the interpolated top of the layer
    pval = layer.ptop;
    dval = interp_pressure(layer.ptop, pressure, data_arr, num_levs);
    if (dval < min_val) {
        min_val = dval;
        pr_min = pval;
    }

    // if not nullptr
    if (pres_of_min) *pres_of_min = pr_min;
    return min_val;
}


float min_value(HeightLayer layer,     const float* height,
                const float* data_arr, int num_levs,
                float* hght_of_min) {
#ifndef NO_QC
    if ((layer.zbot == MISSING) || (layer.ztop == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search.
    LayerIndex layer_idx = get_layer_index(layer, height, num_levs);

    // start with the interpolated bottom of the layer
    float min_val = interp_height(layer.zbot, height, data_arr, num_levs);
    float ht_min = layer.zbot;

    float zval, dval;
    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((height[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        zval = height[k];
        dval = data_arr[k];

        if (dval < min_val) {
            min_val = dval;
            ht_min = zval;
        }
    }

    // now check the interpolated top of the layer
    zval = layer.ztop;
    dval = interp_height(layer.ztop, height, data_arr, num_levs);
    if (dval < min_val) {
        min_val = dval;
        ht_min = zval;
    }

    // if not nullptr
    if (hght_of_min) *hght_of_min = ht_min;
    return min_val;
}


float mean_value(PressureLayer layer,   const float* pressure,
                 const float* data_arr, int num_levs) {
#ifndef NO_QC
    if ((layer.pbot == MISSING) || (layer.ptop == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search.
    LayerIndex layer_idx = get_layer_index(layer, pressure, num_levs);

    // start with interpolated bottom layer
    float data_sum = interp_pressure(layer.pbot, pressure, data_arr, num_levs);
    float count = 1.0;

    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((pressure[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        data_sum += data_arr[k]; 
        count += 1.0;
    }

    // finish with interpolated top layer
    data_sum += interp_pressure(layer.ptop, pressure, data_arr, num_levs);
    count += 1.0;

    return data_sum / count;
}


float mean_value(HeightLayer layer,     const float* height, 
                 const float* data_arr, int num_levs) {
#ifndef NO_QC
    if ((layer.zbot == MISSING) || (layer.ztop == MISSING)) {
        return MISSING;
    }
#endif

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search.
    LayerIndex layer_idx = get_layer_index(layer, height, num_levs);

    // start with interpolated bottom layer
    float data_sum = interp_height(layer.zbot, height, data_arr, num_levs);
    float count = 1.0;

    for (int k = layer_idx.kbot; k <= layer_idx.ktop; k++) {
#ifndef NO_QC
        if ((height[k] == MISSING) || (data_arr[k] == MISSING)) {
            continue;
        }
#endif
        data_sum += data_arr[k]; 
        count += 1.0;
    }

    // finish with interpolated top layer
    data_sum += interp_height(layer.ztop, height, data_arr, num_levs);
    count += 1.0;

    return data_sum / count;
}


} // end namespace sharp


namespace sharp::exper {

} // end namespace sharp::exper



