/**
 * \file
 * \brief Routines for linear interpolation of vertical atmospheric profiles 
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */

#include <cmath>
#include <functional>

#include <SHARPlib/interp.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/algorithms.h>

namespace sharp {

float interp_height(float height_val, const float* height_arr, 
                    const float* data_arr, int num_levs) noexcept {
#ifndef NO_QC
    if (height_val == MISSING)
        return MISSING;
    // If the height value is beyond the top of the profile, 
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[num_levs-1]) || (height_val < height_arr[0]))
        return MISSING;
#endif

    auto comp = std::less<float>();
    int idx_bot = lower_bound(height_arr, num_levs, height_val, comp);
    int idx_top = upper_bound(height_arr, num_levs, height_val, comp);

    if ((height_arr[idx_bot] > height_val) && (idx_bot > 0))
        idx_bot -= 1;

#ifndef NO_QC
    for (; idx_bot > 0; --idx_bot) {
        if (data_arr[idx_bot] != MISSING)
            break;
    }

    for (; idx_top < num_levs; ++idx_top) {
        if (data_arr[idx_top] != MISSING)
            break;
    }
#endif

    float height_bot = height_arr[idx_bot];
    float height_top = height_arr[idx_top];
    float data_bot = data_arr[idx_bot];
    float data_top = data_arr[idx_top];
    
#ifndef NO_QC
    // if we didn't manage to find data, or our
    // profile data is missing, return missing
    if ((data_bot == MISSING) || (data_top == MISSING))
        return MISSING;
#endif

    // normalize the distance between values
    // to a range of 0-1 for the lerp routine
    float dz_norm = (height_val - height_bot) / (height_top - height_bot); 

    // return the linear interpolation
    return lerp(data_bot, data_top, dz_norm);
}


float interp_pressure(float pressure_val, const float* pressure_arr,
                      const float* data_arr, int num_levs) noexcept {
#ifndef NO_QC
    if (pressure_val == MISSING)
        return MISSING;
    // If the pressure value is beyond the top of the profile, 
    // or below the surface, we can't reasonably extrapolate
    if ((pressure_val < pressure_arr[num_levs-1]) || 
        (pressure_val > pressure_arr[0])) {
        
        return MISSING;
    }
#endif

    auto comp = std::greater<float>();
    int idx_bot = lower_bound(pressure_arr, num_levs, pressure_val, comp);
    int idx_top = upper_bound(pressure_arr, num_levs, pressure_val, comp);

    if ((pressure_arr[idx_bot] < pressure_val) && (idx_bot > 0))
        idx_bot -= 1;

#ifndef NO_QC
    for (; idx_bot > 0; --idx_bot) {
        if (data_arr[idx_bot] != MISSING)
            break;
    }

    for (; idx_top < num_levs; ++idx_top) {
        if (data_arr[idx_top] != MISSING)
            break;
    }
#endif

    float pressure_bot = pressure_arr[idx_bot];
    float pressure_top = pressure_arr[idx_top];
    float data_bot = data_arr[idx_bot];
    float data_top = data_arr[idx_top];
    
#ifndef NO_QC
    // if we didn't manage to find data, or our
    // profile data is missing, return missing
    if ((data_bot == MISSING) || (data_top == MISSING))
        return MISSING;
#endif

    // In order to linearly interpolate pressure properly, distance needs
    // to be calculated in log10(pressure) coordinates and normalized
    // between 0 and 1 for the lerp routine.
    float dp_norm = (std::log10(pressure_bot) - std::log10(pressure_val)) / 
                    (std::log10(pressure_bot) - std::log10(pressure_top)); 

    // return the linear interpolation
    return lerp(data_bot, data_top, dp_norm);
}


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


