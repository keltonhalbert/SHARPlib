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

#include "interp.h"
#include "constants.h"

namespace sharp {

float interp_height(float height_val, const float* height_arr, 
                    const float* data_arr, int num_levs) {

    // If the height value is beyond the top of the profile, 
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[num_levs-1]) || (height_val < height_arr[0]))
        return MISSING;

    // Find the height levels that our request lies between
    float height_bot = MISSING;
    float height_top = MISSING;
    float data_bot = MISSING;
    float data_top = MISSING;
    for (int k = 0; k < num_levs-1; k++) {
        height_bot = height_arr[k];
        height_top = height_arr[k+1];
        
        // If we have found our height levels, 
        // load the data levels and exit the loop
        if ((height_bot <= height_val) && (height_top >= height_val)) {
            data_bot = data_arr[k];
            data_top = data_arr[k+1];
            break;
        }
    }

    // if we didn't manage to find data, or our
    // profile data is missing, return missing
    if ((data_bot == MISSING) || (data_top == MISSING))
        return MISSING;

    // normalize the distance between values
    // to a range of 0-1 for the lerp routine
    float dz_norm = (height_val - height_bot) / (height_top - height_bot); 

    // return the linear interpolation
    return lerp(data_bot, data_top, dz_norm);
}


float interp_pressure(float pressure_val, const float* pressure_arr,
                      const float* data_arr, int num_levs) {

    // If the pressure value is beyond the top of the profile, 
    // or below the surface, we can't reasonably extrapolate
    if ((pressure_val < pressure_arr[num_levs-1]) || 
        (pressure_val > pressure_arr[0])) {
        
        return MISSING;
    }

    // Find the pressure levels that our request lies between
    float pressure_bot = MISSING;
    float pressure_top = MISSING;
    float data_bot = MISSING;
    float data_top = MISSING;
    for (int k = 0; k < num_levs-1; k++) {
        pressure_bot = pressure_arr[k];
        pressure_top = pressure_arr[k+1];
        
        // If we have found our pressure levels, 
        // load the data levels and exit the loop
        if ((pressure_bot >= pressure_val) && (pressure_top <= pressure_val)) {
            data_bot = data_arr[k];
            data_top = data_arr[k+1];
            break;
        }
    }

    // if we didn't manage to find data, or our
    // profile data is missing, return missing
    if ((data_bot == MISSING) || (data_top == MISSING))
        return MISSING;

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


