// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 

#include "interp.h"
#include "constants.h"

namespace sharp {

float interp_height(const float& height_val, const float* height_arr, 
                    const float* data_arr, const int& num_levs) {

    // If the height value is beyond the top of the profile, 
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[num_levs-1]) || (height_val < height_arr[0]))
        return MISSING;

    // Find the height levels that our request lies between
    float height_lo = MISSING;
    float height_hi = MISSING;
    float data_lo = MISSING;
    float data_hi = MISSING;
    for (int k = 0; k < num_levs-1; k++) {
        height_lo = height_arr[k];
        height_hi = height_arr[k+1];
        
        // If we have found our height levels, 
        // load the data levels and exit the loop
        if ((height_lo <= height_val) && (height_hi >= height_val)) {
            data_lo = data_arr[k];
            data_hi = data_arr[k+1];
            break;
        }
    }

    // if we didn't manage to find data, or our
    // profile data is missing, return missing
    if ((data_lo == MISSING) || (data_hi == MISSING))
        return MISSING;

    // normalize the distance between values
    // to a range of 0-1 for the lerp routine
    float dz_norm = (height_val - height_lo) / (height_hi - height_lo); 

    // return the linear interpolation
    return lerp(data_lo, data_hi, dz_norm);

}

/*
float interp_pressure() {

}
*/

}
