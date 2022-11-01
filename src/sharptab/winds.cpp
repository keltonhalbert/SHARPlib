// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
#include <cmath>

#include "constants.h"
#include "winds.h"

namespace sharp {


float u_component(float wind_speed, float wind_direction) {
#ifndef NO_QC
    if ((wind_direction == MISSING) || (wind_speed == MISSING))
        return MISSING;
#endif
    wind_direction *= (PI / 180.0);
    return -1.0 * wind_speed * std::sin(wind_direction);
}


float v_component(float wind_speed, float wind_direction) {
#ifndef NO_QC
    if ((wind_direction == MISSING) || (wind_speed == MISSING))
        return MISSING;
#endif
    wind_direction *= (PI / 180.0);
    return -1.0 * wind_speed * std::cos(wind_direction);
}


float vector_angle(float u_comp, float v_comp) {
#ifndef NO_QC
    if ((u_comp == MISSING) || (v_comp == MISSING))
        return MISSING;
#endif
    if ((u_comp == 0) && (v_comp == 0))
        return 0;

    float wind_direction = std::atan2( -1*u_comp, -1*v_comp) * (180.0 / PI);
    if (wind_direction < 0) wind_direction += 360.0;
    if (wind_direction < TOL) wind_direction = 0;
    return wind_direction; 
}


float vector_magnitude(float u_comp, float v_comp) {
#ifndef NO_QC
    if ((u_comp == MISSING) || (v_comp == MISSING))
        return MISSING;
#endif
    return std::sqrt((u_comp * u_comp) + (v_comp * v_comp));
}


float vector_magnitude_precise(float u_comp, float v_comp) {
#ifndef NO_QC
    if ((u_comp == MISSING) || (v_comp == MISSING))
        return MISSING;
#endif
    return std::hypot(u_comp, v_comp);
}


WindVector components_to_vector(float u_comp, float v_comp) {
    float wind_speed = vector_magnitude(u_comp, v_comp);
    float wind_direction = vector_angle(u_comp, v_comp);

    return {wind_speed, wind_direction};
}


WindComponents vector_to_components(float wind_speed, float wind_direction) {
#ifndef NO_QC
    if ((wind_direction == MISSING) || (wind_speed == MISSING))
        return {MISSING, MISSING}; 
#endif
    wind_direction *= (PI / 180.0);
    float u_comp = -1.0 * wind_speed * std::sin(wind_direction);
    float v_comp = -1.0 * wind_speed * std::cos(wind_direction);

    return {u_comp, v_comp};
}


}
