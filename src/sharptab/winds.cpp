/**
 * \file
 * \brief Routines used to compute kinematic attributes and indices of vertical sounding profiles
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

#include "constants.h"
#include "interp.h"
#include "utils.h"
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


WindVector components_to_vector(const WindComponents& comp) {
    WindVector vec = components_to_vector(comp.u, comp.v); 
    return {vec.speed, vec.direction};
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


WindComponents vector_to_components(const WindVector& vect) {
    WindComponents comp = vector_to_components(vect.speed, vect.direction); 
    return {comp.u, comp.v};
}


WindComponents mean_wind(const PressureLayer& layer, 
                         const float* pres,   const float* u_wind, 
                         const float* v_wind, int num_levs) {
#ifndef NO_QC
    if ((layer.pbot == MISSING) || (layer.ptop == MISSING))
        return {MISSING, MISSING};
#endif

    // bounds checking our layer
    // requests to the available data
    float pbot, ptop;
    if (layer.pbot > pres[0])
        pbot = pres[0];
    else
        pbot = layer.pbot;

    if (layer.ptop < pres[num_levs-1])
        ptop = pres[num_levs-1];
    else
        ptop = layer.ptop;


    float pres_level = pbot;
    float u_sum = 0;
    float v_sum = 0;
    float weight = 0;
    while(pres_level >= ptop) {
        u_sum += interp_pressure(pres_level, pres, u_wind, num_levs) * pres_level; 
        v_sum += interp_pressure(pres_level, pres, v_wind, num_levs) * pres_level; 
        weight += pres_level;
        pres_level -= layer.dp;
    }

    float mean_u = u_sum / weight;
    float mean_v = v_sum / weight;

    return {mean_u, mean_v}; 
}


WindComponents mean_wind_npw(const PressureLayer& layer,
                             const float* pres,   const float* u_wind,
                             const float* v_wind, int num_levs) {
#ifndef NO_QC
    if ((layer.pbot == MISSING) || (layer.ptop == MISSING))
        return {MISSING, MISSING};
#endif

    // bounds checking our layer
    // requests to the available data
    float pbot, ptop;
    if (layer.pbot > pres[0])
        pbot = pres[0];
    else
        pbot = layer.pbot;

    if (layer.ptop < pres[num_levs-1])
        ptop = pres[num_levs-1];
    else
        ptop = layer.ptop;


    float pres_level = pbot;
    float u_sum = 0;
    float v_sum = 0;
    float weight = 0;
    while(pres_level >= ptop) {
        u_sum += interp_pressure(pres_level, pres, u_wind, num_levs);
        v_sum += interp_pressure(pres_level, pres, v_wind, num_levs);
        weight += 1; 
        pres_level -= layer.dp;
    }

    float mean_u = u_sum / weight;
    float mean_v = v_sum / weight;

    return {mean_u, mean_v}; 

}


WindComponents mean_wind(float pressure_bot,  float pressure_top,
                         const float* pres,   const float* u_wind, 
                         const float* v_wind, int num_levs) {

    PressureLayer layer = {pressure_bot, pressure_top, 1};
    return mean_wind(layer, pres, u_wind, v_wind, num_levs);

}


WindComponents mean_wind_nwp(float pressure_bot,  float pressure_top,
                             const float* pres,   const float* u_wind, 
                             const float* v_wind, int num_levs) {

    PressureLayer layer = {pressure_bot, pressure_top, 1};
    return mean_wind_npw(layer, pres, u_wind, v_wind, num_levs);

}


}
