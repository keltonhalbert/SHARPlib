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
    float wind_speed = vector_magnitude(comp.u, comp.v);
    float wind_direction = vector_angle(comp.u, comp.v);

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


WindComponents vector_to_components(const WindVector& vect) {
#ifndef NO_QC
    if ((vect.speed == MISSING) || (vect.direction == MISSING))
        return {MISSING, MISSING}; 
#endif
    float wind_direction = vect.direction * (PI / 180.0);
    float u_comp = -1.0 * vect.speed * std::sin(wind_direction);
    float v_comp = -1.0 * vect.speed * std::cos(wind_direction);

    return {u_comp, v_comp};
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


WindComponents wind_shear(const PressureLayer& layer, 
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

    float u_bot = interp_pressure(pbot, pres, u_wind, num_levs);
    float u_top = interp_pressure(ptop, pres, u_wind, num_levs);

    float v_bot = interp_pressure(pbot, pres, v_wind, num_levs);
    float v_top = interp_pressure(ptop, pres, v_wind, num_levs);

#ifndef NO_QC
    if ((u_bot == MISSING) || (v_bot == MISSING) ||
        (u_top == MISSING) || (v_top == MISSING)) {

        return {MISSING, MISSING}; 
    }
#endif

    return {u_top - u_bot, v_top - v_bot};

}


// TO-DO: Right now it is unclear whether the base usings of
// u_wind and v_wind should be knots or m/s. Need to clarify
// design here and update documentation accordingly. 
float helicity(const HeightLayer& layer_agl, 
               const WindComponents& storm_motion,
               const float* height, const float* u_wind, 
               const float* v_wind, int num_levs) {
#ifndef NO_QC
    if ((storm_motion.u == MISSING) || (storm_motion.v == MISSING)) {
        return MISSING;
    }
    if ((layer_agl.zbot == MISSING) || (layer_agl.ztop == MISSING)) {
        return MISSING;
    }
#endif

    // get the height in MSL by adding the surface height
    float hbot = height[0] + layer_agl.zbot;
    float htop = height[0] + layer_agl.ztop;

    // bounds check the integration
    if (hbot < height[0]) {
        hbot = height[0];
    }
    if (htop > height[num_levs-1]) {
        htop = height[num_levs-1];
    }

    // find the lowest observation, but not if it
    // is exactly the bottom of our layer (it gets
    // interpolated).
    int lower_idx = 0;
    while (height[lower_idx] <= hbot) {
        lower_idx++;
    }

    // find the highest observation, but not if it
    // is exactly the top of our layer (it gets 
    // interpolated)
    int upper_idx = num_levs - 1;
    while (height[upper_idx] >= htop) {
        upper_idx--;
    }

    // Get the interpolated bottom of the layer to integrate
    // and convert to storm-relative winds
    float sru_bot = interp_height(hbot, height, u_wind, num_levs);
    float srv_bot = interp_height(hbot, height, v_wind, num_levs);
    sru_bot -= storm_motion.u;
    srv_bot -= storm_motion.v;

    // will get set in first loop iter
    float sru_top;
    float srv_top;
    float lyrh = 0.0;
    float phel = 0.0;
    float nhel = 0.0;
    for (int k = lower_idx; k <= upper_idx; k++) {
#ifndef NO_QC
        if ((u_wind[k] != MISSING) && (v_wind[k] != MISSING)) {
            continue;
        }
#endif
        // top of layer storm relative winds
        sru_top = u_wind[k] - storm_motion.u;
        srv_top = v_wind[k] - storm_motion.v;

        // integrate layer
        lyrh = (sru_top * srv_bot) - (sru_bot * srv_top); 
        if (lyrh > 0.0) {
            phel += lyrh;
        }
        else {
            nhel += lyrh;
        }
        // set the top to be the bottom
        // of the next layer
        sru_bot = sru_top;
        srv_bot = srv_top;
    }

    // Get the interpolated top of the layer to integrate
    // and convert to storm-relative winds
    sru_top = interp_height(htop, height, u_wind, num_levs);
    srv_top = interp_height(htop, height, v_wind, num_levs);
    sru_top -= storm_motion.u;
    srv_top -= storm_motion.v;
    
    // integrate the final layer
    lyrh = (sru_top * srv_bot) - (sru_bot * srv_top); 
    if (lyrh > 0.0) {
        phel += lyrh;
    }
    else {
        nhel += lyrh;
    }

    return (phel + nhel);
}


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


