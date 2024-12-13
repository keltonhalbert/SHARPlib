/**
 * \file
 * \brief Routines used to compute kinematic attributes and indices<!--
 * --> of vertical sounding profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/winds.h>

#include <cmath>

namespace sharp {

float u_component(float wind_speed, float wind_direction) {
#ifndef NO_QC
    if ((wind_direction == MISSING) || (wind_speed == MISSING)) return MISSING;
#endif
    wind_direction *= (PI / 180.0);
    return -1.0 * wind_speed * std::sin(wind_direction);
}

float v_component(float wind_speed, float wind_direction) {
#ifndef NO_QC
    if ((wind_direction == MISSING) || (wind_speed == MISSING)) return MISSING;
#endif
    wind_direction *= (PI / 180.0);
    return -1.0 * wind_speed * std::cos(wind_direction);
}

float vector_angle(float u_comp, float v_comp) {
#ifndef NO_QC
    if ((u_comp == MISSING) || (v_comp == MISSING)) return MISSING;
#endif
    if ((u_comp == 0) && (v_comp == 0)) return 0;

    float wind_direction = std::atan2(-1.0 * u_comp, -1.0 * v_comp) * (180.0 / PI);
    if (wind_direction < 0) wind_direction += 360.0;
    if (wind_direction < TOL) wind_direction = 0;
    return wind_direction;
}

float vector_magnitude(float u_comp, float v_comp) {
#ifndef NO_QC
    if ((u_comp == MISSING) || (v_comp == MISSING)) return MISSING;
#endif
    return std::sqrt((u_comp * u_comp) + (v_comp * v_comp));
}

float vector_magnitude_precise(float u_comp, float v_comp) {
#ifndef NO_QC
    if ((u_comp == MISSING) || (v_comp == MISSING)) return MISSING;
#endif
    return std::hypot(u_comp, v_comp);
}

WindVector components_to_vector(float u_comp, float v_comp) {
    const float wind_speed = vector_magnitude(u_comp, v_comp);
    const float wind_direction = vector_angle(u_comp, v_comp);

    return {wind_speed, wind_direction};
}

WindVector components_to_vector(WindComponents comp) {
    const float wind_speed = vector_magnitude(comp.u, comp.v);
    const float wind_direction = vector_angle(comp.u, comp.v);

    return {wind_speed, wind_direction};
}

WindComponents vector_to_components(float wind_speed, float wind_direction) {
#ifndef NO_QC
    if ((wind_direction == MISSING) || (wind_speed == MISSING))
        return {MISSING, MISSING};
#endif
    wind_direction *= (PI / 180.0);
    const float u_comp = -1.0 * wind_speed * std::sin(wind_direction);
    const float v_comp = -1.0 * wind_speed * std::cos(wind_direction);

    return {u_comp, v_comp};
}

WindComponents vector_to_components(WindVector vect) {
#ifndef NO_QC
    if ((vect.speed == MISSING) || (vect.direction == MISSING))
        return {MISSING, MISSING};
#endif
    const float wind_direction = vect.direction * (PI / 180.0);
    const float u_comp = -1.0 * vect.speed * std::sin(wind_direction);
    const float v_comp = -1.0 * vect.speed * std::cos(wind_direction);

    return {u_comp, v_comp};
}

WindComponents mean_wind(PressureLayer layer, const float pressure[],
                         const float u_wind[], const float v_wind[],
                         const int N, const bool weighted) {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING))
        return {MISSING, MISSING};
#endif

    if (layer.bottom > pressure[0]) layer.bottom = pressure[0];
    if (layer.top < pressure[N - 1]) layer.top = pressure[N - 1];

    float pr_lvl = layer.bottom;
    float u_sum = 0;
    float v_sum = 0;
    float weight = 0;
    while (pr_lvl >= layer.top) {
        float w = (weighted) ? pr_lvl : 1.0f;
        u_sum += interp_pressure(pr_lvl, pressure, u_wind, N) * w;
        v_sum += interp_pressure(pr_lvl, pressure, v_wind, N) * w;
        weight += w;
        pr_lvl += layer.delta;
    }

    const float mean_u = u_sum / weight;
    const float mean_v = v_sum / weight;
    return {mean_u, mean_v};
}

}  // end namespace sharp
