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
#include <cstddef>

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

    float wind_direction =
        std::atan2(-1.0 * u_comp, -1.0 * v_comp) * (180.0 / PI);
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
                         const std::ptrdiff_t N, const bool weighted) {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING))
        return {MISSING, MISSING};
#endif

    if (layer.bottom > pressure[0]) layer.bottom = pressure[0];
    if (layer.top < pressure[N - 1]) layer.top = pressure[N - 1];

    LayerIndex idx = get_layer_index(layer, pressure, N);

    float u_lyr_bot = interp_pressure(layer.bottom, pressure, u_wind, N);
    float v_lyr_bot = interp_pressure(layer.bottom, pressure, v_wind, N);
    float u_lyr_top = interp_pressure(layer.top, pressure, u_wind, N);
    float v_lyr_top = interp_pressure(layer.top, pressure, v_wind, N);

#ifndef NO_QC
    if ((u_lyr_bot == MISSING) || (v_lyr_bot == MISSING) ||
        (u_lyr_top == MISSING) || (v_lyr_top == MISSING)) {
        return {MISSING, MISSING};
    }
#endif

    float u_integ = 0.0f;
    float v_integ = 0.0f;
    float w_integ = 0.0f;

    float p_prev = layer.bottom;
    float u_prev = u_lyr_bot;
    float v_prev = v_lyr_bot;

    for (std::ptrdiff_t k = idx.kbot; k <= idx.ktop; ++k) {
#ifndef NO_QC
        if ((u_wind[k] == MISSING) || (v_wind[k] == MISSING)) {
            continue;
        }
#endif
        const float p_curr = pressure[k];
        const float u_curr = u_wind[k];
        const float v_curr = v_wind[k];
        const float dp = p_prev - p_curr;

        if (weighted) {
            u_integ += 0.5f * (u_prev * p_prev + u_curr * p_curr) * dp;
            v_integ += 0.5f * (v_prev * p_prev + v_curr * p_curr) * dp;
            w_integ += 0.5f * (p_prev + p_curr) * dp;
        } else {
            u_integ += 0.5f * (u_prev + u_curr) * dp;
            v_integ += 0.5f * (v_prev + v_curr) * dp;
            w_integ += dp;
        }

        p_prev = p_curr;
        u_prev = u_curr;
        v_prev = v_curr;
    }

    {
        const float dp = p_prev - layer.top;

        if (weighted) {
            u_integ += 0.5f * (u_prev * p_prev + u_lyr_top * layer.top) * dp;
            v_integ += 0.5f * (v_prev * p_prev + v_lyr_top * layer.top) * dp;
            w_integ += 0.5f * (p_prev + layer.top) * dp;
        } else {
            u_integ += 0.5f * (u_prev + u_lyr_top) * dp;
            v_integ += 0.5f * (v_prev + v_lyr_top) * dp;
            w_integ += dp;
        }
    }

    return {u_integ / w_integ, v_integ / w_integ};
}

/// @cond DOXYGEN_IGNORE

template WindComponents max_wind<PressureLayer>(PressureLayer lyr,
                                                const float coordinate[],
                                                const float u_wind[],
                                                const float v_wind[],
                                                const std::ptrdiff_t N);

template WindComponents max_wind<HeightLayer>(HeightLayer lyr,
                                              const float coordinate[],
                                              const float u_wind[],
                                              const float v_wind[],
                                              const std::ptrdiff_t N);

template WindComponents wind_shear<PressureLayer>(PressureLayer layer,
                                                  const float coord[],
                                                  const float u_wind[],
                                                  const float v_wind[],
                                                  const std::ptrdiff_t N);

template WindComponents wind_shear<HeightLayer>(HeightLayer layer,
                                                const float coord[],
                                                const float u_wind[],
                                                const float v_wind[],
                                                const std::ptrdiff_t N);

template float helicity<PressureLayer>(
    PressureLayer layer, WindComponents storm_motion, const float coord[],
    const float u_wind[], const float v_wind[], const std::ptrdiff_t N);

template float helicity<HeightLayer>(HeightLayer layer,
                                     WindComponents storm_motion,
                                     const float coord[], const float u_wind[],
                                     const float v_wind[],
                                     const std::ptrdiff_t N);

/// @endcond

}  // end namespace sharp
