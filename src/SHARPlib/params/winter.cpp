/**
 * \file
 * \brief Winter weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Prediction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/params/winter.h>

#include <cstddef>

namespace sharp {

PressureLayer dendritic_layer(const float pressure[], const float temperature[],
                              const std::ptrdiff_t N) {
    float dgz_bot = MISSING;
    float dgz_top = MISSING;
    constexpr float minus_17 = -17.0f + ZEROCNK;
    constexpr float minus_12 = -12.0f + ZEROCNK;

    auto interpolate_log_p = [](float p1, float p2, float w) {
        if (w <= 0.0f) return p1;
        if (w >= 1.0f) return p2;
        return p1 * std::pow(p2 / p1, w);
    };

    float last_p = MISSING;
    float last_t = MISSING;
    std::ptrdiff_t idx = N - 1;
    while (idx >= 0) {
        last_p = pressure[idx];
        last_t = temperature[idx];
        --idx;
        if (last_p != MISSING && last_t != MISSING) break;
    }

    for (; idx >= 0; --idx) {
        float p_curr = pressure[idx];
        float t_curr = temperature[idx];

        if (p_curr == MISSING || t_curr == MISSING) continue;

        float t1 = last_t, t2 = t_curr;
        float p1 = last_p, p2 = p_curr;
        float w_in = 0.0f;
        float w_out = 1.0f;
        float dt = t2 - t1;

        if (std::abs(dt) > 1e-5f) {
            float inv_dt = 1.0f / dt;
            float w17 = (minus_17 - t1) * inv_dt;
            float w12 = (minus_12 - t1) * inv_dt;
            w_in = std::max(0.0f, std::min(w17, w12));
            w_out = std::min(1.0f, std::max(w17, w12));
        } else {
            if (t1 >= minus_17 && t1 <= minus_12) {
                w_in = 0.0f;
                w_out = 1.0f;
            } else {
                w_in = 1.0f;
                w_out = 0.0f;
            }
        }

        if (w_in < w_out) {
            if (dgz_top == MISSING) {
                dgz_top = interpolate_log_p(p1, p2, w_in);
            }

            dgz_bot = interpolate_log_p(p1, p2, w_out);

            if (w_out < 1.0f) {
                break;
            }
        } else if (dgz_top != MISSING) {
            break;
        }

        last_p = p_curr;
        last_t = t_curr;
    }

    return {dgz_bot, dgz_top};
}

float snow_squall_parameter(const float wetbulb_2m, const float mean_relh_0_2km,
                            const float delta_thetae_0_2km,
                            const float mean_wind_0_2km) {
    if (wetbulb_2m > ZEROCNK + 1) return 0.0f;
    const float relh_term = (mean_relh_0_2km - 0.60f) / 0.15f;
    const float thetae_term = (4.0f - delta_thetae_0_2km) / 4.0f;
    if ((relh_term < 0) || (thetae_term < 0)) return 0.0f;

    return relh_term * thetae_term * (mean_wind_0_2km / 9.0f);
}

}  // namespace sharp
