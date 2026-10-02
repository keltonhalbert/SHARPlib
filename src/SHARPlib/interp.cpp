/**
 * \file
 * \brief Routines for linear interpolation of vertical atmospheric profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/qc.h>

#include <cmath>
#include <cstddef>
#include <functional>

namespace sharp {

#ifndef NO_QC
static bool widen_bracket_past_missing(const float coord_val,
                                       const float coord_arr[],
                                       const float data_arr[],
                                       const std::ptrdiff_t N,
                                       std::ptrdiff_t& idx_bot,
                                       std::ptrdiff_t& idx_top, float& result) {
    const bool bot_missing = is_missing(data_arr[idx_bot]);
    const bool top_missing = is_missing(data_arr[idx_top]);
    if (!bot_missing && !top_missing) return false;

    if (!bot_missing && (coord_val == coord_arr[idx_bot])) {
        result = data_arr[idx_bot];
        return true;
    }
    if (!top_missing && (coord_val == coord_arr[idx_top])) {
        result = data_arr[idx_top];
        return true;
    }

    for (; idx_bot > 0; --idx_bot) {
        if (!is_missing(data_arr[idx_bot])) break;
    }

    for (; idx_top < N - 1; ++idx_top) {
        if (!is_missing(data_arr[idx_top])) break;
    }

    // in the case the data are still missing at this point,
    // return a missing value
    if (is_missing(data_arr[idx_bot]) || is_missing(data_arr[idx_top])) {
        result = MISSING;
        return true;
    }
    return false;
}
#endif

float interp_height(const float height_val, const float height_arr[],
                    const float data_arr[], const std::ptrdiff_t N) {
    if (N < 1) return MISSING;
#ifndef NO_QC
    if (is_missing(height_val)) return MISSING;
#endif
    // If the height value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[N - 1]) || (height_val < height_arr[0]))
        return MISSING;

    if (N == 1) {
#ifndef NO_QC
        if (is_missing(data_arr[0])) return MISSING;
#endif
        return data_arr[0];
    }

    static constexpr auto comp = std::less<float>();
    std::ptrdiff_t idx_top = upper_bound(height_arr, N, height_val, comp);
    std::ptrdiff_t idx_bot = idx_top - 1;

#ifndef NO_QC
    float bridged = MISSING;
    if (widen_bracket_past_missing(height_val, height_arr, data_arr, N, idx_bot,
                                   idx_top, bridged))
        return bridged;
#endif

    const float height_bot = height_arr[idx_bot];
    const float height_top = height_arr[idx_top];
    const float data_bot = data_arr[idx_bot];
    const float data_top = data_arr[idx_top];

    // normalize the distance between values
    // to a range of 0-1 for the lerp routine
    const float dz_norm = (height_val - height_bot) / (height_top - height_bot);

    // return the linear interpolation
    return lerp(data_bot, data_top, dz_norm);
}

float interp_pressure(const float pressure_val, const float pressure_arr[],
                      const float data_arr[], const std::ptrdiff_t N) {
    if (N < 1) return MISSING;
#ifndef NO_QC
    if (is_missing(pressure_val)) return MISSING;
#endif
    // If the pressure value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((pressure_val < pressure_arr[N - 1]) ||
        (pressure_val > pressure_arr[0])) {
        return MISSING;
    }

    if (N == 1) {
#ifndef NO_QC
        if (is_missing(data_arr[0])) return MISSING;
#endif
        return data_arr[0];
    }

    static constexpr auto comp = std::greater<float>();
    std::ptrdiff_t idx_top = upper_bound(pressure_arr, N, pressure_val, comp);
    std::ptrdiff_t idx_bot = idx_top - 1;

#ifndef NO_QC
    float bridged = MISSING;
    if (widen_bracket_past_missing(pressure_val, pressure_arr, data_arr, N,
                                   idx_bot, idx_top, bridged))
        return bridged;
#endif

    const float pressure_bot = pressure_arr[idx_bot];
    const float pressure_top = pressure_arr[idx_top];
    const float data_bot = data_arr[idx_bot];
    const float data_top = data_arr[idx_top];

    // In order to linearly interpolate pressure properly, distance needs
    // to be calculated in log10(pressure) coordinates and normalized
    // between 0 and 1 for the lerp routine.
    const float dp_norm =
        (std::log10(pressure_bot) - std::log10(pressure_val)) /
        (std::log10(pressure_bot) - std::log10(pressure_top));

    // return the linear interpolation
    return lerp(data_bot, data_top, dp_norm);
}

float find_first_pressure(const float data_val, const float pressure_arr[],
                          const float data_arr[], const std::ptrdiff_t N) {
    std::ptrdiff_t k_start = 0;
#ifndef NO_QC
    if (is_missing(data_val)) {
        return MISSING;
    }
    for (; k_start < N; ++k_start) {
        if (!is_missing(data_arr[k_start])) break;
    }
#endif
    if ((k_start < N) && (data_arr[k_start] == data_val))
        return pressure_arr[k_start];

    for (std::ptrdiff_t k = k_start + 1; k < N; ++k) {
        float val0 = data_arr[k_start];
        float val1 = data_arr[k];
#ifndef NO_QC
        if (is_missing(val1)) continue;
#endif
        if (val1 == data_val) return pressure_arr[k];

        if ((data_val - val0) * (data_val - val1) < 0) {
            const float logp_bot = std::log10(pressure_arr[k_start]);
            const float logp_top = std::log10(pressure_arr[k]);

            const float d_norm = (data_val - val0) / (val1 - val0);

            return std::pow(10, lerp(logp_bot, logp_top, d_norm));
        }

        k_start = k;
    }

    return MISSING;
}

float find_first_height(const float data_val, const float height_arr[],
                        const float data_arr[], const std::ptrdiff_t N) {
    std::ptrdiff_t k_start = 0;
#ifndef NO_QC
    if (is_missing(data_val)) {
        return MISSING;
    }
    for (; k_start < N; ++k_start) {
        if (!is_missing(data_arr[k_start])) break;
    }
#endif
    if ((k_start < N) && (data_arr[k_start] == data_val))
        return height_arr[k_start];

    for (std::ptrdiff_t k = k_start + 1; k < N; ++k) {
        float val0 = data_arr[k_start];
        float val1 = data_arr[k];
#ifndef NO_QC
        if (is_missing(val1)) continue;
#endif
        if (val1 == data_val) return height_arr[k];

        // will have a negative sign if levels straddle
        // the point being searched for
        if ((data_val - val0) * (data_val - val1) < 0) {
            const float hght_bot = height_arr[k_start];
            const float hght_top = height_arr[k];

            const float d_norm = (data_val - val0) / (val1 - val0);

            return lerp(hght_bot, hght_top, d_norm);
        }

        k_start = k;
    }

    return MISSING;
}

}  // end namespace sharp
