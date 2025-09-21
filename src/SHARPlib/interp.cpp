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

#include <cmath>
#include <cstddef>
#include <functional>

namespace sharp {

float interp_height(const float height_val, const float height_arr[],
                    const float data_arr[], const std::ptrdiff_t N) {
#ifndef NO_QC
    if (height_val == MISSING) return MISSING;
    if (std::isnan(height_val)) return MISSING;
#endif
    // If the height value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[N - 1]) || (height_val < height_arr[0]))
        return MISSING;

    static constexpr auto comp = std::less<float>();
    std::ptrdiff_t idx_top = upper_bound(height_arr, N, height_val, comp);
    std::ptrdiff_t idx_bot = idx_top - 1;

#ifndef NO_QC
    for (; idx_bot > 0; --idx_bot) {
        if (data_arr[idx_bot] != MISSING) break;
    }

    for (; idx_top < N - 1; ++idx_top) {
        if (data_arr[idx_top] != MISSING) break;
    }

    // in the case the data are still missing at this point,
    // return a missing value
    if ((data_arr[idx_bot] == MISSING) || (data_arr[idx_top] == MISSING))
        return MISSING;
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
#ifndef NO_QC
    if (pressure_val == MISSING) return MISSING;
    if (std::isnan(pressure_val)) return MISSING;
#endif
    // If the pressure value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((pressure_val < pressure_arr[N - 1]) ||
        (pressure_val > pressure_arr[0])) {
        return MISSING;
    }

    static constexpr auto comp = std::greater<float>();
    std::ptrdiff_t idx_top = upper_bound(pressure_arr, N, pressure_val, comp);
    std::ptrdiff_t idx_bot = idx_top - 1;

#ifndef NO_QC
    for (; idx_bot > 0; --idx_bot) {
        if (data_arr[idx_bot] != MISSING) break;
    }

    for (; idx_top < N - 1; ++idx_top) {
        if (data_arr[idx_top] != MISSING) break;
    }

    // in the case the data are still missing at this point,
    // return a missing value
    if ((data_arr[idx_bot] == MISSING) || (data_arr[idx_top] == MISSING))
        return MISSING;
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
    if (data_val == MISSING) {
        return MISSING;
    }
    for (; k_start < N; ++k_start) {
        if (data_arr[k_start] != MISSING) break;
    }
#endif

    for (std::ptrdiff_t k = k_start + 1; k < N; ++k) {
        float val0 = data_arr[k_start];
        float val1 = data_arr[k];
#ifndef NO_QC
        if (val1 == MISSING) continue;
#endif
        if (val0 == data_val) return pressure_arr[k_start];
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
    if (data_val == MISSING) {
        return MISSING;
    }
    for (; k_start < N; ++k_start) {
        if (data_arr[k_start] != MISSING) break;
    }
#endif

    for (std::ptrdiff_t k = k_start + 1; k < N; ++k) {
        float val0 = data_arr[k_start];
        float val1 = data_arr[k];
#ifndef NO_QC
        if (val1 == MISSING) continue;
#endif
        if (val0 == data_val) return height_arr[k_start];
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
