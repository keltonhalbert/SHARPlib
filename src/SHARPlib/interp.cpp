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
#include <functional>

namespace sharp {

float interp_height(float height_val, const float height_arr[],
                    const float data_arr[], const std::ptrdiff_t N) {
#ifndef NO_QC
    if (height_val == MISSING) return MISSING;
    // If the height value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[N - 1]) || (height_val < height_arr[0]))
        return MISSING;
#endif

    static constexpr auto comp = std::less<float>();
    std::ptrdiff_t idx_top = upper_bound(height_arr, N, height_val, comp);
    std::ptrdiff_t idx_bot = idx_top - 1;

#ifndef NO_QC
    for (; idx_bot > 0; --idx_bot) {
        if (data_arr[idx_bot] != MISSING) break;
    }

    for (; idx_top < N; ++idx_top) {
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

float interp_pressure(float pressure_val, const float pressure_arr[],
                      const float data_arr[], const std::ptrdiff_t N) {
#ifndef NO_QC
    if (pressure_val == MISSING) return MISSING;
    // If the pressure value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((pressure_val < pressure_arr[N - 1]) ||
        (pressure_val > pressure_arr[0])) {
        return MISSING;
    }
#endif

    static constexpr auto comp = std::greater<float>();
    std::ptrdiff_t idx_top = upper_bound(pressure_arr, N, pressure_val, comp);
    std::ptrdiff_t idx_bot = idx_top - 1;

#ifndef NO_QC
    for (; idx_bot > 0; --idx_bot) {
        if (data_arr[idx_bot] != MISSING) break;
    }

    for (; idx_top < N; ++idx_top) {
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

}  // end namespace sharp
