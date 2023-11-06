/**
 * \file
 * \brief Routines for linear interpolation of vertical atmospheric profiles
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

#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/thermo.h>

#include <cmath>
#include <functional>

namespace sharp {

float interp_height(float height_val, const float height_arr[],
                    const float data_arr[], const int N) noexcept {
#ifndef NO_QC
    if (height_val == MISSING) return MISSING;
    // If the height value is beyond the top of the profile,
    // or below the surface, we can't reasonably extrapolate
    if ((height_val > height_arr[N - 1]) || (height_val < height_arr[0]))
        return MISSING;
#endif

    static constexpr auto comp = std::less<float>();
    int idx_top = upper_bound(height_arr, N, height_val, comp);
    int idx_bot = idx_top - 1;

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
                      const float data_arr[], const int N) noexcept {
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
    int idx_top = upper_bound(pressure_arr, N, pressure_val, comp);
    int idx_bot = idx_top - 1;

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

float interp_hghtlevel(float data_val, const float data_arr[],
                      const float height_arr[], const int N) noexcept {
#ifndef NO_QC
    if (data_val == MISSING) return MISSING;
#endif

    static constexpr auto comp = std::greater<float>();
    int idx_top = lower_bound(data_arr, N, data_val, comp);
    int idx_bot = idx_top - 1;

#ifndef NO_QC
    if (idx_bot < 0) {
        return MISSING;
    }
    for (; idx_bot > 0; --idx_bot) {
        if (height_arr[idx_bot] != MISSING) break;
    }

    for (; idx_top < N; ++idx_top) {
        if (height_arr[idx_top] != MISSING) break;
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
    const float dv_norm = (data_val - data_bot) / (data_top - data_bot);

    // return the linear interpolation
    return lerp(height_bot, height_top, dv_norm);
}

float interp_wbzh(const float p_arr[], const float t_arr[], const float td_arr[],
                      const float height_arr[], const int N) noexcept {

    int i;
    float temwet1, temwet2;
    float dv_norm = MISSING;

    temwet2 = wetbulb(p_arr[0],t_arr[0],td_arr[0]);
    for (i=1; i<N; i++) {
	temwet1 = temwet2;
	temwet2 = wetbulb(p_arr[i],t_arr[i],td_arr[i]);
	if(temwet1*temwet2<0) {
		dv_norm = (0.0 - temwet1) / (temwet2 - temwet1);
		return lerp(height_arr[i-1],height_arr[i],dv_norm);
	}
    }

    // in the case the data are still missing at this point,
    // return a missing value
    return MISSING;
}

}  // end namespace sharp

