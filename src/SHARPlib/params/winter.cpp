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
    float bottom =
        find_first_pressure(-12.0f + ZEROCNK, pressure, temperature, N);
    float top = find_first_pressure(-17.0f + ZEROCNK, pressure, temperature, N);

    if (top == MISSING) {
        return {MISSING, MISSING};
    }

    if (bottom == MISSING) {
        bottom = pressure[0];
    }

    return {bottom, top};
}

float snow_squall_parameter(const float mean_relh_0_2km,
                            const float delta_thetae_0_2km,
                            const float mean_wind_0_2km) {
    return ((mean_relh_0_2km - .60f) / .15f) *
           ((4.0f - delta_thetae_0_2km) / 4.0f) * (mean_wind_0_2km / 9.0f);
}

}  // namespace sharp
