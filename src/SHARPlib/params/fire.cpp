/**
 * \file
 * \brief Fire weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params/fire.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>
#define FMT_HEADER_ONLY
#include <fmt/core.h>

#include <algorithm>
#include <cmath>

namespace sharp {

float equilibrium_moisture_content(float temperature, float rel_humidity) {
#ifndef NO_QC
    if ((temperature == MISSING) || (rel_humidity == MISSING)) {
        return MISSING;
    }
#endif

    rel_humidity = rel_humidity * 100;
    temperature = (temperature - ZEROCNK) * (9.f / 5.f) + 32.f;

    float emc;
    if (rel_humidity < 10.f) {
        emc = 0.03229f + (0.281073f * rel_humidity) -
              (0.000578f * rel_humidity * temperature);

    } else if ((rel_humidity >= 10.f) && (rel_humidity <= 50.f)) {
        emc = 2.22749f + (0.160107f * rel_humidity) - (0.01478f * temperature);

    } else {
        emc = 21.0606f + (0.005565f * (rel_humidity * rel_humidity)) -
              (0.00035f * rel_humidity * temperature) -
              (0.483199f * rel_humidity);
    }

    return emc / 100.0f;
}

float fosberg_fire_index(float temperature, float rel_humidity,
                         float wind_speed) {
#ifndef NO_QC
    if ((temperature == MISSING) || (rel_humidity == MISSING) ||
        (wind_speed == MISSING)) {
        return MISSING;
    }
#endif
    // m/s to mph
    wind_speed = wind_speed * 2.237;
    float emc =
        equilibrium_moisture_content(temperature, rel_humidity) * 100.0f;

    emc = emc / 30.0f;
    const float eta = 1.0f - (2.0f * emc) + (1.5f * (std::pow(emc, 2))) -
                      (0.5f * (std::pow(emc, 3)));

    const float fwwi =
        eta * std::sqrt(1.0f + std::pow(wind_speed, 2)) / 0.3002f;
    return std::min(fwwi, 100.f);
}

float pft_plume_potential_temperature(float theta_env, float beta) {
#ifndef NO_QC
    if (theta_env == MISSING) {
        return MISSING;
    }
#endif
    return (beta + 1.0f) * theta_env;
}

float pft_plume_mixratio(float theta_env, float mixr_env, float beta,
                         float phi) {
#ifndef NO_QC
    if ((theta_env == MISSING) || (mixr_env == MISSING)) {
        return MISSING;
    }
#endif
    float spfh_env = specific_humidity(mixr_env);
    float spfh_plume = std::max(0.0f, spfh_env + (beta * phi * theta_env));
    return mixratio(spfh_plume);
}

}  // namespace sharp
