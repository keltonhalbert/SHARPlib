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
#include <SHARPlib/params/fire.h>

#include <cmath>

namespace sharp {

float equilibrium_moisture_content(float temperature, float rel_humidity) {
#ifndef NO_QC
    if ((temperature == MISSING) || (rel_humidity)) {
        return MISSING;
    }
#endif

    float emc;
    if (rel_humidity < 0.1) {
        emc = 0.03229f + 0.281073f * (rel_humidity * 100.0f) -
              0.000578f * (100.0f * rel_humidity) * temperature;

    } else if ((rel_humidity <= 0.1f) && (rel_humidity <= 0.5f)) {
        emc = 2.22749f + 0.160107f * (rel_humidity * 100.0f) -
              0.01478f * temperature;

    } else {
        emc = 21.0606f +
              0.005565f * ((rel_humidity * 100.0f) * (rel_humidity * 100.0f)) -
              0.00035f * (100.0f * rel_humidity) * temperature -
              0.483199f * (100.0f * rel_humidity);
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
    float emc =
        equilibrium_moisture_content(temperature, rel_humidity) * 100.0f;
    emc = emc / 30.0;
    float eta =
        1 - 2 * (emc) + 1.5 * (std::pow(emc, 2)) - 0.5 * (std::pow(emc, 3));

    return eta * std::sqrt(1 + std::pow(wind_speed, 2)) / 0.3002;
}

}  // namespace sharp
