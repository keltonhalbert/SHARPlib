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

float rainfall_velocity(float rho, float rho_0, float rainwater_mixratio) {
#ifndef NO_QC
    if ((rho == MISSING) || (rho_0 == MISSING) ||
        (rainwater_mixratio == MISSING)) {
        return MISSING;
    }
#endif

    float term1 = 0.001f * rho * rainwater_mixratio;
    float term2 = rho_0 / rho;
    return 36.34 * powf(term1, 0.1364f) * powf(term2, 0.5f);
}

float rainfall_efficiency(const float pressure[], const float height[],
                          const float temperature[], const float mixr[],
                          std::ptrdiff_t N, const Parcel& pcl,
                          float rainwater_mixratio) {
    if (pcl.lcl_pressure == MISSING) return MISSING;
    PressureLayer subcloud_layer = {pressure[0], pcl.lcl_pressure};
    HeightLayer subloud_layer_hght =
        pressure_layer_to_height(subcloud_layer, pressure, height, N);

    const float rho0 = density(pressure[0], temperature[0]);
    const float rl0 = rainwater_mixratio;

    float pres_top = pcl.lcl_pressure;
    float tmpk_top =
        interp_pressure(subcloud_layer.top, pressure, temperature, N);
    float rho_top = density(subcloud_layer.top, tmpk_top);

    const float dt = 10;  // seconds
    float V = rainfall_velocity(rho_top, rho0, rainwater_mixratio);
    float hght_bot = subloud_layer_hght.top - V * dt;

    while (hght_bot > height[0]) {
        float pres_bot = interp_height(hght_bot, height, pressure, N);
        float tmpk_bot = interp_height(hght_bot, height, temperature, N);
        float mixr_bot = interp_height(hght_bot, height, mixr, N);
        float rho_bot = density(pres_bot, tmpk_bot);

        float rho_bar = (rho_top + rho_bot) / 2.0;
        float pres_bar = (pres_top + pres_bot) / 2.0;

        float evap_rate = rainfall_evaporation_rate(
            pres_bar, tmpk_bot, rho_bar, mixr_bot, rainwater_mixratio);

        rainwater_mixratio = rainwater_mixratio + (evap_rate * dt);
        V = rainfall_velocity(rho_bar, rho0, rainwater_mixratio);
        hght_bot = hght_bot - V * dt;

        pres_top = pres_bot;
        tmpk_top = tmpk_bot;
        rho_top = rho_bot;
    }

    return rainwater_mixratio / rl0;
}

}  // namespace sharp
