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

float pyrocumulonimbus_firepower_threshold(
    PressureLayer mix_layer, const float pressure[], const float height[],
    const float temperature[], const float mixratio[], const float virtemp[],
    const float uwin[], const float vwin[], const float potential_temperature[],
    float pcl_vtmpk_arr[], float pcl_buoy_arr[], std::ptrdiff_t N, Parcel* pcl,
    float beta_incr, float phi) {
    lifter_cm1 lifter;
    lifter.ma_type = adiabat::pseudo_liq;

    float mean_theta =
        layer_mean(mix_layer, pressure, potential_temperature, N);
    float mean_mixr = layer_mean(mix_layer, pressure, mixratio, N);
    WindComponents mean_uv =
        mean_wind(mix_layer, pressure, uwin, vwin, N, false);
    float mean_wspd = sharp::vector_magnitude(mean_uv.u, mean_uv.v);

    float z_fc = sharp::MISSING;
    float delta_theta = 0.0;

    float beta_max = 0.1;
    float beta_val = 0.0;
    Parcel fire_pcl;
    while (beta_val <= beta_max) {
        float eql_tmpk = 999.0;
        float plume_theta =
            pft_plume_potential_temperature(mean_theta, beta_val);
        float plume_mixr =
            pft_plume_mixratio(mean_theta, mean_mixr, beta_val, phi);

        float plume_tmpk = theta(THETA_REF_PRESSURE, plume_theta, pressure[0]);
        float plume_dwpk = temperature_at_mixratio(plume_mixr, pressure[0]);

        fire_pcl = Parcel(pressure[0], plume_tmpk, plume_dwpk, LPL::USR);
        fire_pcl.lift_parcel(lifter, pressure, pcl_vtmpk_arr, N);
        buoyancy(pcl_vtmpk_arr, virtemp, pcl_buoy_arr, N);
        fire_pcl.cape_cinh(pressure, height, pcl_buoy_arr, N);

        if (fire_pcl.eql_pressure != MISSING) {
            eql_tmpk = interp_pressure(fire_pcl.eql_pressure, pressure,
                                       temperature, N);
        }

        if ((fire_pcl.cinh >= -sharp::TOL) && (fire_pcl.cape > 0) &&
            (eql_tmpk <= 253.15)) {
            z_fc = interp_pressure(fire_pcl.lfc_pressure, pressure, height, N);
            float theta_fc = interp_pressure(fire_pcl.lfc_pressure, pressure,
                                             potential_temperature, N);
            delta_theta = theta_fc - mean_theta;
            break;
        }

        beta_val += beta_incr;
    }

    if (z_fc == MISSING) return MISSING;
    z_fc = z_fc - height[0];

    constexpr float beta_prime = 0.4;
    constexpr float alpha_prime = 0.32;
    const float big_const =
        (PI * CP_DRYAIR) *
        std::pow((beta_prime / (1 + (alpha_prime * beta_prime))), 2.0f);

    float pres_pl_c = ((fire_pcl.lfc_pressure - fire_pcl.pres) /
                       (1 + (alpha_prime * beta_prime))) +
                      pressure[0];
    float theta_pl_c =
        sharp::interp_pressure(pres_pl_c, pressure, potential_temperature, N);

    float rho = (pres_pl_c / (sharp::RDGAS * theta_pl_c)) *
                std::pow(sharp::THETA_REF_PRESSURE / pres_pl_c, sharp::ROCP);

    float PFT =
        big_const * rho * std::pow(z_fc, 2.0f) * mean_wspd * delta_theta;

    if (pcl) *pcl = fire_pcl;

    return PFT;
}

}  // namespace sharp
