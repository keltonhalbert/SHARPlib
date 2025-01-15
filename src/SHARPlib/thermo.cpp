/**
 * \file
 * \brief Thermodynamic routines that do <!--
 * -->not directly involve parcel based ascent
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/thermo.h>

#include <algorithm>
#include <cmath>

namespace sharp {

struct lifter_wobus;
struct lifter_cm1;

float wobf(const float temperature) {
#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif
    float pol;
    const float x = temperature - ZEROCNK - 20.0f;
    if (x <= 0.0f) {
        pol = 1.0f +
              x * (-8.841660499999999e-03f +
                   x * (1.4714143e-04f +
                        x * (-9.671989000000001e-07f +
                             x * (-3.2607217e-08f + x * (-3.8598073e-10f)))));
        pol = pol * pol;
        return (15.13f / (pol * pol)) + ZEROCNK;
    } else {
        pol = x * (4.9618922e-07f +
                   x * (-6.1059365e-09f +
                        x * (3.9401551e-11f +
                             x * (-1.2588129e-13f + x * (1.6688280e-16f)))));
        pol = 1.0f + x * (3.6182989e-03f + x * (-1.3603273e-05f + pol));
        pol = pol * pol;
        return (29.93f / (pol * pol) + 0.96f * x - 14.8f) + ZEROCNK;
    }
}

float vapor_pressure(const float pressure, const float temperature) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) return MISSING;
#endif
    const float tmpc = temperature - ZEROCNK;
    const float es = 611.2f * std::exp(17.67f * tmpc / (temperature - 29.65f));
    // for extremely cold temperatures
    return std::min(es, pressure * 0.5f);
}

float vapor_pressure_ice(const float pressure, const float temperature) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) return MISSING;
#endif
    const float tmpc = temperature - ZEROCNK;
    const float es =
        611.2f * std::exp(21.8745584f * tmpc / (temperature - 7.66f));
    // for extremely cold temperatures
    return std::min(es, pressure * 0.5f);
}

float lcl_temperature(const float temperature, const float dewpoint) {
#ifndef NO_QC
    if ((temperature == MISSING) || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif
    static constexpr float c1 = 56.0f;
    static constexpr float c2 = 800.0f;

    const float term_1 = 1.0 / (dewpoint - c1);
    const float term_2 = std::log(temperature / dewpoint) / c2;
    return (1.0 / (term_1 + term_2)) + c1;
}

float temperature_at_mixratio(const float wv_mixratio, const float pressure) {
#ifndef NO_QC
    if ((wv_mixratio == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif
    float es = (wv_mixratio / EPSILON) * pressure / 100.0f /
               (1.0 + (wv_mixratio / EPSILON));
    // for extremely cold temperatures
    es = std::min(es, pressure * 0.5f);
    const float el = std::log(es);
    return ZEROCNK + (243.5f * el - 440.8f) / (19.48f - el);
}

float theta_level(const float potential_temperature, const float temperature) {
#ifndef NO_QC
    if ((potential_temperature == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif
    static constexpr float CPOR = CP_DRYAIR / RDGAS;
    return THETA_REF_PRESSURE /
           std::pow((potential_temperature / temperature), CPOR);
}

float theta(float pressure, float temperature, float ref_pressure) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING) ||
        (ref_pressure == MISSING)) {
        return MISSING;
    }
#endif
    return (temperature * std::pow(ref_pressure / pressure, ROCP));
}

float mixratio(float q) {
#ifndef NO_QC
    if (q == MISSING) return MISSING;
#endif
    return q / (1.0 - q);
}

float mixratio(float pressure, float temperature) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif

    const float e = vapor_pressure(pressure, temperature);
    return (EPSILON * e) / (pressure - e);
}

float mixratio_ice(float pressure, float temperature) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif
    const float e = vapor_pressure_ice(pressure, temperature);
    return (EPSILON * e) / (pressure - e);
}

float specific_humidity(float rv) {
#ifndef NO_QC
    if (rv == MISSING) return MISSING;
#endif
    return rv / (1.0 + rv);
}

float virtual_temperature(float temperature, float qv, float ql, float qi) {
#ifndef NO_QC
    if (qv == MISSING) {
        return temperature;
    } else if (temperature == MISSING) {
        return MISSING;
    } else if ((ql == MISSING) || (qi == MISSING)) {
        ql = 0.0f;
        qi = 0.0f;
    }
#endif
    return (temperature * ((1.0 + (qv / EPSILON)) / (1.0 + qv + ql + qi)));
}

float saturated_lift(float pressure, float theta_sat, const float converge) {
#ifndef NO_QC
    if ((pressure == MISSING) || (theta_sat == MISSING)) {
        return MISSING;
    }
#endif

    if (std::fabs(pressure - THETA_REF_PRESSURE) <= converge) return theta_sat;

    const float pwrp = std::pow(pressure / THETA_REF_PRESSURE, ROCP);
    // get the temperature
    float t1 = theta_sat * pwrp;
    float e1 = wobf(t1) - wobf(theta_sat);
    float rate = 1.0f;
    float eor = 999;
    float t2;
    int condition = 0;
    // Testing the original showed that only
    // 5 or so iterations are needed, but
    // double that just in case. It'll exit
    // early if it converges anyway.
    for (int iter = 0; iter < 10; ++iter) {
        t2 = t1 - e1 * rate;
        float e2 = (t2) / pwrp;
        e2 = e2 + wobf(t2) - wobf(e2) - theta_sat;

        eor = e2 * rate;
        rate = (t2 - t1) / (e2 - e1);
        t1 = t2;
        e1 = e2;
        condition |= (std::fabs(eor) <= converge);
        if (condition) break;
    }
    return t2 - eor;
}

float wetlift(float pressure, float temperature, float lifted_pressure) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING) ||
        (lifted_pressure == MISSING)) {
        return MISSING;
    }
#endif

    // parcels potential temperature
    const float pcl_theta = theta(pressure, temperature, THETA_REF_PRESSURE);
    // some Wobus voodo
    const float woth = wobf(pcl_theta);
    const float wott = wobf(temperature);
    // This is the wet bulb potential temperature
    const float pcl_thetaw = pcl_theta - woth + wott;
    // get the temperature that crosses the moist adiabat at
    // this pressure level
    return saturated_lift(lifted_pressure, pcl_thetaw);
}

float _solve_cm1(float& pcl_pres_next, float& pcl_pi_next, float& pcl_t_next,
                 float& pcl_rv_next, float& pcl_rl_next, float& pcl_ri_next,
                 const float pcl_pres_prev, const float pcl_t_prev,
                 const float pcl_theta_prev, const float pcl_rv_prev,
                 const float pcl_rl_prev, const float pcl_ri_prev,
                 const float rv_total, const bool ascending,
                 const bool ice = false, const float converge = 0.0002f) {
    // first guess - use the theta of the parcel
    // before lifting, and update the first guess
    // accordingly
    float pcl_theta_last = pcl_theta_prev;
    float pcl_theta_next = pcl_theta_prev;
    bool not_converged = true;

    while (not_converged) {
        pcl_t_next = pcl_theta_last * pcl_pi_next;

        float fliq = 1.0;
        float fice = 0.0;
        if (ice) {
            fliq = std::max(
                std::min((pcl_t_next - 233.15f) / (ZEROCNK - 233.15f), 1.0f),
                0.0f);
            fice = 1.0f - fliq;
        }

        const float rv_term = fliq * mixratio(pcl_pres_next, pcl_t_next) +
                              fice * mixratio_ice(pcl_pres_next, pcl_t_next);
        if (ascending) {
            pcl_rv_next = std::min(rv_total, rv_term);
            pcl_ri_next = std::max(fice * (rv_total - pcl_rv_next), 0.0f);
            pcl_rl_next = std::max(rv_total - pcl_rv_next - pcl_ri_next, 0.0f);
        } else {
            pcl_rv_next = rv_term;
            pcl_rl_next = fliq * (rv_total - pcl_rv_next);
            pcl_ri_next = fice * (rv_total - pcl_rv_next);
        }

        const float tbar = 0.5f * (pcl_t_prev + pcl_t_next);
        const float rvbar = 0.5f * (pcl_rv_prev + pcl_rv_next);
        const float rlbar = 0.5f * (pcl_rl_prev + pcl_rl_next);
        const float ribar = 0.5f * (pcl_ri_prev + pcl_ri_next);

        const float LHV = LV1 - LV2 * tbar;
        const float LHS = LS1 - LS2 * tbar;

        const float RM = RDGAS + RVGAS * rvbar;
        const float CPM =
            CP_DRYAIR + CP_VAPOR * rvbar + CP_LIQUID * rlbar + CP_ICE * ribar;
        const float term =
            LHV * (pcl_rl_next - pcl_rl_prev) / (CPM * tbar) +
            LHS * (pcl_ri_next - pcl_ri_prev) / (CPM * tbar) +
            (RM / CPM - ROCP) * std::log(pcl_pres_next / pcl_pres_prev);

        pcl_theta_next = pcl_theta_prev * std::exp(term);

        if (std::abs(pcl_theta_next - pcl_theta_last) > converge) {
            pcl_theta_last =
                pcl_theta_last + 0.3 * (pcl_theta_next - pcl_theta_last);
        } else {
            not_converged = false;
        }
    }
    return pcl_theta_next;
}

float moist_adiabat_cm1(float pressure, float temperature, float new_pressure,
                        float& rv_total, float& rv, float& rl, float& ri,
                        const float pres_incr, const float converge,
                        const adiabat ma_type) {
    // set up solver variables
    const bool ice = (ma_type >= adiabat::pseudo_ice) ? true : false;
    float dp = new_pressure - pressure;
    const bool ascending = std::signbit(dp);
    const int n_iters =
        (std::abs(dp) < pres_incr) ? 1 : 1 + (int)(std::abs(dp) / pres_incr);
    dp = (n_iters == 1) ? dp : dp / (float)n_iters;

    // Start by setting the "top" variables.
    float pcl_theta_next = theta(pressure, temperature);
    float pcl_pres_next = pressure;
    float pcl_pi_next = std::pow(pressure / THETA_REF_PRESSURE, ROCP);
    float pcl_t_next = pcl_theta_next * pcl_pi_next;
    float pcl_rv_next = rv;
    float pcl_rl_next = rl;
    float pcl_ri_next = ri;

    // Iterate the required number of times to reach the new pressure
    // level from the old one in increments of dp
    for (int i = 0; i < n_iters; ++i) {
        float pcl_theta_prev = pcl_theta_next;
        float pcl_pres_prev = pcl_pres_next;
        float pcl_rv_prev = pcl_rv_next;
        float pcl_rl_prev = pcl_rl_next;
        float pcl_ri_prev = pcl_ri_next;
        float pcl_t_prev = pcl_t_next;

        pcl_pres_next += dp;
        pcl_pi_next = std::pow(pcl_pres_next / THETA_REF_PRESSURE, ROCP);

        // call the iterative solver to get the new parcel
        // theta accounting for liquid (and ice if enabled)
        pcl_theta_next = _solve_cm1(
            pcl_pres_next, pcl_pi_next, pcl_t_next, pcl_rv_next, pcl_rl_next,
            pcl_ri_next, pcl_pres_prev, pcl_t_prev, pcl_theta_prev, pcl_rv_prev,
            pcl_rl_prev, pcl_ri_prev, rv_total, ascending, ice, converge);
        /*printf("%f %f\t%f\t%f %f\n", pcl_pres_prev, pcl_theta_prev,*/
        /*       pcl_pres_next, pcl_rv_next, rv_total);*/

        if ((ma_type == adiabat::pseudo_liq) ||
            (ma_type == adiabat::pseudo_ice)) {
            rv_total = pcl_rv_next;
            pcl_rl_next = 0.0f;
            pcl_ri_next = 0.0f;
        }
    }
    pcl_t_next = pcl_theta_next * pcl_pi_next;
    rv = pcl_rv_next;
    rl = pcl_rl_next;
    ri = pcl_ri_next;
    return pcl_t_next;
}

void drylift(float pressure, float temperature, float dewpoint,
             float& pressure_at_lcl, float& temperature_at_lcl) {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        pressure_at_lcl = MISSING;
        temperature_at_lcl = MISSING;
        return;
    }
#endif

    // theta is constant from parcel level to LCL
    const float pcl_theta = theta(pressure, temperature, THETA_REF_PRESSURE);

    temperature_at_lcl = lcl_temperature(temperature, dewpoint);
    pressure_at_lcl = theta_level(pcl_theta, temperature_at_lcl);

    if (pressure_at_lcl > pressure) pressure_at_lcl = pressure;
    return;
}

float wetbulb(lifter_wobus lifter, float pressure, float temperature,
              float dewpoint);

float wetbulb(lifter_cm1 lifter, float pressure, float temperature,
              float dewpoint);

float theta_wetbulb(lifter_wobus lifter, float pressure, float temperature,
                    float dewpoint);

float theta_wetbulb(lifter_cm1 lifter, float pressure, float temperature,
                    float dewpoint);

float thetae(float pressure, float temperature, float dewpoint) {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;
    const float mixr = mixratio(pressure, dewpoint);
    const float vappres = vapor_pressure(pressure, temperature);

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    const float theta_lcl_dry =
        theta(pressure - vappres, temperature) *
        std::pow((temperature / temperature_at_lcl), (0.28f * mixr));
    const float thetae =
        theta_lcl_dry * std::exp(mixr * (1.0 + 0.448 * mixr) *
                                 (3036.0 / temperature_at_lcl - 1.78));
    return thetae;
}

float lapse_rate(HeightLayer layer_agl, const float height[],
                 const float temperature[], const std::ptrdiff_t N) {
#ifndef NO_QC
    if ((layer_agl.bottom == MISSING) || (layer_agl.top == MISSING)) {
        return MISSING;
    }
#endif

    // convert from agl to msl
    layer_agl.bottom += height[0];
    layer_agl.top += height[0];

    // bounds check the height layer
    if (layer_agl.bottom < height[0]) {
        layer_agl.bottom = height[0];
    }
    if (layer_agl.top > height[N - 1]) {
        layer_agl.top = height[N - 1];
    }

    // lower and upper temperature
    const float tmpc_l =
        interp_height(layer_agl.bottom, height, temperature, N);
    const float tmpc_u = interp_height(layer_agl.top, height, temperature, N);
#ifndef NO_QC
    if ((tmpc_l == MISSING) || (tmpc_u == MISSING)) {
        return MISSING;
    }
#endif

    // dT/dz, positive (definition of lapse rate), in km
    const float dz = layer_agl.top - layer_agl.bottom;
    return ((tmpc_u - tmpc_l) / dz) * -1000.0f;
}

float lapse_rate(PressureLayer layer, const float pressure[],
                 const float height[], const float temperature[],
                 const std::ptrdiff_t N) {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    // bounds check the pressure layer
    if (layer.bottom > pressure[0]) {
        layer.bottom = pressure[0];
    }
    if (layer.top < pressure[N - 1]) {
        layer.top = pressure[N - 1];
    }

    HeightLayer h_layer =
        pressure_layer_to_height(layer, pressure, height, N, true);

    return lapse_rate(h_layer, height, temperature, N);
}

float lapse_rate_max(HeightLayer layer_agl, const float depth,
                     const float height[], const float temperature[],
                     const std::ptrdiff_t N, HeightLayer* max_lyr) {
    float max_lr = MISSING;
    for (float z = layer_agl.bottom; z <= (layer_agl.top - depth);
         z += layer_agl.delta) {
        HeightLayer lyr = {z, z + depth};
        float lr = lapse_rate(lyr, height, temperature, N);
        if (lr > max_lr) {
            max_lr = lr;
            if (max_lyr) {
                max_lyr->bottom = lyr.bottom;
                max_lyr->top = lyr.top;
            }
        }
    }
    return max_lr;
}

float lapse_rate_max(PressureLayer layer, const float depth,
                     const float pressure[], const float height[],
                     const float temperature[], const std::ptrdiff_t N,
                     PressureLayer* max_lyr) {
    float max_lr = MISSING;
    for (float p = layer.bottom; p >= (layer.top + depth); p += layer.delta) {
        PressureLayer lyr = {p, p - depth};
        float lr = lapse_rate(lyr, pressure, height, temperature, N);
        if (lr > max_lr) {
            max_lr = lr;
            if (max_lyr) {
                max_lyr->bottom = lyr.bottom;
                max_lyr->top = lyr.top;
            }
        }
    }
    return max_lr;
}

float buoyancy(float pcl_temperature, float env_temperature) {
    return GRAVITY * (pcl_temperature - env_temperature) / (env_temperature);
}

void buoyancy(const float pcl_temperature[], const float env_temperature[],
              float buoy_arr[], std::ptrdiff_t N) {
    for (std::ptrdiff_t k = 0; k < N; ++k) {
        buoy_arr[k] = buoyancy(pcl_temperature[k], env_temperature[k]);
    }
}

float moist_static_energy(float height_agl, float temperature,
                          float specific_humidity) {
#ifndef NO_QC
    if ((height_agl == MISSING) || (temperature == MISSING) ||
        (specific_humidity == MISSING)) {
        return MISSING;
    }
#endif

    return (CP_DRYAIR * temperature) + (EXP_LV * specific_humidity) +
           (GRAVITY * height_agl);
}

float buoyancy_dilution_potential(float temperature, float mse_bar,
                                  float saturation_mse) {
#ifndef NO_QC
    if ((temperature == MISSING) || (mse_bar == MISSING) ||
        (saturation_mse == MISSING)) {
        return MISSING;
    }
#endif
    return -1.0f * (GRAVITY / (CP_DRYAIR * temperature)) *
           (mse_bar - saturation_mse);
}

}  // end namespace sharp
