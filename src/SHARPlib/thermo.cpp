// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC.
#include <SHARPlib/constants.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#include <cmath>

namespace sharp {

float wobf(float temperature) noexcept {
#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif
    float pol;
    const float x = temperature - ZEROCNK - 20.0f;
    if (x <= 0.0f) {
       pol = 1.0f + x * (-8.841660499999999e-03f + x * ( 1.4714143e-04f \
                 + x * (-9.671989000000001e-07f + x * (-3.2607217e-08f \
                 + x * (-3.8598073e-10f)))));
       pol = pol * pol;
       return (15.13f / (pol * pol))+ZEROCNK;
    }
    else {
       pol = x * (4.9618922e-07f + x * (-6.1059365e-09f + \
             x * (3.9401551e-11f + x * (-1.2588129e-13f + \
             x * (1.6688280e-16f)))));
       pol = 1.0f + x * (3.6182989e-03f + x * (-1.3603273e-05f + pol));
       pol = pol * pol;
       return (29.93f / (pol * pol) + 0.96f * x - 14.8f)+ZEROCNK;
    }
}

float vapor_pressure(float pressure, float temperature) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) return MISSING;
#endif
    const float tmpc = temperature - ZEROCNK;
    const float es =
        611.2f * std::exp(17.67f * tmpc / (temperature - 29.65f));
	// for extremely cold temperatures
    return std::min(es, pressure * 0.5f);
}

float vapor_pressure_ice(float pressure, float temperature) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) return MISSING;
#endif
    const float tmpc = temperature - ZEROCNK;
    const float es =
        611.2f * std::exp(21.8745584f * tmpc / (temperature - 7.66f));
	// for extremely cold temperatures
    return std::min(es, pressure * 0.5f);
}

float lcl_temperature(float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif
    static constexpr float c1 = 56.0f;
    static constexpr float c2 = 800.0f;

    const float term_1 = 1.0f / (dewpoint - c1);
    const float term_2 = std::log(temperature / dewpoint) / c2;
    return (1.0f / (term_1 + term_2)) + c1;
}

float temperature_at_mixratio(float wv_mixratio, float pressure) noexcept {
#ifndef NO_QC
    if ((wv_mixratio == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif
	float es = (wv_mixratio / EPSILON) * pressure / 100.0f / (1.0f + (wv_mixratio / EPSILON));
	// for extremely cold temperatures
	es = std::min(es, pressure * 0.5f);
    const float el = std::log(es);
    return ZEROCNK + (243.5f * el - 440.8f) / (19.48f - el);
}

float theta_level(float potential_temperature, float temperature) noexcept {
#ifndef NO_QC
    if ((potential_temperature == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif
    static constexpr float CPOR = 1.0f / ROCP;
    return THETA_REF_PRESSURE /
           std::pow((potential_temperature / temperature), CPOR);
}

float theta(float pressure, float temperature, float ref_pressure) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING) ||
        (ref_pressure == MISSING)) {
        return MISSING;
    }
#endif
    return (temperature * std::pow(ref_pressure / pressure, ROCP));
}

float mixratio(float pressure, float temperature) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif

    const float e = vapor_pressure(pressure, temperature);
    return (EPSILON*e)/(pressure-e);
}

float mixratio_ice(float pressure, float temperature) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif
    const float e = vapor_pressure_ice(pressure, temperature);
    return (EPSILON*e)/(pressure-e);
}

float specific_humidity(float q) noexcept {
#ifndef NO_QC
    if (q == MISSING) return MISSING;
#endif
    return (1.0f - q) * q;
}

float virtual_temperature(float temperature, float qv, float ql,
                          float qi) noexcept {
#ifndef NO_QC
    if (qv == MISSING) {
        return temperature;
    } else if (temperature == MISSING) {
        return MISSING;
    }
    else if ((ql == MISSING) || (qi == MISSING)) {
        ql = 0.0f;
        qi = 0.0f;
    }
#endif
    return (temperature * ((1.0f + (qv / EPSILON)) / (1.0f + qv + ql + qi)));
}

float saturated_lift(float pressure, float theta_sat,
                     const float converge) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (theta_sat == MISSING)) {
        return MISSING;
    }
#endif

    if (std::fabs(pressure - THETA_REF_PRESSURE) <= converge)
        return theta_sat;

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

float wetlift(float pressure, float temperature,
              float lifted_pressure) noexcept {
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

[[nodiscard]] float moist_adiabat_cm1(float pressure, float temperature,
                                      float new_pressure, float& qv_total,
                                      float& qv, float&ql, float& qi,
                                      const float pres_incr,
                                      const float converge,
                                      const adiabat ma_type) {

    qv_total = mixratio(pressure, temperature);

    bool ice = (ma_type >= adiabat::pseudo_ice) ? true : false;
    // set up increment variables
    float dp = std::abs(pressure - new_pressure);
    int n_iters = (dp < pres_incr) ? 1 : 1 + (int)(dp / pres_incr);
    dp = (dp < pres_incr) ? dp : dp / (float)n_iters;

    // Start by setting the "top" variables.
    float pcl_theta_hi =
        theta(pressure, temperature, THETA_REF_PRESSURE);
    float pcl_pres_hi = pressure;
    float pcl_pi_hi = std::pow(pressure / THETA_REF_PRESSURE, ROCP);
    float pcl_t_hi = pcl_theta_hi * pcl_pi_hi;
    float pcl_qv_hi = mixratio(pcl_pres_hi, pcl_t_hi);
    float pcl_ql_hi = 0.0f;
    float pcl_qi_hi = 0.0f;

    for (int i = 0; i < n_iters; ++i ) {
        bool not_converged = true;

        float pcl_theta_lo = pcl_theta_hi; 
        float pcl_pres_lo = pcl_pres_hi;
        float pcl_qv_lo = pcl_qv_hi;
        float pcl_ql_lo = pcl_ql_hi;
        float pcl_qi_lo = pcl_qi_hi;
        float pcl_t_lo = pcl_t_hi; 

        pcl_pres_hi = pcl_pres_hi - dp;
        pcl_pi_hi = std::pow(pcl_pres_hi / THETA_REF_PRESSURE, ROCP);
        float pcl_theta_last = pcl_theta_lo;

        while(not_converged) {
            pcl_t_hi = pcl_theta_last * pcl_pi_hi;

            float fliq = 1.0;
            float fice = 0.0;
            if (ice) {
                fliq = std::max(
                    std::min((pcl_t_hi - 233.15f) / (ZEROCNK - 233.15f), 1.0f),
                    0.0f);
                fice = 1.0f - fliq;
            }

            float qv_term = fliq * mixratio(pcl_pres_hi, pcl_t_hi) +
                            fice * mixratio_ice(pcl_pres_hi, pcl_t_hi);
            pcl_qv_hi = std::min(qv_total, qv_term); 
            pcl_qi_hi = std::max(fice*(qv_total - pcl_qv_hi), 0.0f);
            pcl_ql_hi = std::max(qv_total - pcl_qv_hi - pcl_qi_hi, 0.0f);

            float tbar = 0.5f*(pcl_t_lo + pcl_t_hi);
            float qvbar = 0.5f*(pcl_qv_lo + pcl_qv_hi);
            float qlbar = 0.5f*(pcl_ql_lo + pcl_ql_hi);
            float qibar = 0.5f*(pcl_qi_lo + pcl_qi_hi);

            float LHV = LV1 - LV2 * tbar;
            float LHS = LS1 - LS2 * tbar;

            float RM = RDGAS + RVGAS * qvbar;
            float CPM = CP_DRYAIR + CP_VAPOR * qvbar + CP_LIQUID * qlbar +
                        CP_ICE * qibar;
            float term = LHV*(pcl_ql_hi-pcl_ql_lo)/(CPM*tbar)
                        +LHS*(pcl_qi_hi-pcl_qi_lo)/(CPM*tbar)
                        +(RM/CPM-ROCP)*std::log(pcl_pres_hi / pcl_pres_lo); 

            pcl_theta_hi = pcl_theta_lo * std::exp(term);

            if (std::abs(pcl_theta_hi - pcl_theta_last) > converge) {
                pcl_theta_last =
                    pcl_theta_last + 0.3 * (pcl_theta_hi - pcl_theta_last);
            }
            else {
                not_converged = false;
            }
        }
        if ((ma_type == adiabat::pseudo_liq) ||
            (ma_type == adiabat::pseudo_ice)) {
            qv_total = pcl_qv_hi;
            pcl_ql_hi = 0.0f;
            pcl_qi_hi = 0.0f;
        }
    }
    pcl_t_hi = pcl_theta_hi*pcl_pi_hi;
    qv = pcl_qv_hi;
    ql = pcl_ql_hi;
    qi = pcl_qi_hi;
    return pcl_t_hi; 
}


void drylift(float pressure, float temperature, float dewpoint,
             float& pressure_at_lcl, float& temperature_at_lcl) noexcept {
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

float lifted(float pressure, float temperature, float dewpoint,
             float lifted_pressure) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, lifted_pressure);
}

float wetbulb(float pressure, float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, pressure);
}

float theta_wetbulb(float pressure, float temperature,
                    float dewpoint) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, THETA_REF_PRESSURE);
}

float thetae(float pressure, float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;
    static constexpr float lift_top = 10000.0f; // 100 hPa but its Pa

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    // Lift a saturated parcel to 100 mb
    const float lifted_temperature =
        wetlift(pressure_at_lcl, temperature_at_lcl, lift_top);
    // Return the potential temperature of the 100 hPa value
    return theta(lift_top, lifted_temperature, THETA_REF_PRESSURE);
}

float lapse_rate(HeightLayer layer_agl, const float height[],
                 const float temperature[], const int N) noexcept {
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
                 const int N) noexcept {
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

float buoyancy(float pcl_temperature, float env_temperature) noexcept {
    return GRAVITY * (pcl_temperature - env_temperature) /
           (env_temperature);
}

float moist_static_energy(float height_agl, float temperature,
                          float specific_humidity) noexcept {
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

