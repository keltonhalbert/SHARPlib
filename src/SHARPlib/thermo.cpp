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
#include <iostream>

namespace sharp {

float wobf(float temperature) {
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

float vapor_pressure(float pressure, float temperature) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) return MISSING;
#endif
    const float tmpc = temperature - ZEROCNK;
    const float es = 611.2f * std::exp(17.67f * tmpc / (temperature - 29.65f));
    // for extremely cold temperatures
    return std::min(es, pressure * 0.5f);
}

float vapor_pressure_ice(float pressure, float temperature) {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) return MISSING;
#endif
    const float tmpc = temperature - ZEROCNK;
    const float es =
        611.2f * std::exp(21.8745584f * tmpc / (temperature - 7.66f));
    // for extremely cold temperatures
    return std::min(es, pressure * 0.5f);
}

float lcl_temperature(float temperature, float dewpoint) {
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

float temperature_at_mixratio(float wv_mixratio, float pressure) {
#ifndef NO_QC
    if ((wv_mixratio == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif
    float es = (wv_mixratio / EPSILON) * pressure / 100.0f /
               (1.0f + (wv_mixratio / EPSILON));
    // for extremely cold temperatures
    es = std::min(es, pressure * 0.5f);
    const float el = std::log(es);
    return ZEROCNK + (243.5f * el - 440.8f) / (19.48f - el);
}

float theta_level(float potential_temperature, float temperature) {
#ifndef NO_QC
    if ((potential_temperature == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif
    static constexpr float CPOR = 1.0f / ROCP;
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
    return q / (1.0f - q);
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
    return rv / (1.0f + rv);
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
    return (temperature * ((1.0f + (qv / EPSILON)) / (1.0f + qv + ql + qi)));
}

float density_temperature(float temperature, float qv, float qt) {
#ifndef NO_QC
    if (qv == MISSING) {
        return temperature;
    } else if (temperature == MISSING) {
        return MISSING;
    } else if ((qt == MISSING)) {
        qt = qv;
    }
#endif
    return (temperature * (1 - qt + qv / EPSILON));
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

float _solve_cm1(float& pcl_pres_hi, float& pcl_pi_hi, float& pcl_t_hi,
                 float& pcl_rv_hi, float& pcl_rl_hi, float& pcl_ri_hi,
                 float pcl_pres_lo, float pcl_t_lo, float pcl_theta_lo,
                 float pcl_rv_lo, float pcl_rl_lo, float pcl_ri_lo,
                 float rv_total, const bool ice = false,
                 const float converge = 0.0002f) {
    // first guess - use the theta of the parcel
    // before lifting, and update the first guess
    // accordingly
    float pcl_theta_last = pcl_theta_lo;
    float pcl_theta_hi = pcl_theta_lo;
    bool not_converged = true;
    int counter = 0;

    while (not_converged) {
        pcl_t_hi = pcl_theta_last * pcl_pi_hi;

        float fliq = 1.0;
        float fice = 0.0;
        if (ice) {
            fliq = std::max(
                std::min((pcl_t_hi - 233.15f) / (ZEROCNK - 233.15f), 1.0f),
                0.0f);
            fice = 1.0f - fliq;
        }

        const float rv_term = fliq * mixratio(pcl_pres_hi, pcl_t_hi) +
                              fice * mixratio_ice(pcl_pres_hi, pcl_t_hi);
        pcl_rv_hi = std::min(rv_total, rv_term);
        pcl_ri_hi = std::max(fice * (rv_total - pcl_rv_hi), 0.0f);
        pcl_rl_hi = std::max(rv_total - pcl_rv_hi - pcl_ri_hi, 0.0f);

        const float tbar = 0.5f * (pcl_t_lo + pcl_t_hi);
        const float rvbar = 0.5f * (pcl_rv_lo + pcl_rv_hi);
        const float rlbar = 0.5f * (pcl_rl_lo + pcl_rl_hi);
        const float ribar = 0.5f * (pcl_ri_lo + pcl_ri_hi);

        const float LHV = LV1 - LV2 * tbar;
        const float LHS = LS1 - LS2 * tbar;

        const float RM = RDGAS + RVGAS * rvbar;
        const float CPM =
            CP_DRYAIR + CP_VAPOR * rvbar + CP_LIQUID * rlbar + CP_ICE * ribar;
        const float term =
            LHV * (pcl_rl_hi - pcl_rl_lo) / (CPM * tbar) +
            LHS * (pcl_ri_hi - pcl_ri_lo) / (CPM * tbar) +
            (RM / CPM - ROCP) * std::log(pcl_pres_hi / pcl_pres_lo);

        pcl_theta_hi = pcl_theta_lo * std::exp(term);
        counter += 1;

        if (std::abs(pcl_theta_hi - pcl_theta_last) > converge) {
            pcl_theta_last =
                pcl_theta_last + 0.3 * (pcl_theta_hi - pcl_theta_last);
        } else {
            not_converged = false;
        }
    }
    return pcl_theta_hi;
}

float moist_adiabat_cm1(float pressure, float temperature, float new_pressure,
                        float& rv_total, float& rv, float& rl, float& ri,
                        const float pres_incr, const float converge,
                        const adiabat ma_type) {
    // This is the total water vapor mixing ratio at the beginning
    // of parcel ascent, which is most likely the LCL or LPL value
    rv_total = mixratio(pressure, temperature);

    // set up solver variables
    const bool ice = (ma_type >= adiabat::pseudo_ice) ? true : false;
    float dp = std::abs(pressure - new_pressure);
    const int n_iters = (dp < pres_incr) ? 1 : 1 + (int)(dp / pres_incr);
    dp = (dp < pres_incr) ? dp : dp / (float)n_iters;

    // Start by setting the "top" variables.
    float pcl_theta_hi = theta(pressure, temperature, THETA_REF_PRESSURE);
    float pcl_pres_hi = pressure;
    float pcl_pi_hi = std::pow(pressure / THETA_REF_PRESSURE, ROCP);
    float pcl_t_hi = pcl_theta_hi * pcl_pi_hi;
    float pcl_rv_hi = mixratio(pcl_pres_hi, pcl_t_hi);
    float pcl_rl_hi = 0.0f;
    float pcl_ri_hi = 0.0f;

    // Iterate the required number of times to reach the new pressure
    // level from the old one in increments of dp
    for (int i = 0; i < n_iters; ++i) {
        float pcl_theta_lo = pcl_theta_hi;
        float pcl_pres_lo = pcl_pres_hi;
        float pcl_rv_lo = pcl_rv_hi;
        float pcl_rl_lo = pcl_rl_hi;
        float pcl_ri_lo = pcl_ri_hi;
        float pcl_t_lo = pcl_t_hi;

        pcl_pres_hi = pcl_pres_hi - dp;
        pcl_pi_hi = std::pow(pcl_pres_hi / THETA_REF_PRESSURE, ROCP);

        // call the iterative solver to get the new parcel
        // theta accounting for liquid (and ice if enabled)
        pcl_theta_hi = _solve_cm1(pcl_pres_hi, pcl_pi_hi, pcl_t_hi, pcl_rv_hi,
                                  pcl_rl_hi, pcl_ri_hi, pcl_pres_lo, pcl_t_lo,
                                  pcl_theta_lo, pcl_rv_lo, pcl_rl_lo, pcl_ri_lo,
                                  rv_total, ice, converge);

        if ((ma_type == adiabat::pseudo_liq) ||
            (ma_type == adiabat::pseudo_ice)) {
            rv_total = pcl_rv_hi;
            pcl_rl_hi = 0.0f;
            pcl_ri_hi = 0.0f;
        }
    }
    pcl_t_hi = pcl_theta_hi * pcl_pi_hi;
    rv = pcl_rv_hi;
    rl = pcl_rl_hi;
    ri = pcl_ri_hi;
    return pcl_t_hi;
}

/**
 * \author Amelia Urquhart - OU-SoM
 * 
 * Helper method for the Peters et al 2022 saturated lapse rate.
 */
float ice_fraction(float temperature, float warmest_mixed_phase_temp,
                    float coldest_mixed_phase_temp) {
    if(temperature >= warmest_mixed_phase_temp) {
        return 0;
    } else if (temperature <= coldest_mixed_phase_temp) {
        return 1;
    } else {
        return (1 / (coldest_mixed_phase_temp - warmest_mixed_phase_temp))
            * (temperature - warmest_mixed_phase_temp);
    }
}

/**
 * \author Amelia Urquhart - OU-SoM
 * 
 * Helper method for the Peters et al 2022 saturated lapse rate.
 */
float deriv_ice_fraction(float temperature, float warmest_mixed_phase_temp,
                    float coldest_mixed_phase_temp) {
    if(temperature >= warmest_mixed_phase_temp) {
        return 0;
    } else if (temperature <= coldest_mixed_phase_temp) {
        return 0;
    } else {
        return (1 / (coldest_mixed_phase_temp - warmest_mixed_phase_temp));
    }
}

float saturation_mixing_ratio(float pressure, float temperature, int ice_flag,
                float warmest_mixed_phase_temp,
                float coldest_mixed_phase_temp) {

    float term_1, term_2, r_sat;
    
    if(ice_flag == 0) {
        term_1 = (CP_VAPOR - CP_LIQUID) / RVGAS;
        term_2 = (EXP_LV - T_TRIP * (CP_VAPOR - CP_LIQUID)) / RVGAS;

        // Saturation vapor pressure with respect to liquid I think
		float esl = std::exp((temperature - T_TRIP) * term_2 / (temperature * 
                T_TRIP)) * VAPPRES_REF * std::pow(temperature / T_TRIP, term_1);

        r_sat = EPSILON * esl / (pressure - esl);

        return r_sat;
    } else if(ice_flag == 1) {
        float omega = ice_fraction(temperature, warmest_mixed_phase_temp, coldest_mixed_phase_temp);

        float r_sat_liquid = saturation_mixing_ratio(pressure, temperature, 0, warmest_mixed_phase_temp, coldest_mixed_phase_temp);
        float r_sat_ice = saturation_mixing_ratio(pressure, temperature, 2, warmest_mixed_phase_temp, coldest_mixed_phase_temp);

        r_sat = (1 - omega) * r_sat_liquid + omega * r_sat_ice;

        return r_sat;
    } else if(ice_flag == 2) {
        term_1 = (CP_VAPOR - CP_ICE) / RVGAS;
        term_2 = (EXP_LV - T_TRIP * (CP_VAPOR - CP_ICE)) / RVGAS;

        // Saturation vapor pressure with respect to liquid I think
		float esi = std::exp((temperature - T_TRIP) * term_2 / (temperature * 
                T_TRIP)) * VAPPRES_REF * std::pow(temperature / T_TRIP, term_1);

        r_sat = EPSILON * esi / (pressure - esi);

        return r_sat;
    } else {
        return saturation_mixing_ratio(pressure, temperature, 1, warmest_mixed_phase_temp, coldest_mixed_phase_temp);
    }
}

float saturated_adiabatic_lapse_rate_peters_et_al(float temperature, 
                                        float qt, 
                                        float pressure, 
                                        float temperature_env,
                                        float qv_env,
                                        float entrainment_rate,
                                        float precip_rate,
                                        float qt_entrainment,
                                        float warmest_mixed_phase_temp,
                                        float coldest_mixed_phase_temp) {
    float omega = ice_fraction(temperature, warmest_mixed_phase_temp, 
        coldest_mixed_phase_temp);
    float d_omega = deriv_ice_fraction(temperature, warmest_mixed_phase_temp, 
        coldest_mixed_phase_temp);

    float q_vsl = (1 - qt) * 
        saturation_mixing_ratio(pressure, temperature, 0, 
            warmest_mixed_phase_temp, coldest_mixed_phase_temp);
    float q_vsi = (1 - qt) * 
        saturation_mixing_ratio(pressure, temperature, 2, 
            warmest_mixed_phase_temp, coldest_mixed_phase_temp);

    // Computes water vapor mass fraction in the parcel. Is more efficient and 
    // robust than keeping track externally
    float qv = (1 - omega) * q_vsl + omega * q_vsi;

    float temperature_entrainment = 
        -entrainment_rate * (temperature - temperature_env);
    float qv_entrainment = -entrainment_rate * (qv - qv_env);

    if(qt_entrainment == MISSING) {
        qt_entrainment = -entrainment_rate * (qt - qv_env) 
            - precip_rate * (qt - qv);
    }

    float q_condensate = qt - qv;

    float ql = q_condensate * (1 - omega);
    float qi = q_condensate * omega;

    float cp_moist_air = (1 - qt) * CP_DRYAIR + qv * CP_VAPOR + ql * CP_LIQUID
        + qi * CP_ICE;

    float density_temperature_parcel = 
        density_temperature(temperature, qv, qt);
    float density_temperature_env = 
        density_temperature(temperature_env, qv_env, qv_env);

    float buoyancy = GRAVITY 
        * (density_temperature_parcel - density_temperature_env) 
        / density_temperature_env;

    float L_v = EXP_LV + (temperature - T_TRIP) * (CP_VAPOR - CP_LIQUID);
    float L_i = EXP_LF + (temperature - T_TRIP) * (CP_LIQUID - CP_ICE);

    float L_s = L_v + omega * L_i;

    float Q_vsl = q_vsl / (EPSILON - EPSILON * qt + qv);
    float Q_vsi = q_vsi / (EPSILON - EPSILON * qt + qv);

    float Q_M = (1 - omega) * (q_vsl) / (1 - Q_vsl) + omega * (q_vsi) 
        / (1 - Q_vsi);
    float L_M = (1 - omega) * L_v * (q_vsl) 
        / (1 - Q_vsl) + omega * (L_v + L_i) * (q_vsi) / (1 - Q_vsi);
    // Moist air gas constant for environmental air
    float R_m_env = (1 - qv_env) * RDGAS + qv_env * RVGAS;

    float term_1 = buoyancy;
    float term_2 = GRAVITY;
    float term_3 = ((L_s * Q_M) / (R_m_env * temperature_env)) * GRAVITY;
    float term_4 = (cp_moist_air - L_i * (qt - qv) * d_omega) * temperature_entrainment;
    float term_5 = L_s * (qv_entrainment + qv / (1 - qt) * qt_entrainment);

    float term_6 = cp_moist_air;
    float term_7 = (L_i * (qt - qv) - L_s * (q_vsi - q_vsl)) * d_omega;
    float term_8 = (L_s * L_M) / (RVGAS * temperature * temperature);

    float dT_dz = -(term_1 + term_2 + term_3 - term_4 - term_5) 
        / (term_6 - term_7 + term_8);

    return dT_dz;
}

float moist_adiabat_peters_et_al(float pressure, float temperature,
                                      float new_pressure, float& qv, float& qt,
									  Profile* prof,
                                      const float pres_incr,
                                      float entrainment_rate,
                                      const ascent_type ma_type,
                                      const float warmest_mixed_phase_temp,
                                      const float coldest_mixed_phase_temp) {
    // Used to keep track of the index used for interpolation of the 
    // environmental profile. Removing the need to search for these each
    // iteration should save some computation time.
    int profile_index_0 = MISSING;

    float parcel_temperature = temperature;
    float parcel_qv = qv;
    float parcel_qt = qt;

    // Keeps track of changes to qt caused by entrainment and precipitation
    float dqt_dz = MISSING;

    // Ensures that the entrainment rate is zero if using non-entraining ascent
    // logic
    if(ma_type == ascent_type::adiab_nonentr
        || ma_type == ascent_type::pseudo_nonentr) {
        entrainment_rate = 0;
    }

    while(pressure > new_pressure){
        float target_pressure = pressure - pres_incr;

        if(target_pressure < new_pressure) {
            target_pressure = new_pressure;
        }

        bool need_to_find_index = false;

        // Runs on the first iteration only, sets the profile index
        if(profile_index_0 == MISSING) {
            need_to_find_index = true;
        } else {
            int profile_index_1 = std::min(profile_index_0 + 1, prof->NZ - 1);

            float pres_bottom_of_layer = prof->pres[profile_index_0];
            float pres_top_of_layer = prof->pres[profile_index_1];

            if(pres_bottom_of_layer != pres_top_of_layer 
                && pressure < pres_top_of_layer) {
                need_to_find_index = true;
            }
        }

        // If the interpolation index needs to be found, search for it
        if(need_to_find_index) {   
            // Checks if parcel is already above profile, in which case it 
            // skips the search
            if(pressure < prof->pres[prof->NZ - 1]) {
                profile_index_0 = prof->NZ - 1;
            }

            for(int i = 0; i < prof->NZ - 1; i++) {
                float pres_bottom_of_layer = prof->pres[i];
                float pres_top_of_layer = prof->pres[i + 1];

                if(pressure <= pres_bottom_of_layer 
                    && pressure > pres_top_of_layer) {
                    profile_index_0 = i;
                    break; // End search, index has been found
                }
            }
        }

        // Since the lifter does not automatically initialize qv or qt,
        // this initializes them if necessary
        if(qv == MISSING) {
            qv = saturation_mixing_ratio(pressure, temperature, 1, warmest_mixed_phase_temp, coldest_mixed_phase_temp);
            parcel_qv = qv;
        }
        if(qt == MISSING) {
            qt = saturation_mixing_ratio(pressure, temperature, 1, warmest_mixed_phase_temp, coldest_mixed_phase_temp);
            parcel_qt = qt;
        }

        // Environmental interpolation scheme used here is very simple, but
        // this can be changed later.
        int profile_index_1 = std::min(profile_index_0 + 1, prof->NZ - 1);
        float env_temperature = 
            (prof->tmpk[profile_index_0] + prof->tmpk[profile_index_1])/2;
        float qv_0 = specific_humidity(prof->mixr[profile_index_0]);
        float qv_1 = specific_humidity(prof->mixr[profile_index_1]);
        float env_qv = (qv_0 + qv_1)/2;

        float dz = RDGAS * sharp::virtual_temperature(env_temperature, env_qv) 
            / GRAVITY * std::log(pressure / target_pressure); // hypsometric

        float precip_rate; // Units: m^-1
        if(ma_type == ascent_type::pseudo_entr 
            || ma_type == ascent_type::pseudo_nonentr) {
            precip_rate = 1 / dz; // sets correct behavior for pseudoadiabatic
        } else {
            precip_rate = 0; // sets correct behavior for irrev-adiabatic
        }

        float dT_dz = saturated_adiabatic_lapse_rate_peters_et_al(
            parcel_temperature, parcel_qt, pressure, env_temperature,
            env_qv, entrainment_rate, precip_rate, dqt_dz,
            warmest_mixed_phase_temp, coldest_mixed_phase_temp);

        float new_temperature = parcel_temperature + dT_dz * dz;

        float new_parcel_qv = (1 - parcel_qt)
            * saturation_mixing_ratio(target_pressure, new_temperature, 1,
                warmest_mixed_phase_temp, coldest_mixed_phase_temp);

        float new_parcel_qt;
        if(ma_type == ascent_type::pseudo_entr 
            || ma_type == ascent_type::pseudo_nonentr) {
            dqt_dz = (new_parcel_qv - parcel_qv) / dz;
            new_parcel_qt = parcel_qv;
        } else {
            dqt_dz = -entrainment_rate * (parcel_qt - env_qv) 
                - precip_rate * (parcel_qt - parcel_qv);
            new_parcel_qt = parcel_qt + dqt_dz * dz;
        }

        parcel_temperature = new_temperature;
        parcel_qv = new_parcel_qv;
        parcel_qt = new_parcel_qt;

        pressure = target_pressure;
    }

    // Updates qv and qt in the lifter before returning the new temperature
    qv = parcel_qv;
    qt = parcel_qt;
    return parcel_temperature;
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

float lifted(float pressure, float temperature, float dewpoint,
             float lifted_pressure) {
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

float wetbulb(float pressure, float temperature, float dewpoint) {
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

float theta_wetbulb(float pressure, float temperature, float dewpoint) {
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

float thetae(float pressure, float temperature, float dewpoint) {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;
    static constexpr float lift_top = 10000.0f;  // 100 hPa but its Pa

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
                 const float temperature[], const int N) {
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
                 const float height[], const float temperature[], const int N) {
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

float buoyancy(float pcl_temperature, float env_temperature) {
    return GRAVITY * (pcl_temperature - env_temperature) / (env_temperature);
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
