/**
 * \file
 * \brief Routines used to computed derived sounding parameters from<!--
 * --> vertical atmospheric profiles.
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#include "SHARPlib/params/convective.h"

#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

#include <utility>

#include "SHARPlib/interp.h"

namespace sharp {

PressureLayer effective_inflow_layer(
    sharp::lifter_wobus &lifter, const float pressure[], const float height[],
    const float temperature[], const float dewpoint[],
    const float virtemp_arr[], float pcl_vtmpk_arr[], float buoy_arr[],
    const std::ptrdiff_t N, const float cape_thresh, const float cinh_thresh,
    Parcel *mupcl);

PressureLayer effective_inflow_layer(
    sharp::lifter_cm1 &lifter, const float pressure[], const float height[],
    const float temperature[], const float dewpoint[],
    const float virtemp_arr[], float pcl_vtmpk_arr[], float buoy_arr[],
    const std::ptrdiff_t N, const float cape_thresh, const float cinh_thresh,
    Parcel *mupcl);

WindComponents storm_motion_bunkers(
    const float pressure[], const float height[], const float u_wind[],
    const float v_wind[], const std::ptrdiff_t N,
    HeightLayer mean_wind_layer_agl, HeightLayer wind_shear_layer_agl,
    const bool leftMover, const bool pressureWeighted) {
    constexpr float deviation = 7.5;  // deviation from mean wind in m/s

    PressureLayer mw_lyr = height_layer_to_pressure(mean_wind_layer_agl,
                                                    pressure, height, N, true);

    WindComponents layer_mean_wind = {MISSING, MISSING};
    layer_mean_wind =
        mean_wind(mw_lyr, pressure, u_wind, v_wind, N, pressureWeighted);

    // The shear is computed by finding the 500m deep
    // mean winds at the top and bottom of the wind_shear_layer
    // and then differencing the two. Means are not pressure weighted
    // for the non-parcel based method.
    HeightLayer h_layer_lo = {wind_shear_layer_agl.bottom,
                              wind_shear_layer_agl.bottom + 500.0f};
    HeightLayer h_layer_hi = {wind_shear_layer_agl.top - 500.0f,
                              wind_shear_layer_agl.top};

    PressureLayer p_layer_lo =
        height_layer_to_pressure(h_layer_lo, pressure, height, N, true);
    PressureLayer p_layer_hi =
        height_layer_to_pressure(h_layer_hi, pressure, height, N, true);

    WindComponents winds_lo =
        mean_wind(p_layer_lo, pressure, u_wind, v_wind, N, false);
    WindComponents winds_hi =
        mean_wind(p_layer_hi, pressure, u_wind, v_wind, N, false);

    const float shear_u = winds_hi.u - winds_lo.u;
    const float shear_v = winds_hi.v - winds_lo.v;
    const float mag = vector_magnitude(shear_u, shear_v);

    float storm_u = MISSING;
    float storm_v = MISSING;

    if (leftMover) {
        storm_u = layer_mean_wind.u - ((deviation / mag) * shear_v);
        storm_v = layer_mean_wind.v + ((deviation / mag) * shear_u);
    } else {
        storm_u = layer_mean_wind.u + ((deviation / mag) * shear_v);
        storm_v = layer_mean_wind.v - ((deviation / mag) * shear_u);
    }
    return {storm_u, storm_v};
}

[[nodiscard]] WindComponents storm_motion_bunkers(
    const float pressure[], const float height[], const float u_wind[],
    const float v_wind[], const std::ptrdiff_t N, PressureLayer eff_infl_lyr,
    const Parcel &mupcl, const bool leftMover) {
    HeightLayer shr_layer = {0, 6000.0};
    HeightLayer dflt_mw_lyr = {0.0, 6000.0};

    if (mupcl.eql_pressure == MISSING) {
        return storm_motion_bunkers(pressure, height, u_wind, v_wind, N,
                                    dflt_mw_lyr, shr_layer, leftMover, false);
    }

    const float eql_pres = mupcl.eql_pressure;
    if ((eff_infl_lyr.bottom == MISSING) || (eff_infl_lyr.top == MISSING)) {
        return storm_motion_bunkers(pressure, height, u_wind, v_wind, N,
                                    dflt_mw_lyr, shr_layer, leftMover, false);
    }

    HeightLayer eil_hght =
        pressure_layer_to_height(eff_infl_lyr, pressure, height, N, true);

    float eql_ht = interp_pressure(eql_pres, pressure, height, N);
    // get AGL
    eql_ht -= height[0];
    const float htop = 0.65 * (eql_ht - eil_hght.bottom);

    if ((htop < 3000.0f) || (eil_hght.bottom > htop)) {
        return storm_motion_bunkers(pressure, height, u_wind, v_wind, N,
                                    dflt_mw_lyr, shr_layer, leftMover, false);
    }

    HeightLayer mw_layer = {eil_hght.bottom, htop};
    return storm_motion_bunkers(pressure, height, u_wind, v_wind, N, mw_layer,
                                shr_layer, leftMover, true);
}

[[nodiscard]] std::pair<WindComponents, WindComponents> mcs_motion_corfidi(
    const float pressure[], const float height[], const float u_wind[],
    const float v_wind[], const std::ptrdiff_t N) {
    const float pres_sfc = pressure[0];

    WindComponents cloud_layer_mean;
    if (pres_sfc < 85000.0f) {
        cloud_layer_mean =
            mean_wind({pres_sfc, 30000.0}, pressure, u_wind, v_wind, N, false);
    } else {
        cloud_layer_mean =
            mean_wind({85000.0, 30000.0}, pressure, u_wind, v_wind, N, false);
    }

    HeightLayer low_layer = {0, 1500.0};  // agl
    PressureLayer low_layer_pres =
        height_layer_to_pressure(low_layer, pressure, height, N, true);

    WindComponents low_level_mean =
        mean_wind(low_layer_pres, pressure, u_wind, v_wind, N, false);

    WindComponents upshear = {cloud_layer_mean.u - low_level_mean.u,
                              cloud_layer_mean.v - low_level_mean.v};
    WindComponents downshear = {cloud_layer_mean.u + upshear.u,
                                cloud_layer_mean.v + upshear.v};

    return std::make_pair(upshear, downshear);
}

WindComponents effective_bulk_wind_difference(
    const float pressure[], const float height[], const float u_wind[],
    const float v_wind[], const std::ptrdiff_t N,
    sharp::PressureLayer effective_inflow_lyr,
    const float equilibrium_level_pressure) {
    if ((equilibrium_level_pressure == MISSING) ||
        (effective_inflow_lyr.bottom == MISSING))
        return {MISSING, MISSING};

    sharp::HeightLayer eil_hght =
        pressure_layer_to_height(effective_inflow_lyr, pressure, height, N);
    float eql_hght =
        interp_pressure(equilibrium_level_pressure, pressure, height, N);

    float depth = 0.5f * (eql_hght - eil_hght.bottom);
    sharp::HeightLayer ebwd_lyr = {eil_hght.bottom, depth};

    return sharp::wind_shear(ebwd_lyr, height, u_wind, v_wind, N);
}

float entrainment_cape(const float pressure[], const float height[],
                       const float temperature[], const float mse_arr[],
                       const float u_wind[], const float v_wind[],
                       const std::ptrdiff_t N, Parcel *pcl) {
    // if cape is zero, we get a divide by zero issue later.
    // there can "technically" be LFC/EL without CAPE because of very,
    // very shallow buoyancy near zero when searching for LFC/EL.
    // To-Do: maybe refactor LFC/EL search?
    if ((pcl->lfc_pressure == MISSING) || (pcl->lfc_pressure == MISSING) ||
        (pcl->cape == 0)) {
        return 0.0;
    }

    float *mse_diff = new float[N];
    float *mse_bar = new float[N];

    // compute MSE_bar
    mse_bar[0] = mse_arr[0];
    const float psfc = pressure[0];
    for (std::ptrdiff_t k = 1; k < N; ++k) {
        PressureLayer mn_lyr = {psfc, pressure[k]};
        mse_bar[k] = layer_mean(mn_lyr, pressure, mse_arr, N);
    }

    // compute MSE_star
    const float hsfc = height[0];
    for (std::ptrdiff_t k = 0; k < N; ++k) {
        const float tmpk = temperature[k];
        const float rsat = mixratio(pressure[k], tmpk);
        const float qsat = specific_humidity(rsat);
        const float height_agl = height[k] - hsfc;
        const float mse_star = moist_static_energy(height_agl, tmpk, qsat);
        mse_diff[k] = buoyancy_dilution_potential(tmpk, mse_bar[k], mse_star);
    }

    // compute NCAPE
    PressureLayer plyr = {pcl->lfc_pressure, pcl->eql_pressure};
    HeightLayer hlyr = pressure_layer_to_height(plyr, pressure, height, N);
    const float NCAPE = integrate_layer_trapz(hlyr, mse_diff, height, N);

    // Compute the Bunkers non-parcel based storm motion
    WindComponents strm_mtn =
        storm_motion_bunkers(pressure, height, u_wind, v_wind, N, {0.0, 6000.0},
                             {0.0, 6000.0}, false, false);

    // get the mean 0-1km storm relative wind
    HeightLayer layer = {hsfc + 0.0f, hsfc + 1000.0f};
    LayerIndex layer_idx = get_layer_index(layer, height, N);

    // loop from the surface to the last level before 1km AGL.
    float V_sr_mean = 0.0;
    int count = 0;
    for (std::ptrdiff_t k = 0; k < layer_idx.ktop + 1; ++k) {
#ifndef NO_QC
        if ((u_wind[k] == MISSING) || (v_wind[k] == MISSING)) {
            continue;
        }
#endif
        V_sr_mean +=
            vector_magnitude(u_wind[k] - strm_mtn.u, v_wind[k] - strm_mtn.v);
        count += 1;
    }
    const float u_1km = interp_height(layer.top, height, u_wind, N);
    const float v_1km = interp_height(layer.top, height, v_wind, N);

    V_sr_mean += vector_magnitude(u_1km - strm_mtn.u, v_1km - strm_mtn.v);
    V_sr_mean = V_sr_mean / (count + 1);

    // now for all of the ECAPE nonsense
    constexpr float L = 120.0;
    const float H =
        interp_pressure(pcl->eql_pressure, pressure, height, N) - hsfc;
    constexpr float sigma = 1.6;
    constexpr float alpha = 0.8;
    const float pitchfork = (VKSQ * (alpha * alpha) * (PI * PI) * L) /
                            (PRANDTL * (sigma * sigma) * H);
    const float V_sr_tilde = V_sr_mean / std::sqrt(2.0f * pcl->cape);
    const float V_sr_tilde_sq = V_sr_tilde * V_sr_tilde;
    const float N_tilde = NCAPE / pcl->cape;

    const float term1 = pitchfork / V_sr_tilde_sq;
    const float term2 = 1.0f + pitchfork + term1 * N_tilde;
    const float term3 = 4.0f * term1 * (1.0f - pitchfork * N_tilde);
    const float sqrt_term = term2 * term2 + term3;

    // in the case of a negative solution,
    // set ECAPE to 0;
    if (sqrt_term < 0) {
        return 0;
    }

    const float E_tilde = V_sr_tilde_sq + (-1.0f - pitchfork - (term1)*N_tilde +
                                           std::sqrt(sqrt_term)) /
                                              (2.0f * term1);
    delete[] mse_diff;
    delete[] mse_bar;

    return E_tilde * pcl->cape;
}

float energy_helicity_index(float cape, float helicity) {
#ifndef NO_QC
    if ((cape == MISSING) || (helicity == MISSING)) {
        return MISSING;
    }
#endif
    return (cape * helicity) / 160000.0f;
}

float supercell_composite_parameter(float mu_cape, float eff_srh,
                                    float eff_shear) {
#ifndef NO_QC
    if ((mu_cape == MISSING) || (eff_srh == MISSING) ||
        (eff_shear == MISSING)) {
        return MISSING;
    }
#endif

    if (eff_shear > 20.0) {
        eff_shear = 20.0;
    } else if (eff_shear < 10.0) {
        eff_shear = 0.0;
    }

    const float mu_cape_term = mu_cape / 1000.0f;
    const float eff_srh_term = eff_srh / 50.0f;
    const float eff_shear_term = eff_shear / 20.0f;

    return mu_cape_term * eff_srh_term * eff_shear_term;
}

float significant_tornado_parameter(Parcel pcl, float lcl_hght_agl,
                                    float storm_relative_helicity,
                                    float bulk_wind_difference) {
    float cinh_term, lcl_term, shear_term, srh_term, cape_term;
    if (std::isnan(pcl.cinh)) return 0.0;

    if (pcl.cinh >= -50.0)
        cinh_term = 1.0;
    else if (pcl.cinh < -200.0)
        cinh_term = 0.0;
    else
        cinh_term = ((200.0 + pcl.cinh) / 150.0);

    // units of comparisons are meters
    if (lcl_hght_agl < 1000.0)
        lcl_term = 1.0;
    else if (lcl_hght_agl > 2000.0)
        lcl_term = 0.0;
    else
        lcl_term = ((2000.0 - lcl_hght_agl) / 1000.0);

    // units of comparisons are m/s
    if (bulk_wind_difference > 30)
        shear_term = 1.5;
    else if (bulk_wind_difference < 12.5)
        shear_term = 0.0;
    else
        shear_term = bulk_wind_difference / 20.0;

    if (storm_relative_helicity == MISSING)
        srh_term = 0.0;
    else
        srh_term = storm_relative_helicity / 150.0;

    cape_term = pcl.cape / 1500.0;

    float stp = cape_term * cinh_term * lcl_term * srh_term * shear_term;
    if (stp < 0) stp = 0;
    return stp;
}

float significant_hail_parameter(const sharp::Parcel &mu_pcl,
                                 float lapse_rate_700_500mb, float tmpk_500mb,
                                 float freezing_lvl_agl, float shear_0_6km) {
    float mu_mixr =
        mixratio(mu_pcl.pres, mu_pcl.dwpk) * 1000.0f;  // convert to g/kg

    // restrict the inputs to specific ranges...
    shear_0_6km = std::min(shear_0_6km, 27.0f);
    shear_0_6km = std::max(shear_0_6km, 7.0f);

    mu_mixr = std::min(mu_mixr, 13.6f);
    mu_mixr = std::max(mu_mixr, 11.0f);

    tmpk_500mb = std::min(tmpk_500mb, -5.5f + ZEROCNK);

    float ship = (mu_pcl.cape * mu_mixr * lapse_rate_700_500mb *
                  ((tmpk_500mb - ZEROCNK) * -1) * shear_0_6km) /
                 42000000.0f;

    if (mu_pcl.cape < 1300.0f) {
        ship = ship * (mu_pcl.cape / 1300.0f);
    }

    if (lapse_rate_700_500mb < 5.8f) {
        ship = ship * (lapse_rate_700_500mb / 5.8f);
    }

    if (freezing_lvl_agl < 2400.0f) {
        ship = ship * (freezing_lvl_agl / 2400.0f);
    }

    return ship;
}

float precipitable_water(PressureLayer layer, const float pressure[],
                         const float mixing_ratio[], const std::ptrdiff_t N) {
    float pwat =
        integrate_layer_trapz(layer, mixing_ratio, pressure, N, 0, false) /
        (GRAVITY * RHO_LWAT);
    // convert from meters to millimeters
    return pwat * 1000.0f;
}

PressureLayer hail_growth_layer(const float pressure[],
                                const float temperature[],
                                const std::ptrdiff_t N) {
    sharp::PressureLayer hgz = {sharp::MISSING, sharp::MISSING};

    hgz.bottom =
        find_first_pressure(-10.0f + ZEROCNK, pressure, temperature, N);
    hgz.top = find_first_pressure(-30.0f + ZEROCNK, pressure, temperature, N);

    // do_stuff
    return hgz;
}

}  // end namespace sharp
