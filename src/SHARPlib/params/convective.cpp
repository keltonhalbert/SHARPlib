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

#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/lifters.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

namespace sharp {

PressureLayer effective_inflow_layer(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float virtemp_arr[], float buoy_arr[],
    const int N, const float cape_thresh, const float cinh_thresh,
    Parcel *mupcl) {
    // TO-DO: At some point, this will need to be
    // templated or generalized to take other parcel
    // lifters once things progress to that level...
    static constexpr lifter_wobus lifter;

    int eff_kbot = 0;
    float eff_pbot = MISSING;
    float eff_ptop = MISSING;
    float sfc_hght = height[0];
    Parcel maxpcl;

    // search for the effective inflow bottom
    for (int k = 0; k < N; ++k) {
#ifndef NO_QC
        if ((temperature[k] == MISSING) || (dewpoint[k] == MISSING)) {
            continue;
        }
#endif
        Parcel effpcl;
        effpcl.pres = pressure[k];
        effpcl.tmpk = temperature[k];
        effpcl.dwpk = dewpoint[k];
        // We don't want to lift every single profile...
        if (height[k] - sfc_hght > 4000.0) break;

        lift_parcel(lifter, pressure, virtemp_arr, buoy_arr, N, &effpcl);
        cape_cinh(pressure, height, buoy_arr, N, &effpcl);

        if (effpcl.cape > maxpcl.cape) maxpcl = effpcl;
        if ((effpcl.cape >= cape_thresh) && (effpcl.cinh >= cinh_thresh)) {
            eff_pbot = effpcl.pres;
            eff_kbot = k;
            break;
        }
    }

    if (eff_pbot == MISSING) return {MISSING, MISSING};

    for (int k = eff_kbot + 1; k < N; ++k) {
#ifndef NO_QC
        if ((temperature[k] == MISSING) || (dewpoint[k] == MISSING)) {
            continue;
        }
#endif
        Parcel effpcl;
        effpcl.pres = pressure[k];
        effpcl.tmpk = temperature[k];
        effpcl.dwpk = dewpoint[k];
        lift_parcel(lifter, pressure, virtemp_arr, buoy_arr, N, &effpcl);
        cape_cinh(pressure, height, buoy_arr, N, &effpcl);

        if (effpcl.cape > maxpcl.cape) maxpcl = effpcl;
        if ((effpcl.cape < cape_thresh) || (effpcl.cinh < cinh_thresh)) {
            eff_ptop = effpcl.pres;
            break;
        }
    }

    if (eff_ptop == MISSING) return {MISSING, MISSING};
    *mupcl = maxpcl;
    return {eff_pbot, eff_ptop};
}

WindComponents storm_motion_bunkers(const float pressure[],
                                    const float height[], const float u_wind[],
                                    const float v_wind[], const int N,
                                    HeightLayer mean_wind_layer_agl,
                                    HeightLayer wind_shear_layer_agl,
                                    const bool leftMover = false,
                                    const bool pressureWeighted = false) {
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
    const float v_wind[], const int N, PressureLayer eff_infl_lyr,
    const Parcel *mupcl, const bool leftMover = false) {
    HeightLayer shr_layer = {0, 6000.0};
    HeightLayer dflt_mw_lyr = {0.0, 6000.0};

    if (mupcl->eql_pressure == MISSING) {
        return storm_motion_bunkers(pressure, height, u_wind, v_wind, N,
                                    dflt_mw_lyr, shr_layer, leftMover, false);
    }

    const float eql_pres = mupcl->eql_pressure;
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

float entrainment_cape(const float pressure[], const float height[],
                       const float temperature[], const float mse_arr[],
                       const float u_wind[], const float v_wind[], const int N,
                       Parcel *pcl, const float inflow_bottom = 0.0f, 
                       const float inflow_top = 1000.0f,
                       const bool left_mover = false) {
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
    for (int k = 1; k < N; ++k) {
        PressureLayer mn_lyr = {psfc, pressure[k]};
        mse_bar[k] = layer_mean(mn_lyr, pressure, mse_arr, N);
    }

    // compute MSE_star
    const float hsfc = height[0];
    for (int k = 0; k < N; ++k) {
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
                             {0.0, 6000.0}, left_mover, false);

    // get the mean inflow layer storm relative wind (default: 0-1km AGL)
    HeightLayer layer = {hsfc + inflow_bottom, hsfc + inflow_top};
    LayerIndex layer_idx = get_layer_index(layer, height, N);

    // loop from the surface to the last level before top of inflow layer.
    float V_sr_mean = 0.0;
    int count = 0;
    for (int k = 0; k < layer_idx.ktop + 1; ++k) {
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
    constexpr float sigma = 1.1; // Amelia: John ended up changing this from 
                                 // 1.6 to 1.1 between the preprint and the
                                 // final paper
    constexpr float alpha = 0.8;
    const float pitchfork = (VKSQ * (alpha * alpha) * (PI * PI) * L) /
                            (4 * PRANDTL * (sigma * sigma) * H);
                            // Amelia: I emailed John about this and he said
                            // there should be a four in the denominator
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



float entrainment_rate(const float pressure[], const float height[],
                       const float temperature[], const float mse_arr[],
                       const float u_wind[], const float v_wind[], const int N,
                       Parcel *pcl, const float inflow_bottom = 0.0f, 
                       const float inflow_top = 1000.0f,
                       const bool left_mover = false) {
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
    for (int k = 1; k < N; ++k) {
        PressureLayer mn_lyr = {psfc, pressure[k]};
        mse_bar[k] = layer_mean(mn_lyr, pressure, mse_arr, N);
    }

    // compute MSE_star
    const float hsfc = height[0];
    for (int k = 0; k < N; ++k) {
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
                             {0.0, 6000.0}, left_mover, false);

    // get the mean inflow layer storm relative wind (default: 0-1km AGL)
    HeightLayer layer = {hsfc + inflow_bottom, hsfc + inflow_top};
    LayerIndex layer_idx = get_layer_index(layer, height, N);

    // loop from the surface to the last level before 1km AGL.
    float V_sr_mean = 0.0;
    int count = 0;
    for (int k = 0; k < layer_idx.ktop + 1; ++k) {
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
    constexpr float sigma = 1.1; // Amelia: John ended up changing this from 
                                 // 1.6 to 1.1 between the preprint and the
                                 // final paper
    constexpr float alpha = 0.8;
    const float pitchfork = (VKSQ * (alpha * alpha) * (PI * PI) * L) /
                            (4 * PRANDTL * (sigma * sigma) * H);
                            // Amelia: I emailed John about this and he said
                            // there should be a four in the denominator
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

    float cape = pcl->cape;
    float ecape = E_tilde * cape;

    float E_tilde_no_vsr = E_tilde - V_sr_tilde_sq;

    float entrainment_rate = ((2 * (1 - E_tilde)) / (E_tilde + N_tilde)) / H;

    return entrainment_rate;
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
    if (pcl.cape == MISSING) return MISSING;

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

}  // end namespace sharp

namespace sharp::exper {}  // end namespace sharp::exper
