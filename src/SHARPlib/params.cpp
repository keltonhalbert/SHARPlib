/**
 * \file
 * \brief Routines used to computed derived sounding parameters from vertical atmospheric profiles.  
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */


#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>
#include <SHARPlib/profile.h>
#include <SHARPlib/parcel.h>

namespace sharp {


PressureLayer effective_inflow_layer(Profile *prof, float cape_thresh, 
                                     float cinh_thresh) noexcept {
    // TO-DO: At some point, this will need to be
    // templated or generalized to take other parcel 
    // lifters once things progress to that level...
    lifter_wobus lifter;

    // create our parcel objects
    Parcel mupcl;
    Parcel sbpcl;
    Parcel pcl;

    // find the most unstable parcel
    define_parcel(prof, &mupcl, LPL::MU);
    integrate_parcel<lifter_wobus>(lifter, prof, &mupcl);

    // find the sfc parcel to check against the MU
    // parcel, since sometimes the max Theta-E parcel
    // has less CAPE than the surface parcel
    define_parcel(prof, &sbpcl, LPL::SFC);
    integrate_parcel<lifter_wobus>(lifter, prof, &sbpcl);

    if (sbpcl.cape > mupcl.cape) {
        pcl = sbpcl;
    }
    else {
        pcl = mupcl;
    }

    float cape = pcl.cape;
    float cinh = pcl.cinh;

    if ((cape < cape_thresh) || (cinh < cinh_thresh)) {
        return {MISSING, MISSING};
    }

    int eff_kbot = 0;
    float eff_pbot = MISSING;
    float eff_ptop = MISSING;

    // search for the effective inflow bottom
    for (int k = 0; k <= prof->NZ-1; k++) {
#ifndef NO_QC
        if ((prof->tmpc[k] == MISSING) || (prof->dwpc[k] == MISSING)) {
            continue;
        }
#endif
        Parcel effpcl;
        effpcl.pres = prof->pres[k];
        effpcl.tmpc = prof->tmpc[k];
        effpcl.dwpc = prof->dwpc[k]; 
        integrate_parcel<lifter_wobus>(lifter, prof, &effpcl);

        if ((effpcl.cape >= cape_thresh) && (effpcl.cinh >= cinh_thresh)) {
            eff_pbot = effpcl.pres;
            eff_kbot = k;
            break;
        }
    }

    if (eff_pbot == MISSING) return {MISSING, MISSING};

    for (int k = eff_kbot+1; k <= prof->NZ-1; k++) {
#ifndef NO_QC
        if ((prof->tmpc[k] == MISSING) || (prof->dwpc[k] == MISSING)) {
            continue;
        }
#endif

        Parcel effpcl;
        effpcl.pres = prof->pres[k];
        effpcl.tmpc = prof->tmpc[k];
        effpcl.dwpc = prof->dwpc[k]; 
        integrate_parcel<lifter_wobus>(lifter, prof, &effpcl);

        if ((effpcl.cape < cape_thresh) || (effpcl.cinh < cinh_thresh)) {
            eff_ptop = effpcl.pres;
            break;
        }
    }

    return {eff_pbot, eff_ptop};
}

WindComponents storm_motion_bunkers(Profile *prof, 
                HeightLayer mean_wind_layer_agl, 
                HeightLayer wind_shear_layer_agl, 
                bool leftMover, bool pressureWeighted) noexcept {

    float deviation = 7.5; // deviation from mean wind in m/s

    PressureLayer mw_lyr = height_layer_to_pressure(
                mean_wind_layer_agl, prof->pres,
                prof->hght, prof->NZ, true);

    WindComponents layer_mean_wind = {MISSING, MISSING};
    if (pressureWeighted) {
        layer_mean_wind = mean_wind(mw_lyr, prof->pres, 
                prof->uwin, prof->vwin, prof->NZ); 
    }
    else {
        layer_mean_wind = mean_wind_npw(mw_lyr, prof->pres, 
                prof->uwin, prof->vwin, prof->NZ); 
    }

    // The shear is computed by finding the 500m deep
    // mean winds at the top and bottom of the wind_shear_layer
    // and then differencing the two. Means are not pressure weighted
    // for the non-parcel based method.  
    HeightLayer h_layer_lo = {wind_shear_layer_agl.zbot, 
                              wind_shear_layer_agl.zbot + 500.0f};
    HeightLayer h_layer_hi = {wind_shear_layer_agl.ztop - 500.0f,
                              wind_shear_layer_agl.ztop};

    PressureLayer p_layer_lo = height_layer_to_pressure(
                h_layer_lo, prof->pres, prof->hght,
                prof->NZ, true);
    PressureLayer p_layer_hi = height_layer_to_pressure(
                h_layer_hi, prof->pres, prof->hght,
                prof->NZ, true);

    WindComponents winds_lo = mean_wind_npw(
                p_layer_lo, prof->pres, prof->uwin, prof->vwin, prof->NZ);
    WindComponents winds_hi = mean_wind_npw(
                p_layer_hi, prof->pres, prof->uwin, prof->vwin, prof->NZ);

    float shear_u = winds_hi.u - winds_lo.u;
    float shear_v = winds_hi.v - winds_lo.v;
    float mag = vector_magnitude(shear_u, shear_v);

    float storm_u = MISSING;
    float storm_v = MISSING;

    if (leftMover) {
        storm_u = layer_mean_wind.u - ( (deviation / mag) * shear_v); 
        storm_v = layer_mean_wind.v + ( (deviation / mag) * shear_u); 
    }
    else {
        storm_u = layer_mean_wind.u + ( (deviation / mag) * shear_v); 
        storm_v = layer_mean_wind.v - ( (deviation / mag) * shear_u); 
    }

    return {storm_u, storm_v};
}


WindComponents storm_motion_bunkers(Profile* prof, bool leftMover) noexcept {

    Parcel pcl;
    Parcel sbpcl;
    Parcel mupcl;
    lifter_wobus lifter;

    define_parcel(prof, &sbpcl, LPL::SFC);
    define_parcel(prof, &mupcl, LPL::MU);

    integrate_parcel<lifter_wobus>(lifter, prof, &sbpcl);
    integrate_parcel<lifter_wobus>(lifter, prof, &mupcl);

    if (sbpcl.cape > mupcl.cape) {
        pcl = sbpcl;
    }
    else {
        pcl = mupcl;
    }

    if (pcl.eql_pressure == MISSING) {
        return storm_motion_bunkers(prof, 
                {0.0, 6000.0}, {0.0, 6000.0}, leftMover, false);
    }

    // set up the layers
    PressureLayer eil = effective_inflow_layer(prof, 100.0, -250.0);
    float eql_pres = pcl.eql_pressure;

    if ((eil.pbot == MISSING) || (eil.ptop == MISSING)) {
        return storm_motion_bunkers(prof, 
                {0.0, 6000.0}, {0.0, 6000.0}, leftMover, false);
    }

    HeightLayer eil_hght = pressure_layer_to_height(
            eil, prof->pres, prof->hght, prof->NZ, true); 

    float eql_ht = interp_pressure(eql_pres, prof->pres, prof->hght, prof->NZ);
    // get AGL
    eql_ht -= prof->hght[0];
    float htop = 0.65*(eql_ht - eil_hght.zbot);

    if ((htop < 3000.0) || (eil_hght.zbot > htop)) {
        return storm_motion_bunkers(prof, 
                {0.0, 6000.0}, {0.0, 6000.0}, leftMover, false);
    }

    HeightLayer mw_layer = {eil_hght.zbot, htop};
    HeightLayer shr_layer = {0, 6000.0};

    return storm_motion_bunkers(prof, mw_layer, shr_layer, leftMover, true);
}


float entrainment_cape(Profile* prof, Parcel *pcl) noexcept {

    float *mse_star = new float[prof->NZ];
    float *mse_bar = new float[prof->NZ];

    // compute MSE_bar
    mse_bar[0] = prof->moist_static_energy[0];
    for (int k = 1; k < prof->NZ; k++) {
        float lyr_bot = prof->moist_static_energy[0];
        float lyr_top = MISSING;
        float pbot = prof->pres[0];
        float ptop = MISSING;
        float dp = 0.0;

        float lyr_mean = 0.0;
        float weight = 0.0;
        for (int iz = 1; iz <= k; iz++) {
            if (prof->moist_static_energy[iz] == MISSING) {
                continue;
            }

            lyr_top = prof->moist_static_energy[iz];
            ptop = prof->pres[iz];
            dp = pbot - ptop;
            
            lyr_mean += ((lyr_top + lyr_bot) / 2.0) * dp;
            weight += dp;

            lyr_bot = lyr_top;
            pbot = ptop;
        }

        mse_bar[k] = lyr_mean / weight;
    }

    // compute MSE_star
    for (int k = 0; k < prof->NZ; k++) {
        // default units are g/kg - convert to kg/kg
        float rsat = mixratio(prof->pres[k], prof->tmpc[k]) / 1000.0f;
        float qsat = (1.0 - rsat ) * rsat;
        float height_agl = prof->hght[k] - prof->hght[0];
        //printf("%d %f %f %f %f %f %f\n", k, prof->pres[k], prof->hght[k], prof->tmpc[k], prof->dwpc[k], rsat, qsat);
        mse_star[k] = moist_static_energy(height_agl, prof->tmpc[k] + ZEROCNK, qsat); 
    }

    // compute NCAPE
    float NCAPE = buoyancy_dilution(prof, pcl,mse_star, mse_bar); 
    //printf("NCAPE = %f\n", NCAPE);

    // Compute the Bunkers non-parcel based storm motion
    WindComponents strm_mtn =  storm_motion_bunkers(prof, 
            {0.0, 6000.0}, {0.0, 6000.0}, false, false);
    //printf("Bunkers U = %f, Bunkers V = %f\n", strm_mtn.u, strm_mtn.v);

    // get the mean 0-1km storm relative wind
    HeightLayer layer = {prof->hght[0] + 0.0f, prof->hght[0] + 1000.0f};
    LayerIndex layer_idx = get_layer_index(layer, prof->hght, prof->NZ); 

    // loop from the surface to the last level before 1km AGL. 
    float V_sr_mean = 0.0;
    int count = 0;
    for (int k = 0; k <= layer_idx.ktop; ++k) {
#ifndef NO_QC
        if ((prof->uwin[k] == MISSING) || (prof->vwin[k] == MISSING)) {
            continue;
        }
#endif
        V_sr_mean += vector_magnitude(prof->uwin[k] - strm_mtn.u, 
                                      prof->vwin[k] - strm_mtn.v); 
        count += 1;
    }
    float u_1km = interp_height(layer.ztop, prof->hght, prof->uwin, prof->NZ);
    float v_1km = interp_height(layer.ztop, prof->hght, prof->vwin, prof->NZ);

    V_sr_mean += vector_magnitude(u_1km - strm_mtn.u, v_1km - strm_mtn.v);
    // the +2 is +1 for swapping from zero-index to size, and then
    // the last level of the interpolated top. 
    V_sr_mean = V_sr_mean / ( count + 1);
    //printf("V_sr_mean = %f\n", V_sr_mean);


    // now for all of the ECAPE nonsense
    float L = 120.0;
    float H = interp_pressure(pcl->eql_pressure, prof->pres, prof->hght, prof->NZ) - prof->hght[0];
    float sigma = 1.6;
    float alpha = 0.8;
    float pitchfork = (VKSQ * (alpha*alpha) * (PI*PI) * L) / 
                      (PRANDTL * (sigma*sigma) * H);
    float V_sr_tilde = V_sr_mean / std::sqrt(2.0f*pcl->cape);
    float V_sr_tilde_sq = V_sr_tilde*V_sr_tilde;
    float N_tilde = NCAPE / pcl->cape;

    float term1 = pitchfork / V_sr_tilde_sq;
    float term2 = 1.0f + pitchfork + term1*N_tilde;
    float term3 = 4.0f * term1 * (1.0f - pitchfork * N_tilde);
    float sqrt_term = term2*term2 + term3;

    //printf("H = %f\tpitchfork = %f\n", H, pitchfork);
    //printf("V_sr_tilde = %f\tN_tilde = %f\n", V_sr_tilde, N_tilde);

    // in the case of a negative solution, 
    // set ECAPE to 0;
    if (sqrt_term < 0) {
        return 0;
    }

    float E_tilde = V_sr_tilde_sq + (-1.0f - pitchfork - (term1)*N_tilde + 
                    std::sqrt(sqrt_term) ) / ( 2*term1 );

    delete[] mse_star;
    delete[] mse_bar;

    return E_tilde * pcl->cape;
}

float energy_helicity_index(float cape, float helicity) noexcept {
#ifndef NO_QC
    if ((cape == MISSING) || (helicity == MISSING)) {
        return MISSING;
    }
#endif
    return (cape * helicity) / 160000.0;
}


float supercell_composite_parameter(float mu_cape, float eff_srh,
                                    float eff_shear) noexcept {
#ifndef NO_QC
    if ((mu_cape == MISSING) || (eff_srh == MISSING) || (eff_shear == MISSING)) {
        return MISSING;
    }
#endif

    if (eff_shear > 20.0) {
        eff_shear = 20.0;
    }
    else if (eff_shear < 10.0) {
        eff_shear = 0.0;
    }

    float mu_cape_term = mu_cape / 1000.0;
    float eff_srh_term = eff_srh / 50.0;
    float eff_shear_term = eff_shear / 20.0;

    return mu_cape_term * eff_srh_term * eff_shear_term;
}

float significant_tornado_parameter(Profile* prof, Parcel pcl,
                                    float storm_relative_helicity,
                                    float bulk_wind_difference) noexcept {
    float cinh_term, lcl_term, shear_term, srh_term, cape_term;
    if (pcl.cape == MISSING) return MISSING;

    if (pcl.cinh >= -50.0) cinh_term = 1.0;
    else if (pcl.cinh < -200.0) cinh_term = 0.0;
    else cinh_term = ((200.0 + pcl.cinh)/150.0);


    float lcl_hght = interp_pressure(pcl.lcl_pressure, prof->pres, 
                                            prof->hght, prof->NZ); 
    lcl_hght -= prof->hght[0]; // convert to AGL

    // units of comparisons are meters
    if (lcl_hght < 1000.0) lcl_term = 1.0;
    else if (lcl_hght > 2000.0) lcl_term = 0.0;
    else lcl_term = ((2000.0 - lcl_hght) / 1000.0 );

    // units of comparisons are m/s
    if (bulk_wind_difference > 30.8667) shear_term = 1.5; 
    else if (bulk_wind_difference < 12.8611) shear_term = 0.0;
    else shear_term = bulk_wind_difference / 20.5778;

    if (storm_relative_helicity == MISSING) srh_term = 0.0;
    else srh_term = storm_relative_helicity / 150.0;

    cape_term = pcl.cape / 1500.0;

    float stp = cape_term * cinh_term * lcl_term * srh_term * shear_term;
    if (stp < 0) stp = 0;
    return stp;
}


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


