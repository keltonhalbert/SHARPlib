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
#include <SHARPlib/utils.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>
#include <SHARPlib/profile.h>
#include <SHARPlib/parcel.h>

namespace sharp {


PressureLayer effective_inflow_layer(Profile *prof, float cape_thresh, 
                                     float cinh_thresh) {
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


float entrainment_cape(Profile* prof, Parcel *pcl) {

    float *mse_star = new float[prof->NZ];
    float *mse_bar = new float[prof->NZ];

    // compute MSE_bar
    mse_bar[0] = prof->moist_static_energy[0];
    for (int k = 1; k < prof-> NZ; k++) {
        float mean_mse = 0.0;
        int missing_counter = 0;
        for (int iz = 0; iz <= k; iz ++) {
            if (prof->moist_static_energy[iz] == MISSING) {
                missing_counter += 1;
                continue;
            }
            else {
                mean_mse += prof->moist_static_energy[iz];
            }
        }
        mse_bar[k] = mean_mse / (k + 1 - missing_counter);
        //printf("mse_bar[%d] = %f\n", k, mse_bar[k]);
    }

    // compute MSE_star
    for (int k = 0; k < prof->NZ; k++) {
        // default units are g/kg - convert to kg/kg
        float rsat = mixratio(prof->pres[k], prof->tmpc[k]) / 1000.0f;
        float qsat = (1.0 - rsat ) * rsat;
        //printf("%d %f %f %f %f %f %f\n", k, prof->pres[k], prof->hght[k], prof->tmpc[k], prof->dwpc[k], rsat, qsat);
        mse_star[k] = moist_static_energy(prof->hght[k] - prof->hght[0], 
                                    prof->tmpc[k] + ZEROCNK, qsat); 
        //printf("mse_star[%d] = %f\n", k, mse_star[k]);
        //printf("mse[%d] = %f\n", k, prof->moist_static_energy[k]);
    }

    // compute NCAPE
    float NCAPE = buoyancy_dilution_ncape(prof, pcl,mse_star, mse_bar); 
    //printf("NCAPE = %f\n", NCAPE);

    // Compute the Bunkers non-parcel based storm motion
    WindComponents strm_mtn = storm_motion_bunkers_np(prof->pres, prof->hght,
                                  prof->uwin, prof->vwin, prof->NZ, 0);
    //printf("Bunkers U = %f, Bunkers V = %f\n", strm_mtn.u, strm_mtn.v);

    // get the mean 0-1km storm relative wind
    HeightLayer layer = {prof->hght[0] + 0.0f, prof->hght[0] + 1000.0f};
    LayerIndex layer_idx = get_layer_index(layer, prof->hght, prof->NZ); 

    // loop from the surface to the last level before 1km AGL. 
    float V_sr_mean = 0.0;
    for (int k = 0; k <= layer_idx.ktop; ++k) {
#ifndef NO_QC
        if ((prof->uwin[k] == MISSING) || (prof->vwin[k] == MISSING)) {
            continue;
        }
#endif
        V_sr_mean += vector_magnitude(prof->uwin[k] - strm_mtn.u, 
                                      prof->vwin[k] - strm_mtn.v); 
    }
    float u_1km = interp_height(layer.ztop, prof->hght, prof->uwin, prof->NZ);
    float v_1km = interp_height(layer.ztop, prof->hght, prof->vwin, prof->NZ);

    V_sr_mean += vector_magnitude(u_1km - strm_mtn.u, v_1km - strm_mtn.v);
    // the +2 is +1 for swapping from zero-index to size, and then
    // the last level of the interpolated top. 
    V_sr_mean = V_sr_mean / (layer_idx.ktop + 2);
    //printf("V_sr_mean = %f\n", V_sr_mean);


    // now for all of the ECAPE nonsense
    float L = 120.0;
    float H = interp_pressure(pcl->eql_pressure, prof->pres, prof->hght, prof->NZ);
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

float energy_helicity_index(float cape, float helicity) {
#ifndef NO_QC
    if ((cape == MISSING) || (helicity == MISSING)) {
        return MISSING;
    }
#endif
    return (cape * helicity) / 160000.0;
}


float supercell_composite_parameter(float mu_cape, float eff_srh,
                                                   float eff_shear) {
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


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


