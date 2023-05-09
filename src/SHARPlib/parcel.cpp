
/**
 * \file
 * \brief Routines used for parcel lifting and integration
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-11-09
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/profile.h>
#include <SHARPlib/thermo.h>

namespace sharp {

Parcel::Parcel() {
    // Set LPL values to MISSING
    pres = MISSING;
    tmpc = MISSING;
    dwpc = MISSING;

    // Set derived values to MISSING
    lcl_pressure = MISSING;
    lfc_pressure = MISSING;
    eql_pressure = MISSING;
    mpl_pressure = MISSING;

    cape = 0.0;
    cinh = 0.0;
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the parcel initial values to the surface of the profile
 *
 * \param prof
 * \param pcl
 *
 */
void _sfc(Profile* prof, Parcel* pcl) noexcept {
    pcl->pres = prof->pres[0];
    pcl->tmpc = prof->tmpc[0];
    pcl->dwpc = prof->dwpc[0];
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the parcel initial values to the most unstable parcel level
 *
 * \param prof
 * \param pcl
 *
 */
void _mu(Profile* prof, Parcel* pcl) noexcept {
    // Search for the most unstable parcel in the bottom
    // 400 hPa of the profile
    PressureLayer mu_layer(prof->pres[0], prof->pres[0] - 400.0);

    // layer_max returns the max, and will set the pressure
    // of the max via a pointer to a float.
    layer_max(mu_layer, prof->pres, prof->theta_e, prof->NZ, &(pcl->pres));
    pcl->tmpc = interp_pressure(pcl->pres, prof->pres, prof->tmpc, prof->NZ);
    pcl->dwpc = interp_pressure(pcl->pres, prof->pres, prof->dwpc, prof->NZ);
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the parcel initial values to the bottom 100mb mixed layer
 *
 * \param prof
 * \param pcl
 *
 */
void _ml(Profile* prof, Parcel* pcl) noexcept {
    PressureLayer mix_layer(prof->pres[0], prof->pres[0] - 100.0);

    // get the mean attributes of the lowest 100 hPa
    float mean_mixr = layer_mean(mix_layer, prof->pres, prof->mixr, prof->NZ);
    float mean_thta = layer_mean(mix_layer, prof->pres, prof->theta, prof->NZ);

    // set the parcel attributes
    pcl->pres = prof->pres[0];
    pcl->tmpc = theta(1000.0, mean_thta, prof->pres[0]);
    pcl->dwpc = temperature_at_mixratio(mean_mixr, prof->pres[0]);
}

void define_parcel(Profile* prof, Parcel* pcl, LPL source) noexcept {
    pcl->source = source;

    if (source == LPL::SFC) {
        _sfc(prof, pcl);
        return;
    } else if (source == LPL::FCST) {
        // TO-DO: Write the forecast_surface routine
        return;
    } else if (source == LPL::MU) {
        _mu(prof, pcl);
        return;
    } else if (source == LPL::ML) {
        _ml(prof, pcl);
        return;
    } else if (source == LPL::EIL) {
        // TO-DO: Write the EIL routine
        return;
    } else if (source == LPL::USR) {
        // do nothing - its already been set!
        return;
    } else {
        // TO-DO: probably should raise an error or something
        return;
    }
}

// To-Do: Make the LFC and EL search configurable by the user,
// i.e. base it off of max CAPE, lowest LFC, or highest LFC
void find_lfc_el(Parcel* pcl, float* pres_arr, float* hght_arr, float* buoy_arr,
                 int NZ) noexcept {
    PressureLayer sat_lyr = {pcl->lcl_pressure, pres_arr[NZ - 1]};
    LayerIndex lyr_idx = get_layer_index(sat_lyr, pres_arr, NZ);
    if (lyr_idx.kbot != 0) --lyr_idx.kbot;

    // check if the buoyancy at the LCL in case this is also the LFC,
    float lcl_buoy = interp_pressure(sat_lyr.bottom, pres_arr, buoy_arr, NZ);
    float lfc_pres = (lcl_buoy > 0) ? sat_lyr.bottom : MISSING;
    float eql_pres = MISSING;
    float pos_buoy = 0.0;

    float lfc_pres_last = MISSING;
    float eql_pres_last = MISSING;
    float pos_buoy_last = 0.0;

    for (int k = lyr_idx.kbot; k <= lyr_idx.ktop; ++k) {
        float pbot = pres_arr[k];
        float ptop = pres_arr[k + 1];
        float buoy_bot = buoy_arr[k];
        float buoy_top = buoy_arr[k + 1];
        float dz = hght_arr[k+1] - hght_arr[k];
        float avg_buoy = (buoy_top + buoy_bot) / 2.0;
        if ((buoy_top > 0) && (buoy_bot <= 0)) {
            if (lfc_pres != MISSING) {
                pos_buoy_last = pos_buoy;
                lfc_pres_last = lfc_pres;
                eql_pres_last = eql_pres;
                pos_buoy = 0.0;  // reset
            }
            for (lfc_pres = pbot; lfc_pres >= ptop; lfc_pres -= 1.0) {
                float buoy = interp_pressure(lfc_pres, pres_arr, buoy_arr, NZ);
                if (buoy > 0) break;
            }
        }

        if ((lfc_pres != MISSING) && (avg_buoy > 0.0)) {
            pos_buoy += dz*avg_buoy;
        }

        if ((lfc_pres != MISSING) && (buoy_top < 0) && (buoy_bot >= 0)) {
            for (eql_pres = pbot; eql_pres >= ptop; eql_pres -= 1.0) {
                float buoy = interp_pressure(eql_pres, pres_arr, buoy_arr, NZ);
                if (buoy < 0) break;
            }
            if (pos_buoy_last > pos_buoy) {
                lfc_pres = lfc_pres_last;
                eql_pres = eql_pres_last;
                pos_buoy = pos_buoy_last;
            }
        }
        if ((ptop == pres_arr[NZ - 1]) && (buoy_top > 0)) {
            eql_pres = ptop;
        }
    }
    pcl->lfc_pressure = lfc_pres;
    pcl->eql_pressure = eql_pres;
}

float cinh_below_lcl(Profile* prof, Parcel* pcl, float pres_lcl,
                     float tmpc_lcl) noexcept {
    // get the virtual temoerature of the LCL
    float vtmp_lcl = virtual_temperature(pres_lcl, tmpc_lcl, tmpc_lcl);

    // parcel thetav, or virtual potential temperature,
    // is constant from the LPL to the LCL.
    float pcl_thetav = theta(pres_lcl, vtmp_lcl, 1000.0);

    // Accumulate CINH in the mixing layer below the LCL.
    // This will be done in 10 mb increments and will use the
    // virtual temperature correction.
    float cinh = 0.0;
    for (float pbot = pcl->pres; pbot > pres_lcl; pbot -= 10.0) {
        float ptop = pbot - 10.0;
        // don't accidentally go above the LCL
        if (ptop < pres_lcl) ptop = pres_lcl;

        // _pcb - parcel bottom layer
        // _pct - parcel top layer
        // _enb - environment bottom layer
        // _ent - environment top layer

        // "lift" the thetav of the parcel to a new pressure level
        float vtmp_pcb = theta(1000.0, pcl_thetav, pbot);
        float vtmp_pct = theta(1000.0, pcl_thetav, ptop);

        // get the virtual temperature of the environment
        float vtmp_enb =
            interp_pressure(pbot, prof->pres, prof->vtmp, prof->NZ);
        float vtmp_ent =
            interp_pressure(ptop, prof->pres, prof->vtmp, prof->NZ);

        float buoy_bot = buoyancy(vtmp_pcb, vtmp_enb);
        float buoy_top = buoyancy(vtmp_pct, vtmp_ent);

        float hbot = interp_pressure(pbot, prof->pres, prof->hght, prof->NZ);
        float htop = interp_pressure(ptop, prof->pres, prof->hght, prof->NZ);

        float dz = htop - hbot;

        // integrate using trapezoid method
        float lyre = ((buoy_bot + buoy_top) / 2.0) * dz;
        if (lyre < 0.0) cinh += lyre;
    }
    return cinh;
}

void parcel_wobf(Profile* prof, Parcel* pcl) noexcept {
    constexpr lifter_wobus lifter;
    integrate_parcel(lifter, prof, pcl);
    return;
}

void parcel_wobf_new(Profile* prof, Parcel* pcl) noexcept {
    constexpr lifter_wobus lifter;
    lift_parcel(lifter, prof, pcl);

    find_lfc_el(pcl, prof->pres, prof->hght, prof->buoyancy, prof->NZ);

    if ((pcl->lfc_pressure != MISSING) && (pcl->eql_pressure != MISSING)) {
		PressureLayer lfc_el = {pcl->lfc_pressure, pcl->eql_pressure};
		PressureLayer lpl_lfc = {pcl->pres, pcl->lfc_pressure};
        HeightLayer lfc_el_ht =
            pressure_layer_to_height(lfc_el, prof->pres, prof->hght, prof->NZ);
        HeightLayer lpl_lfc_ht =
            pressure_layer_to_height(lpl_lfc, prof->pres, prof->hght, prof->NZ);
        float CINH = integrate_layer_trapz(lpl_lfc_ht, prof->buoyancy,
                                           prof->hght, prof->NZ, -1);
        float CAPE = integrate_layer_trapz(lfc_el_ht, prof->buoyancy,
                                           prof->hght, prof->NZ, 1);
		pcl->cape = CAPE;
		pcl->cinh = CINH;
    }
    return;
}

// NCAPE/buoyancy dilution from Peters et al. 2022
float buoyancy_dilution(Profile* prof, Parcel* pcl, const float* mse_star,
                        const float* mse_bar) noexcept {
    // this routine requires there to be a defined
    // LFC and EL, since these are the integration bounds
    if ((pcl->lfc_pressure == MISSING) || (pcl->eql_pressure == MISSING)) {
        return MISSING;
    }

    // set our integration layer between the LFC and EL,
    // and get the top and bottom indices of these
    // layers for iteration
    PressureLayer integ_layer = {pcl->lfc_pressure, pcl->eql_pressure};
    LayerIndex integ_idx = get_layer_index(integ_layer, prof->pres, prof->NZ);

    // start with the interpolated bottom layer
    float mse_bar_bot =
        interp_pressure(integ_layer.bottom, prof->pres, mse_bar, prof->NZ);
    float mse_star_bot =
        interp_pressure(integ_layer.bottom, prof->pres, mse_star, prof->NZ);
    float tmpk_bot =
        interp_pressure(integ_layer.bottom, prof->pres, prof->tmpc, prof->NZ) +
        ZEROCNK;
    float zbot =
        interp_pressure(integ_layer.bottom, prof->pres, prof->hght, prof->NZ);

    // initialize the top layer variables
    // to missing - these get set in the loop
    float mse_bar_top = MISSING;
    float mse_star_top = MISSING;
    float tmpk_top = MISSING;
    float ztop = MISSING;
    float NCAPE = 0.0;

    // kbot is always the first level above the pressure
    // value of the bottom layer, and ktop is always
    // the level just below the top layer
    for (int k = integ_idx.kbot; k <= integ_idx.ktop; k++) {
#ifndef NO_QC
        if ((prof->tmpc[k] == MISSING) || (prof->dwpc[k] == MISSING)) {
            continue;
        }
#endif
        // get the top-layer values
        mse_bar_top = mse_bar[k];
        mse_star_top = mse_star[k];
        tmpk_top = prof->tmpc[k] + ZEROCNK;
        ztop = prof->hght[k];

        // compute the integrated quantity
        float mse_diff_bot = -1.0 * (GRAVITY / (CP_DRYAIR * tmpk_bot)) *
                             (mse_bar_bot - mse_star_bot);
        float mse_diff_top = -1.0 * (GRAVITY / (CP_DRYAIR * tmpk_top)) *
                             (mse_bar_top - mse_star_top);
        float dz = ztop - zbot;

        NCAPE += ((mse_diff_bot + mse_diff_top) / 2.0) * dz;

        // set the top of the current layer
        // as the bottom of the next layer
        mse_bar_bot = mse_bar_top;
        mse_star_bot = mse_star_top;
        tmpk_bot = tmpk_top;
        zbot = ztop;
    }

    // finish with the interpolated top layer
    mse_bar_top =
        interp_pressure(integ_layer.top, prof->pres, mse_bar, prof->NZ);
    mse_star_top =
        interp_pressure(integ_layer.top, prof->pres, mse_star, prof->NZ);
    ztop = interp_pressure(integ_layer.top, prof->pres, prof->hght, prof->NZ);
    tmpk_top =
        interp_pressure(integ_layer.top, prof->pres, prof->tmpc, prof->NZ) +
        ZEROCNK;

    // compute the integrated quantity
    float mse_diff_bot = -1.0 * (GRAVITY / (CP_DRYAIR * tmpk_bot)) *
                         (mse_bar_bot - mse_star_bot);
    float mse_diff_top = -1.0 * (GRAVITY / (CP_DRYAIR * tmpk_top)) *
                         (mse_bar_top - mse_star_top);
    float dz = ztop - zbot;

    NCAPE += ((mse_diff_bot + mse_diff_top) / 2.0) * dz;

    return NCAPE;
}

}  // end namespace sharp

namespace sharp::exper {}  // end namespace sharp::exper
