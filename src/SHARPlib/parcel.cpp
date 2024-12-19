/**
 * \file
 * \brief Routines used for parcel lifting and integration
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
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
#include <SHARPlib/thermo.h>

namespace sharp {

Parcel::Parcel() {}

Parcel::Parcel(const float pressure, const float temperature,
               const float dewpoint, const LPL lpl) {
    this->pres = pressure;
    this->tmpk = temperature;
    this->dwpk = dewpoint;
    this->source = lpl;
}

// Explicit template instatiations of template header functions.
template Parcel Parcel::most_unstable_parcel<PressureLayer, lifter_wobus>(
    PressureLayer& search_layer, lifter_wobus& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

template Parcel Parcel::most_unstable_parcel<HeightLayer, lifter_wobus>(
    HeightLayer& search_layer, lifter_wobus& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

template Parcel Parcel::most_unstable_parcel<PressureLayer, lifter_cm1>(
    PressureLayer& search_layer, lifter_cm1& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

template Parcel Parcel::most_unstable_parcel<HeightLayer, lifter_cm1>(
    HeightLayer& search_layer, lifter_cm1& lifter, const float pressure[],
    const float height[], const float temperature[], const float virtemp[],
    const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
    const std::ptrdiff_t N);

template Parcel Parcel::mixed_layer_parcel<PressureLayer>(
    PressureLayer& mix_layer, const float pressure[], const float height[],
    const float pot_temperature[], const float wv_mixratio[],
    const std::ptrdiff_t N);

template Parcel Parcel::mixed_layer_parcel<HeightLayer>(
    HeightLayer& mix_layer, const float pressure[], const float height[],
    const float pot_temperature[], const float wv_mixratio[],
    const std::ptrdiff_t N);

template void Parcel::lift_parcel<lifter_wobus>(lifter_wobus& liftpcl,
                                                const float pressure_arr[],
                                                float pcl_vtmpk_arr[],
                                                const std::ptrdiff_t N);

template void Parcel::lift_parcel<lifter_cm1>(lifter_cm1& liftpcl,
                                              const float pressure_arr[],
                                              float pcl_vtmpk_arr[],
                                              const std::ptrdiff_t N);

void Parcel::find_lfc_el(const float pres_arr[], const float hght_arr[],
                         const float buoy_arr[], const std::ptrdiff_t N) {
    PressureLayer sat_lyr = {this->lcl_pressure, pres_arr[N - 1]};
    LayerIndex lyr_idx = get_layer_index(sat_lyr, pres_arr, N);

    float lyr_bot = 0.0;
    float pos_buoy = 0.0;
    float pos_buoy_last = 0.0;
    float pbot = sat_lyr.bottom;
    float buoy_bot = interp_pressure(sat_lyr.bottom, pres_arr, buoy_arr, N);
    float hbot = interp_pressure(sat_lyr.bottom, pres_arr, hght_arr, N);
    // set the LFC pressure to the LCL if the buoyancy is positive
    float lfc_pres = (buoy_bot > 0) ? sat_lyr.bottom : MISSING;
    float eql_pres = MISSING;
    float lfc_pres_last = MISSING;
    float eql_pres_last = MISSING;

    for (std::ptrdiff_t k = lyr_idx.kbot; k < lyr_idx.ktop + 1; ++k) {
#ifndef NO_QC
        if (buoy_arr[k] == MISSING) continue;
#endif
        const float ptop = pres_arr[k];
        const float htop = hght_arr[k];
        const float buoy_top = buoy_arr[k];
        const float lyr_top = (buoy_top + buoy_bot) / 2.0f;
        // LFC condition
        if ((lyr_bot <= 0) && (lyr_top > 0)) {
            if (lfc_pres != MISSING) {
                pos_buoy_last = pos_buoy;
                lfc_pres_last = lfc_pres;
                eql_pres_last = eql_pres;
                pos_buoy = 0.0;
            }
            for (lfc_pres = pbot - 500; lfc_pres > ptop + 500;
                 lfc_pres -= 100.0) {
                const float buoy =
                    interp_pressure(lfc_pres, pres_arr, buoy_arr, N);
                if (buoy > 0) break;
            }
        }

        // keep track of buoyancy so that we pick the max CAPE layer
        const float condition = ((lfc_pres != MISSING) & (lyr_top > 0));
        pos_buoy += condition * lyr_top * (htop - hbot);
        // EL condition
        if ((lfc_pres != MISSING) && ((lyr_bot >= 0) && (lyr_top < 0))) {
            for (eql_pres = pbot - 500; eql_pres > ptop + 500;
                 eql_pres -= 100.0) {
                const float buoy =
                    interp_pressure(eql_pres, pres_arr, buoy_arr, N);
                if (buoy < 0) break;
            }
            if (pos_buoy_last > pos_buoy) {
                lfc_pres = lfc_pres_last;
                eql_pres = eql_pres_last;
                pos_buoy = pos_buoy_last;
            }
        }
        // If there is no EL, just use the last available level
        if ((k == N - 1) && (lyr_top > 0)) eql_pres = pres_arr[N - 1];
        // set loop variables
        pbot = ptop;
        hbot = htop;
        buoy_bot = buoy_top;
        lyr_bot = lyr_top;
    }
    if (pos_buoy > 0.0f) {
        this->lfc_pressure = lfc_pres;
        this->eql_pressure = eql_pres;
    }
}

void Parcel::cape_cinh(const float pres_arr[], const float hght_arr[],
                       const float buoy_arr[], const ptrdiff_t N) {
    if (this->lcl_pressure == MISSING) return;
    find_lfc_el(pres_arr, hght_arr, buoy_arr, N);
    if ((this->lfc_pressure != MISSING) && (this->eql_pressure != MISSING)) {
        PressureLayer lfc_el = {this->lfc_pressure, this->eql_pressure};
        PressureLayer lpl_lfc = {this->pres, this->lfc_pressure};
        HeightLayer lfc_el_ht =
            pressure_layer_to_height(lfc_el, pres_arr, hght_arr, N);
        HeightLayer lpl_lfc_ht =
            pressure_layer_to_height(lpl_lfc, pres_arr, hght_arr, N);

        this->cinh =
            integrate_layer_trapz(lpl_lfc_ht, buoy_arr, hght_arr, N, -1);
        this->cape = integrate_layer_trapz(lfc_el_ht, buoy_arr, hght_arr, N, 1);
    }
}

}  // end namespace sharp
