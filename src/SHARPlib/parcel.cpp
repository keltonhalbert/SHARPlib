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
#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>

#include <cstddef>

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
    if (this->lcl_pressure <= pres_arr[N - 1]) return;
    PressureLayer sat_lyr = {this->lcl_pressure, pres_arr[N - 1]};
    LayerIndex lyr_idx = get_layer_index(sat_lyr, pres_arr, N);

    float pos_buoy = 0.0;
    float pos_buoy_max = 0.0;
    float pbot = sat_lyr.bottom;
    float buoy_bot = interp_pressure(sat_lyr.bottom, pres_arr, buoy_arr, N);
    float hbot = interp_pressure(sat_lyr.bottom, pres_arr, hght_arr, N);
    // set the LFC pressure to the LCL if the buoyancy is positive
    float lfc_pres = (buoy_bot > 0) ? sat_lyr.bottom : MISSING;
    float eql_pres = MISSING;
    float lfc_pres_final = MISSING;
    float eql_pres_final = MISSING;
    bool in_pos_area = (buoy_bot > 0);

    for (std::ptrdiff_t k = lyr_idx.kbot; k < lyr_idx.ktop + 1; ++k) {
#ifndef NO_QC
        if (buoy_arr[k] == MISSING) continue;
#endif
        const float ptop = pres_arr[k];
        const float htop = hght_arr[k];
        const float buoy_top = buoy_arr[k];
        const float lyr_top = (buoy_top + buoy_bot) / 2.0f;

        if (!in_pos_area && (buoy_bot <= 0) && (buoy_top > 0)) {
            pos_buoy = 0.0;
            eql_pres = MISSING;
            in_pos_area = true;

            const float log_pbot = std::log10(pbot);
            const float log_ptop = std::log10(ptop);
            const float buoy_diff = buoy_top - buoy_bot;
            const float log_pres_diff = log_ptop - log_pbot;
            const float log_lfc_pres =
                log_pbot - buoy_bot * (log_pres_diff / buoy_diff);
            lfc_pres =
                std::min(std::pow(10.0f, log_lfc_pres), this->lcl_pressure);
        }

        const float condition = ((in_pos_area) & (lyr_top > 0));
        pos_buoy += condition * lyr_top * (htop - hbot);

        if (in_pos_area && ((buoy_bot >= 0) && (buoy_top < 0))) {
            in_pos_area = false;
            const float log_pbot = std::log10(pbot);
            const float log_ptop = std::log10(ptop);
            const float buoy_diff = buoy_top - buoy_bot;
            const float log_pres_diff = log_ptop - log_pbot;
            const float log_eql_pres =
                log_pbot - buoy_bot * (log_pres_diff / buoy_diff);
            eql_pres = std::max(std::pow(10.0f, log_eql_pres), pres_arr[N - 1]);
            if (pos_buoy > pos_buoy_max) {
                pos_buoy_max = pos_buoy;
                lfc_pres_final = lfc_pres;
                eql_pres_final = eql_pres;
            }
        }
        pbot = ptop;
        hbot = htop;
        buoy_bot = buoy_top;
    }
    if (in_pos_area) {
        eql_pres = pres_arr[N - 1];
        if (pos_buoy > pos_buoy_max) {
            pos_buoy_max = pos_buoy;
            lfc_pres_final = lfc_pres;
            eql_pres_final = eql_pres;
        }
    }
    if (pos_buoy_max > 0.0f) {
        this->lfc_pressure = lfc_pres_final;
        this->eql_pressure = eql_pres_final;
    }
}

float Parcel::maximum_parcel_level(const float pres_arr[],
                                   const float hght_arr[],
                                   const float buoy_arr[],
                                   const std::ptrdiff_t N) {
    if (this->cape == 0) return MISSING;
    if (this->eql_pressure == MISSING) return MISSING;

    float weights = 0.0;
    float integrated = 0.0f;
    sharp::PressureLayer mpl_search = {this->eql_pressure, pres_arr[N - 1]};
    const sharp::LayerIndex mpl_idx = get_layer_index(mpl_search, pres_arr, N);

    const float lyr_bottom_buoy =
        interp_pressure(mpl_search.bottom, pres_arr, buoy_arr, N);
    const float lyr_bottom_hght =
        interp_pressure(mpl_search.bottom, pres_arr, hght_arr, N);

    integrated =
        _integ_trapz(buoy_arr[mpl_idx.kbot], lyr_bottom_buoy,
                     hght_arr[mpl_idx.kbot], lyr_bottom_hght, weights, false);

    bool mpl_candidate_found = false;
    std::ptrdiff_t k = mpl_idx.kbot;
    for (; k < mpl_idx.ktop; ++k) {
        const float hght_bottom = hght_arr[k];
        const float buoy_bottom = buoy_arr[k];

        const float hght_top = hght_arr[k + 1];
        const float buoy_top = buoy_arr[k + 1];

        const float lyr_energy = _integ_trapz(buoy_top, buoy_bottom, hght_top,
                                              hght_bottom, weights, false);
        integrated += lyr_energy;
        if (std::abs(integrated) > this->cape) {
            mpl_candidate_found = true;
            integrated -= lyr_energy;
            break;
        }
    }

    if (!mpl_candidate_found) {
        this->mpl_pressure = pres_arr[mpl_idx.ktop];
        return this->mpl_pressure;
    }

    const float remaining_energy = this->cape + integrated;

    const float pres_bottom = pres_arr[k];
    const float hght_bottom = hght_arr[k];
    const float buoy_bottom = buoy_arr[k];
    const float pres_top = pres_arr[k + 1];
    const float hght_top = hght_arr[k + 1];
    const float buoy_top = buoy_arr[k + 1];

    const float hght_delta = hght_top - hght_bottom;
    if (std::abs(hght_delta) < 1e-6) {
        this->mpl_pressure = pres_arr[k];
        return pres_arr[k];
    }

    const float m = (buoy_top - buoy_bottom) / hght_delta;

    const float a = m;
    const float b = 2.0f * buoy_bottom;
    const float c = 2.0f * remaining_energy;

    const float discriminant = b * b - 4.0f * a * c;
    if (discriminant < 0) {
        this->mpl_pressure = pres_arr[k];
        return pres_arr[k];
    }
    const float delta_h = (-b + std::sqrt(discriminant)) / (2.0f * a);
    const float travel_fraction = delta_h / hght_delta;

    const float log_p_bottom = std::log(pres_bottom);
    const float log_p_top = std::log(pres_top);

    const float log_mpl_p =
        log_p_bottom + travel_fraction * (log_p_top - log_p_bottom);

    this->mpl_pressure = std::exp(log_mpl_p);
    return this->mpl_pressure;
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

DowndraftParcel::DowndraftParcel() {}

DowndraftParcel::DowndraftParcel(const float pressure, const float temperature,
                                 const float dewpoint) {
    this->pres = pressure;
    this->tmpk = temperature;
    this->dwpk = dewpoint;
}

template void DowndraftParcel::lower_parcel<lifter_wobus>(
    lifter_wobus& liftpcl, const float pressure_arr[], float pcl_tmpk_arr[],
    const std::ptrdiff_t N);

template void DowndraftParcel::lower_parcel<lifter_cm1>(
    lifter_cm1& liftpcl, const float pressure_arr[], float pcl_tmpk_arr[],
    const std::ptrdiff_t N);

void DowndraftParcel::cape_cinh(const float pres_arr[], const float hght_arr[],
                                const float buoy_arr[], const ptrdiff_t N) {
    if ((this->pres == MISSING) || (this->tmpk == MISSING) ||
        (this->dwpk == MISSING)) {
        return;
    }

    sharp::PressureLayer dcape_lyr = {pres_arr[0], this->pres};
    sharp::HeightLayer dcape_hght_lyr =
        pressure_layer_to_height(dcape_lyr, pres_arr, hght_arr, N);

    this->cinh =
        integrate_layer_trapz(dcape_hght_lyr, buoy_arr, hght_arr, N, 1);
    this->cape =
        integrate_layer_trapz(dcape_hght_lyr, buoy_arr, hght_arr, N, -1);
}

}  // end namespace sharp
