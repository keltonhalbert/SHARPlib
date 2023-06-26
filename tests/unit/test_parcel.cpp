#define doctest_config_implement_with_main

#include <sharplib/constants.h>
#include <sharplib/layer.h>
#include <sharplib/parcel.h>
#include <sharplib/profile.h>
#include <sharplib/winds.h>

#include <iostream>
#include <memory>
#include <string>

#include "doctest.h"

auto array_from_range = [](const float bottom, const float top,
                           const std::ptrdiff_t size) {
    auto arr = std::make_unique<float[]>(size);
    float delta = (top - bottom) / static_cast<float>(size);
    for (int k = 0; k < size; ++k) {
        arr[k] = bottom + static_cast<float>(k) * delta;
    }
    return arr;
};

auto tmpk_snd = [](const float tmpk_sfc, const float delta_tmpk_cap,
                   const float delta_tmpk_troposphere, const float hght_pbl_top,
                   const float hght_trop_top, const float height[],
                   const std::ptrdiff_t size) {
    // get the indices corresponding to the pbl top
    // and the tropopause...
    sharp::heightlayer free_troposphere =
        sharp::heightlayer(hght_pbl_top, hght_trop_top);
    sharp::layerindex idx =
        sharp::get_layer_index(free_troposphere, height, size);

    const float dse_pbl = sharp::moist_static_energy(0.0f, tmpk_sfc, 0.0f);
    // use the temperature lapse rate to get a
    // dry static energy rate of change
    const float delta_dse_trop =
        sharp::gravity - (sharp::cp_dryair * delta_tmpk_troposphere);

    auto tmpk_arr = std::make_unique<float[]>(size);
    tmpk_arr[0] = tmpk_sfc;

    // calculate the temperature of the boundary layer
    // profile assuming a constant dry static energy
    for (int k = 1; k < idx.kbot; ++k) {
        tmpk_arr[k] = (dse_pbl - height[k] * sharp::gravity) / sharp::cp_dryair;
    }

    // add a capping inversion to the top of the boundary layer
    // and use that as the starting dry static energy for the
    // free troposphere, which will then be modified by the
    // tropospheric lapse rate.
    const float dse_troposphere_bottom = sharp::moist_static_energy(
        height[idx.kbot], tmpk_arr[idx.kbot - 1] + delta_tmpk_cap, 0.0f);

    // calculate the temperature profile of the free troposphere
    for (int k = idx.kbot; k < idx.ktop + 1; ++k) {
        float dse = dse_troposphere_bottom +
                    delta_dse_trop * (height[k] - hght_pbl_top);
        tmpk_arr[k] = (dse - sharp::gravity * height[k]) / sharp::cp_dryair;
    }

    // temperature above the tropopause is isothermal
    for (int k = idx.ktop + 1; k < size; ++k) {
        tmpk_arr[k] = tmpk_arr[idx.ktop];
    }

    return tmpk_arr;
};

auto pres_dry_snd = [](const float pres_sfc, const float height[],
                       const float tmpk[], const std::ptrdiff_t size) {
    auto pres_arr = std::make_unique<float[]>(size);
    pres_arr[0] = pres_sfc;

    for (std::ptrdiff_t k = 1; k < size; ++k) {
        float tmpk_mean = 0.5f * (tmpk[k] + tmpk[k - 1]);
        float delta_z = height[k] - height[k - 1];
        float delta_p =
            -(sharp::gravity * delta_z) / (sharp::rdgas * tmpk_mean);
        pres_arr[k] = pres_arr[k - 1] * std::exp(delta_p);
    }

    return pres_arr;
};

auto pres_moist_snd = [](const float pres_sfc, const float height[],
                         const float tmpk[], const float mixr[],
                         const std::ptrdiff_t size) {
    auto pres_arr = std::make_unique<float[]>(size);
    pres_arr[0] = pres_sfc;

    for (std::ptrdiff_t k = 1; k < size; ++k) {
        const float tmpk_mean = 0.5f * (tmpk[k] + tmpk[k - 1]);
        const float mixr_mean = 0.5f * (mixr[k] + mixr[k - 1]);
        const float tv_mean = sharp::virtual_temperature(tmpk_mean, mixr_mean);
        const float delta_z = height[k] - height[k - 1];
        const float delta_p =
            -(sharp::gravity * delta_z) / (sharp::rdgas * tv_mean);
        pres_arr[k] = pres_arr[k - 1] * std::exp(delta_p);
    }

    return pres_arr;
};

auto mixr_snd = [](const float pres_sfc, const float mixr_sfc,
                   const float relh_troposphere, const float height[],
                   const float tmpk[], const std::ptrdiff_t size) {
    auto pres_arr = pres_dry_snd(pres_sfc, height, tmpk, size);
    auto mixr_arr = std::make_unique<float[]>(size);
    mixr_arr[0] = mixr_sfc;

    for (std::ptrdiff_t k = 1; k < size; ++k) {
        const float sat_mixr = sharp::mixratio(pres_arr[k], tmpk[k]);
        mixr_arr[k] = relh_troposphere * sat_mixr;
    }

    return mixr_arr;
};

auto dwpk_snd = [](const float pres[], const float mixr[],
                   const std::ptrdiff_t size) {
    auto dwpk_arr = std::make_unique<float[]>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        dwpk_arr[k] = sharp::temperature_at_mixratio(mixr[k], pres[k]);
    }
    return dwpk_arr;
};

auto vtmpk_snd = [](const float tmpk[], const float mixr[],
                    const std::ptrdiff_t size) {
    auto vtmpk_arr = std::make_unique<float[]>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        vtmpk_arr[k] = sharp::virtual_temperature(tmpk[k], mixr[k]);
    }
    return vtmpk_arr;
};

auto theta_snd = [](const float pres[], const float tmpk[],
                    const std::ptrdiff_t size) {
    auto theta_arr = std::make_unique<float[]>(size);
    for (std::ptrdiff_t k = 0; k < size; ++k) {
        theta_arr[k] = sharp::theta(pres[k], tmpk[k]);
    }
    return theta_arr;
};

test_case("testing new parcel definitions") {
    const std::ptrdiff_t n = 5000;
    const float pres_sfc = 100000.0f;
    const float tmpk_sfc = 301.5f;
    const float mixr_sfc = 0.0157f;
    const float relh_troposphere = 0.85f;
    const float hght_sfc = 0.0f;
    const float hght_top = 15000.0f;
    const float hght_tropopause = 12000.0f;
    const float delta_tmpk_cap = 1.0f;
    const float delta_tmpk_trop = 0.00725f;
    const float hght_pbl_top = 850.0f;

    // initialize the analytical sounding using
    // our input parameters
    auto hght = array_from_range(hght_sfc, hght_top, n);
    auto tmpk = tmpk_snd(tmpk_sfc, delta_tmpk_cap, delta_tmpk_trop,
                         hght_pbl_top, hght_tropopause, hght.get(), n);
    auto mixr = mixr_snd(pres_sfc, mixr_sfc, relh_troposphere, hght.get(),
                         tmpk.get(), n);
    auto pres = pres_moist_snd(pres_sfc, hght.get(), tmpk.get(), mixr.get(), n);
    auto dwpk = dwpk_snd(pres.get(), mixr.get(), n);
    auto vtmpk = vtmpk_snd(tmpk.get(), mixr.get(), n);

    // allocate an array for our buoyancy data
    auto buoy = std::make_unique<float[]>(n);

    printf("%f\n", hght[1] - hght[0]);

    sharp::parcel sfc_pcl;
    sfc_pcl.pres = pres[0];
    sfc_pcl.tmpk = tmpk[0];
    sfc_pcl.dwpk = dwpk[0];

    static constexpr sharp::lifter_wobus lifter;
    sharp::lifter_cm1 cm1_pi;
    sharp::lifter_cm1 cm1_pl;
    sharp::lifter_cm1 cm1_ai;
    sharp::lifter_cm1 cm1_al;
    cm1_pi.ma_type = sharp::adiabat::pseudo_ice;
    cm1_pl.ma_type = sharp::adiabat::pseudo_liq;
    cm1_ai.ma_type = sharp::adiabat::adiab_ice;
    cm1_al.ma_type = sharp::adiabat::adiab_liq;

    // lift and integrate the surface parcel
    sharp::lift_parcel(lifter, pres.get(), vtmpk.get(), buoy.get(), n,
                       &sfc_pcl);
    sharp::cape_cinh(pres.get(), hght.get(), buoy.get(), n, &sfc_pcl);

    std::cout << "wobus lifter" << std::endl;
    std::cout << "sfc pcl\t";
    std::cout << "lfc pres: " << sfc_pcl.lfc_pressure << "\t";
    std::cout << "el pres: " << sfc_pcl.eql_pressure << "\t";
    std::cout << "cape: " << sfc_pcl.cape << "\t";
    std::cout << "cinh: " << sfc_pcl.cinh << std::endl;

    sharp::lift_parcel(cm1_pi, pres.get(), vtmpk.get(), buoy.get(), n,
                       &sfc_pcl);
    sharp::cape_cinh(pres.get(), hght.get(), buoy.get(), n, &sfc_pcl);

    std::cout << "cm1 lifter (psuedo adiabatic with ice)" << std::endl;
    std::cout << "sfc pcl\t";
    std::cout << "lfc pres: " << sfc_pcl.lfc_pressure << "\t";
    std::cout << "el pres: " << sfc_pcl.eql_pressure << "\t";
    std::cout << "cape: " << sfc_pcl.cape << "\t";
    std::cout << "cinh: " << sfc_pcl.cinh << std::endl;

    sharp::lift_parcel(cm1_pl, pres.get(), vtmpk.get(), buoy.get(), n,
                       &sfc_pcl);
    sharp::cape_cinh(pres.get(), hght.get(), buoy.get(), n, &sfc_pcl);

    std::cout << "cm1 lifter (psuedo adiabatic with no ice)" << std::endl;
    std::cout << "sfc pcl\t";
    std::cout << "lfc pres: " << sfc_pcl.lfc_pressure << "\t";
    std::cout << "el pres: " << sfc_pcl.eql_pressure << "\t";
    std::cout << "cape: " << sfc_pcl.cape << "\t";
    std::cout << "cinh: " << sfc_pcl.cinh << std::endl;

    sharp::lift_parcel(cm1_ai, pres.get(), vtmpk.get(), buoy.get(), n,
                       &sfc_pcl);
    sharp::cape_cinh(pres.get(), hght.get(), buoy.get(), n, &sfc_pcl);

    std::cout << "cm1 lifter (adiabatic with ice)" << std::endl;
    std::cout << "sfc pcl\t";
    std::cout << "lfc pres: " << sfc_pcl.lfc_pressure << "\t";
    std::cout << "el pres: " << sfc_pcl.eql_pressure << "\t";
    std::cout << "cape: " << sfc_pcl.cape << "\t";
    std::cout << "cinh: " << sfc_pcl.cinh << std::endl;

    sharp::lift_parcel(cm1_al, pres.get(), vtmpk.get(), buoy.get(), n,
                       &sfc_pcl);
    sharp::cape_cinh(pres.get(), hght.get(), buoy.get(), n, &sfc_pcl);

    std::cout << "cm1 lifter (adiabatic with no ice)" << std::endl;
    std::cout << "sfc pcl\t";
    std::cout << "lfc pres: " << sfc_pcl.lfc_pressure << "\t";
    std::cout << "el pres: " << sfc_pcl.eql_pressure << "\t";
    std::cout << "cape: " << sfc_pcl.cape << "\t";
    std::cout << "cinh: " << sfc_pcl.cinh << std::endl;
}
