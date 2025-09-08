/**
 * \file
 * \brief Routines used to computed derived sounding parameters <!--
 * --> from vertical atmospheric profiles.
 * \author
 *     Kelton Halbert                  \n
 *     Email: kelton.halbert@noaa.gov  \n
 * \date 2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef SHARP_PARAMS_H
#define SHARP_PARAMS_H

#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes the Effective Inflow Layer for a given atmospheric sounding.
 *
 * Computes the Effective Inflow Layer, or the layer of the atmosphere believed
 * to be the primary source of inflow for supercell thunderstorms. The Effective
 * Inflow Layer, and its use in computing shear and storm relative helicity, is
 * described by Thompson et al. 2007:
 *     https://www.spc.noaa.gov/publications/thompson/effective.pdf
 *
 * Standard/default values for cape_thresh and cinh_thresh have been
 * experimentally determined to be cape_thresh = 100 J/kg and
 * cinh_thresh = -250 J/kg. If a pointer to a parcel object is passed
 * as the last arguments, the most ustanble parcel found during the
 * EIL search will be returned to that pointer.
 *
 * All arrays are treated as inputs (and expected to be precomputed) except
 * for the buoyancy array. It is more accurately meant to be an empty buffer
 * array that can be used to hold buoyancy values during parcel lifting and
 * integration.
 *
 * \param   lifter          (parcel lifting functor)
 * \param   pressure        (Pa)
 * \param   height          (meters)
 * \param   temperature     (K)
 * \param   dewpoint        (K)
 * \param   virtemp_arr     (K)
 * \param   pcl_vtmpk_arr   (K)
 * \param   buoy_arr        (m/s^2)
 * \param   N               (length of arrays)
 * \param   cape_thresh     (J/kg; default=100.0)
 * \param   cinh_thresh     (J/kg; default=-250.0)
 * \param   mupcl           (store the most unstable parcel; defualt=nullptr)
 *
 * \return  {bottom, top}
 */
template <typename Lifter>
[[nodiscard]] PressureLayer effective_inflow_layer(
    Lifter& lifter, const float pressure[], const float height[],
    const float temperature[], const float dewpoint[],
    const float virtemp_arr[], float pcl_vtmpk_arr[], float buoy_arr[],
    const std::ptrdiff_t N, const float cape_thresh = 100.0,
    const float cinh_thresh = -250.0, Parcel* mupcl = nullptr) {
    int eff_kbot = 0;
    float eff_pbot = MISSING;
    float eff_ptop = MISSING;
    float sfc_hght = height[0];
    Parcel maxpcl;

    // search for the effective inflow bottom
    for (std::ptrdiff_t k = 0; k < N; ++k) {
#ifndef NO_QC
        if ((temperature[k] == MISSING) || (dewpoint[k] == MISSING)) {
            continue;
        }
#endif
        Parcel effpcl(pressure[k], temperature[k], dewpoint[k], LPL::MU);
        // We don't want to lift every single profile...
        if (height[k] - sfc_hght > 4000.0) break;

        effpcl.lift_parcel(lifter, pressure, pcl_vtmpk_arr, N);
        buoyancy(pcl_vtmpk_arr, virtemp_arr, buoy_arr, N);
        effpcl.cape_cinh(pressure, height, buoy_arr, N);

        if (effpcl.cape > maxpcl.cape) maxpcl = effpcl;
        if ((effpcl.cape >= cape_thresh) && (effpcl.cinh >= cinh_thresh)) {
            eff_pbot = effpcl.pres;
            eff_kbot = k;
            break;
        }
    }

    if (eff_pbot == MISSING) {
        if (mupcl) *mupcl = maxpcl;
        return {MISSING, MISSING};
    }

    for (std::ptrdiff_t k = eff_kbot + 1; k < N; ++k) {
#ifndef NO_QC
        if ((temperature[k] == MISSING) || (dewpoint[k] == MISSING)) {
            continue;
        }
#endif
        Parcel effpcl(pressure[k], temperature[k], dewpoint[k], LPL::MU);
        effpcl.lift_parcel(lifter, pressure, pcl_vtmpk_arr, N);
        buoyancy(pcl_vtmpk_arr, virtemp_arr, buoy_arr, N);
        effpcl.cape_cinh(pressure, height, buoy_arr, N);

        if (effpcl.cape > maxpcl.cape) maxpcl = effpcl;
        if ((effpcl.cape < cape_thresh) || (effpcl.cinh < cinh_thresh)) {
            eff_ptop = effpcl.pres;
            break;
        }
    }

    if (mupcl) *mupcl = maxpcl;
    if (eff_ptop == MISSING) return {MISSING, MISSING};
    return {eff_pbot, eff_ptop};
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Estimates supercell storm motion using the Bunkers et al. 2000 method.
 *
 * Estimates supercell storm motion using the Bunkers et al. 2000 method
 * described in the following paper:
 *     https://doi.org/10.1175/1520-0434(2000)015%3C0061:PSMUAN%3E2.0.CO;2
 *
 * This does not use any of the updated methods described by Bunkers et al.
 * 2014, which uses Effective Inflow Layer metricks to get better estimates of
 * storm motion, especially when considering elevated convection.
 *
 * \param   pressure                (Pa)
 * \param   height                  (meters)
 * \param   u_wind                  (m/s)
 * \param   v_wind                  (m/s)
 * \param   N                       (length of arrays)
 * \param   mean_wind_layer_agl     {bottom, top}
 * \param   wind_shear_layer_agl    {bottom, top}
 * \param   leftMover               (default=false)
 * \param   pressureWeighted        (default=false)
 *
 * \return  {storm_u, storm_v}
 */
[[nodiscard]] WindComponents storm_motion_bunkers(
    const float pressure[], const float height[], const float u_wind[],
    const float v_wind[], const std::ptrdiff_t N,
    HeightLayer mean_wind_layer_agl, HeightLayer wind_shear_layer_agl,
    const bool leftMover = false, const bool pressureWeighted = false);

/**
 *  \author Kelton Halbert - NWS Storm Prediction Center
 *
 *  \brief Estimates supercell storm motion by using the<!--
 *  --> Bunkers et al. 2014 method.
 *
 *  Estimates supercell storm motion using the Bunkers et al. 2014 method
 *  described in the following paper:
 *      http://dx.doi.org/10.15191/nwajom.2014.0211
 *
 *  This method is parcel based, using a mean-wind vector defined as the
 *  pressure-weighted mean wind between the Effective Inflow Layer surface
 *  (see sharp::effective_inflow_layer) and 65% of the depth between that
 *  surface and the most-unstable parcel's Equilibrium Level. This method
 *  produces the same storm motion estimate for surface based supercells,
 *  and captures the motion of elevated supercells much better than
 *  the Bunkers 2000 method.
 *
 *  The input parameters eff_infl_lyr and mupcl (effective inflow layer
 *  bounds and the most unstable parcel, respectively) are required to be
 *  be precomputed and passed to this routine. These are expensive
 *  operations that are presumed to be computed at some other point
 *  in the analysis pipeline, so just pass those variables here.
 *
 * \param   pressure        (Pa)
 * \param   height          (meters)
 * \param   u_wind          (m/s)
 * \param   v_wind          (m/s)
 * \param   N               (length of arrays)
 * \param   eff_infl_lyr    {bottom, top}
 * \param   mupcl           (Precomputed parcel)
 * \param   leftMover       (default=false)
 *
 * \return  {storm_u, storm_v}
 */
[[nodiscard]] WindComponents storm_motion_bunkers(
    const float pressure[], const float height[], const float u_wind[],
    const float v_wind[], const std::ptrdiff_t N, PressureLayer eff_infl_lyr,
    const Parcel& mupcl, const bool leftMover = false);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes Entrainment CAPE using a previously lifted parcel.
 *
 * Computes Entrainment CAPE, or ECAPE, as described by Peters et al. 2023,
 * "An analytic formula for entraining CAPE in mid-latitude storm environments".
 *
 * \param   pressure        (Pa)
 * \param   height          (meters)
 * \param   temperature     (degK)
 * \param   mse_arr         ()
 * \param   u_wind          (m/s)
 * \param   v_wind          (m/s)
 * \param   N               (length of arrays)
 * \param   pcl             (Precomputed parcel)
 *
 * \return  ECAPE
 */
[[nodiscard]] float entrainment_cape(const float pressure[],
                                     const float height[],
                                     const float temperature[],
                                     const float mse_arr[],
                                     const float u_wind[], const float v_wind[],
                                     const std::ptrdiff_t N, Parcel* pcl);

[[nodiscard]] float energy_helicity_index(float cape, float helicity);

[[nodiscard]] float supercell_composite_parameter(float mu_cape, float eff_srh,
                                                  float eff_shear);

[[nodiscard]] float significant_tornado_parameter(Parcel pcl,
                                                  float lcl_hght_agl,
                                                  float storm_relative_helicity,
                                                  float bulk_wind_difference);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief compute the Significant Hail Parameter
 *
 * Compute the Significant Hail Parameter, which was developed to
 * delineate between significant (>=2" diameter) and non-significant
 * (<2" diameter) hail environments.
 *
 * The 0-6 km shear is confined to a range of 7-27 m/s, mixing ratio
 * is confined to a range of 11-13.6 g/kg, and the 500mb temperature
 * maximum is capped at -5.5C.
 *
 * After an initial calculation of SHIP, it undergoes the following
 * modifications:
 *  1) if MUCAPE < 1300 J/kg, SHIP = SHIP * (MUCAPE/1300)
 *  2) if 700-500 mb lapse rate < 5.8 C/km, SHIP = SHIP * (lr75/5.8)
 *  3) if the freezing level < 2400 m AGL, SHIP = SHIP * (fzl/2400)
 *
 * It is important to note the SHIP is NOT a forecast hail size.
 * Values of SHIP greater than 1 indicate a favorable environment
 * for significant hail, and values greater than 4 are considered
 * very high.
 *
 * \param   mu_pcl    precomputed Most Unstable Parcel
 * \param   lapse_rate_700_500mb    (K/km)
 * \param   tmpk_500mb              (K)
 * \param   freezing_lvl_agl        (meters)
 * \param   shear_0_6km             (m/s)
 *
 * \return Significant Hail Parameter
 */
[[nodiscard]] float significant_hail_parameter(const sharp::Parcel& mu_pcl,
                                               float lapse_rate_700_500mb,
                                               float tmpk_500mb,
                                               float freezing_lvl_agl,
                                               float shear_0_6km);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes the precipitable water vapor content over a layer
 *
 * Given a sharp::PressureLayer to integrate over, compute the preciptable
 * water from the given pressure and mixing ratio arrays.
 *
 * \param   layer           (Pa)
 * \param   presssure       (Pa)
 * \param   mixing_ratio    (unitless)
 * \param   N               (length of arrays)
 *
 * \return precipitable water (mm)
 */
[[nodiscard]] float precipitable_water(PressureLayer layer,
                                       const float pressure[],
                                       const float mixing_ratio[],
                                       const std::ptrdiff_t N);
/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Get the layer encompassing the lowest hail growth zone (-10 to
 * -30 C)
 *
 * Search for and return the sharp::PressuerLayer of the lowest altitude
 * hail growth zone. If none is found, the top and bottom pressure levels
 * are set to sharp::MISSING.
 *
 * \param    pressure    (Pa)
 * \param    temperature (K)
 *
 * \return   The top and bottom of the hail growth zone (Pa)
 */
[[nodiscard]] PressureLayer hail_growth_layer(const float pressure[],
                                              const float temperature[],
                                              const std::ptrdiff_t N);

}  // end namespace sharp

#endif  // SHARP_PARAMS_H
