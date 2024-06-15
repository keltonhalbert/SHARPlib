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
#ifndef __SHARP_PARAMS_H__
#define __SHARP_PARAMS_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
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
 * \param   pressure        (Pa)
 * \param   height          (meters)
 * \param   temperature     (degK)
 * \param   dewpoint        (degK)
 * \param   virtemp_arr     (degK)
 * \param   buoy_arr        (m/s^2)
 * \param   N               (length of arrays)
 * \param   cape_thresh     (J/kg; default=100.0)
 * \param   cinh_thresh     (J/kg; default=-250.0)
 * \param   mupcl           (store the most unstable parcel; defualt=nullptr)
 *
 * \return  {bottom, top}
 */
[[nodiscard]] PressureLayer effective_inflow_layer(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float virtemp_arr[], float buoy_arr[],
    const int N, const float cape_thresh = 100.0,
    const float cinh_thresh = -250.0, Parcel* mupcl = nullptr);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
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
    const float v_wind[], const int N, HeightLayer mean_wind_layer_agl,
    HeightLayer wind_shear_layer_agl, const bool leftMover = false,
    const bool pressureWeighted = false);

/**
 *  \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
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
    const float v_wind[], const int N, PressureLayer eff_infl_lyr,
    const Parcel* mupcl, const bool leftMover = false);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
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
                                     const int N, Parcel* pcl);

/**
 * \author Amelia Urquhart - OU-SoM
 *
 * \brief Computes Entrainment Rate using a previously lifted parcel.
 *
 * Computes Entrainment Rate based on a formula from John Peter's 
 * ECAPE_FUNCTIONS Python script. Method is copied from Kelton's ECAPE code
 * and altered to compute entrainment rate instead.
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
 * \return  Entrainment Rate (m^-1)
 */
[[nodiscard]] float entrainment_rate(const float pressure[],
                                     const float height[],
                                     const float temperature[],
                                     const float mse_arr[],
                                     const float u_wind[], const float v_wind[],
                                     const int N, Parcel* pcl);

[[nodiscard]] float energy_helicity_index(float cape, float helicity);

[[nodiscard]] float supercell_composite_parameter(float mu_cape, float eff_srh,
                                                  float eff_shear);

[[nodiscard]] float significant_tornado_parameter(Parcel pcl,
                                                  float lcl_hght_agl,
                                                  float storm_relative_helicity,
                                                  float bulk_wind_difference);

}  // end namespace sharp

#endif  // __SHARP_PARAMS_H__
