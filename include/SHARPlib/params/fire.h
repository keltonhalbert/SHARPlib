/**
 * \file
 * \brief Fire weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#ifndef SHARP_PARAMS_FIRE_H
#define SHARP_PARAMS_FIRE_H

#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>

#include <cstddef>

namespace sharp {

/**
 * \author Kelton Halber - NSW Storm Prediction Center
 *
 * \brief Computes the equilibrium moisture content for fuel
 *
 * Compute the equilibrium moisture content for fuel as in
 * Simard (1968).
 *
 * \param temperature   (K)
 * \param rel_humidity  (fraction)
 *
 * \return EMC (fraction)
 */
[[nodiscard]] float equilibrium_moisture_content(float temperature,
                                                 float rel_humidity);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the Fosberg fire-weather index
 *
 * Compute the Fosberg Fire-Weather Index (FWWI) as in
 * Fosberg (1978).
 *
 * \param temperature    (K)
 * \param rel_humidity   (fraction, unitless)
 * \param wind_speed     (m/s)
 *
 * \return Fosberg Fire-Weather Index
 *
 */
[[nodiscard]] float fosberg_fire_index(float temperature, float rel_humidity,
                                       float wind_speed);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes a fire plume parcel's potential temperature.
 *
 * Compute the potential temperature of a hypothetical fire plume, as in
 * Tory et al. 2018. Beta is the plume buoyancy factor, with zero implying
 * complete dilution to the environment. Useful values are deemed to be
 * between 0 and 1e-1, though most realistic fires are with beta < 1e-2.
 * Theta_env is usually some sort of mixed/boundary-layer mean value of
 * potential temperature.
 *
 * References:
 * https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
 *
 * \param theta_env                     (K)
 * \param beta                          (unitless)
 *
 * \return Plume potential temperature  (K)
 */
[[nodiscard]] float pft_plume_potential_temperature(float theta_env,
                                                    float beta);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes a fire plume parcel's mixing ratio.
 *
 * Compute the mixing ratio of a hypothetical fire plume, as in
 * Tory et al. 2018. Beta is the plume buoyancy factor, with zero
 * implying complete dilution to the environment. Useful values are
 * deemed to be between 0 and 1e-1, though most realistic fires are
 * with beta < 1e-2. Theta_env and mixr_env are typically averaged
 * values within some sort of mixed/boundary-layer. Phi is the fire
 * moisture to potential temperature increment ratio.
 *
 * References:
 * Tory et al. 2018:
 * https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
 *
 * \param theta_env             (K)
 * \param mixr_env              (kg/kg)
 * \param beta                  (unitless)
 * \param phi                   (kg/(kg*K))
 *
 * \return Plume mixing ratio   (kg/kg)
 */
[[nodiscard]] float pft_plume_mixratio(float theta_env, float mixr_env,
                                       float beta, float phi);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes the Pyrocumulonimbus Firepower Threshold (PFT)
 *
 * Computes the Pyrocumulonimbus Firepower Threshold (PFT), or the minimum
 * amount of firepower required to generate pyrocumulonimbus clouds for a
 * given atmospheric profile. Reqires a PressureLayer to define a mixing layer
 * used to average values of potential temperature, mixing ratio, and wind
 * speed. The beta increment determines how to vary the plume buoyancy factor,
 * with smaller values resulting in more iteration steps. Phi is the fire
 * moisture to potential temperature increment ratio.
 *
 * Default values for beta_incr and phi are 0.005 and 6.67e-5, respectively.
 *
 * References:
 * Tory et al. 2018:
 * https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
 *
 * Tory et al. 2021:
 * https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml
 *
 * \param lifter                                (e.g. CM1, Wobus)
 * \param mix_layer                             {bottom, top}
 * \param pressure                              (Pa)
 * \param height                                (meters)
 * \param temperature                           (K)
 * \param mixratio                              (kg/kg)
 * \param virtemp                               (K)
 * \param uwin                                  (m/s)
 * \param vwin                                  (m/s)
 * \param potential_temperature                 (K)
 * \param pcl_vtmpk_arr                         (K)
 * \param pcl_buoy_arr                          (K)
 * \param N                                     (length of arrays)
 * \param pcl                                   A fire parcel to return
 * \param phi                                   (kg/(kg*K))
 * \param beta_incr                             (unitless)
 *
 * \return Pyrocumulonimbus Firepower Threshold (W)
 */
template <typename Lifter>
[[nodiscard]] float pyrocumulonimbus_firepower_threshold(
    Lifter& lifter, PressureLayer mix_layer, const float pressure[],
    const float height[], const float temperature[], const float mixratio[],
    const float virtemp[], const float uwin[], const float vwin[],
    const float potential_temperature[], float pcl_vtmpk_arr[],
    float pcl_buoy_arr[], std::ptrdiff_t N, Parcel* pcl = nullptr,
    float phi = 6.67e-5, float beta_incr = 0.005);

}  // namespace sharp

#endif  // SHARP_PARAMS_FIRE_H
