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
 * \param theta_env (K)
 * \param beta      (unitless)
 *
 * \return Plume potential temperature
 */
[[nodiscard]] float pft_plume_potential_temperature(float theta_env,
                                                    float beta = 0.05);

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
 * https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
 *
 * \param theta_env (K)
 * \param mixr_env  (kg/kg)
 * \param beta      (unitless)
 * \param phi       (kg/(kg*K))
 *
 * \return Plume mixing ratio
 */
[[nodiscard]] float pft_plume_mixratio(float theta_env, float mixr_env,
                                       float beta = 0.05, float phi = 6.6e-5);

}  // namespace sharp

#endif  // SHARP_PARAMS_FIRE_H
