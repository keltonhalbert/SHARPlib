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

}  // namespace sharp

#endif  // SHARP_PARAMS_FIRE_H
