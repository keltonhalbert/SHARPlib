/**
 * \file
 * \brief Fire weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Prediction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#ifndef SHARP_PARAMS_FIRE_H
#define SHARP_PARAMS_FIRE_H

#include <SHARPlib/parcel.h>

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
 * \brief Computes an estimated rainfall terminal velocity
 *
 * Computes an estimated rainfall terminal velocity as
 * formulated by Klemp and Wilhelmson 1978.
 *
 * \param   rho (density)       (kg/m^3)
 * \param   rho_0 (density)     (kg/m^3)
 * \param   rainwater_mixratio  (kg/kg)
 *
 * \return  rainfall velocity   (m/s)
 */
[[nodiscard]] float rainfall_velocity(float rho, float rho_0,
                                      float rainwater_mixratio);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes an estimated rainfall efficiency
 *
 * Computes an estimated rainfall efficiency, using forumations of
 * rainfall terminal velocity and evaporation rate as in Klemp and
 * Wilhelmson 1978.
 *
 * Starting from the parcel's lifted condensation level, it integrates
 * downward to the surface while evaporating water content and updating
 * the velocity between steps. Once it reaches the surface, it computes
 * the final fraction of rainwater mixing ratio relative to the starting
 * rainwater mixing ratio.
 *
 * \param   pressure            (Pa)
 * \param   height              (meters)
 * \param   temperature         (K)
 * \param   mixr                (kg/kg)
 * \param   N                   (length of arrays)
 * \param   pcl                 (parcel to get LCL from)
 * \param   rainwater_mixratio  (kg/kg)
 */
float rainfall_efficiency(const float pressure[], const float height[],
                          const float temperature[], const float mixr[],
                          std::ptrdiff_t N, const Parcel& pcl,
                          float rainwater_mixratio);

}  // namespace sharp

#endif  // SHARP_PARAMS_FIRE_H
