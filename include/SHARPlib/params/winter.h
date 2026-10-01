/**
 * \file
 * \brief Winter weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#ifndef SHARP_PARAMS_WINTER_H
#define SHARP_PARAMS_WINTER_H

#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>

#include <cstddef>

namespace sharp {
/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Get the layer encompassing the lowest dendritic growth zone (-12 to
 * -17 C)
 *
 * Search for and return the sharp::PressuerLayer of the lowest altitude
 * dendritic growth zone. If none is found, the top and bottom pressure levels
 * are set to sharp::MISSING.
 *
 * \param   pressure    (Pa)
 * \param   temperature (K)
 * \param   N           (length of arrays)
 *
 * \return   The top and bottom of the dendritic growth zone (Pa)
 */
[[nodiscard]] PressureLayer dendritic_layer(const float pressure[],
                                            const float temperature[],
                                            const std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes the Snow Squall Parameter
 *
 * The Snow Squall Parameter is a non-dimensional parameter that
 * combines several ingredients believed to be beneficial for
 * identifying snow squall environments by identifying the overlap
 * of low-level potential instability, sufficient moisture, and
 * strong low-level winds.
 *
 * References:
 * Banacos et al. 2014:
 * https://www.weather.gov/media/btv/research/Snow%20Squalls%20Forecasting%20and%20Hazard%20Mitigation.pdf
 *
 * \param    wetbulb_2m          (K)
 * \param    mean_relh_0_2km     (fraction)
 * \param    delta_thetae_0_2km  (K)
 * \param    mean_wind_0_2km     (m/s)
 *
 * \return The Snow Squall Parameter
 */
[[nodiscard]] float snow_squall_parameter(const float wetbulb_2m,
                                          const float mean_relh_0_2km,
                                          const float delta_thetae_0_2km,
                                          const float mean_wind_0_2km);

// ===========================================================================
// Precipitation type: the modified Bourgouin method (Birk et al. 2021)
// ===========================================================================

/**
 * \brief Default upper pressure limit for wet-bulb energies (Pa)
 *
 * Wet-bulb melting and refreezing energies come only from levels at
 * pressures at or above this value (250 hPa). This keeps stratospheric
 * temperatures out of the melting energy, and has no effect on realistic
 * tropospheric profiles.
 */
static constexpr float BOURGOUIN_PRESSURE_MIN = 25000.0f;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Wet-bulb melting and refreezing energies of the modified <!--
 * --> Bourgouin method.
 *
 * Energies are areas between the wet-bulb temperature profile and 0 C
 * (Birk et al. 2021, Eq. 1 with the wet-bulb temperature), reported as
 * positive values in J/kg. The near-surface cold layer is the lowest layer
 * colder than 0 C that has a warmer-than-0 C layer above it.
 *
 * Every field defaults to sharp::MISSING.
 *
 * References:
 * Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
 */
struct BourgouinEnergy {
    /**
     * \brief Total wet-bulb melting energy of the column (J/kg)
     *
     * The sum of all wet-bulb energy above 0 C. Used for snow and for
     * freezing rain or rain.
     */
    float melting_energy_total = MISSING;

    /**
     * \brief Wet-bulb melting energy above the near-surface cold layer (J/kg)
     *
     * Used for ice pellets. 0 when there is no near-surface cold layer.
     */
    float melting_energy_aloft = MISSING;

    /**
     * \brief Wet-bulb refreezing energy of the near-surface cold layer <!--
     * --> (J/kg, positive)
     *
     * 0 when there is no near-surface cold layer.
     */
    float refreezing_energy = MISSING;
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Probabilities of the four basic precipitation types.
 *
 * Each field is a fraction in [0, 1]. The probabilities are independent and
 * do not sum to 1. Two types with high probabilities describe a mix of the
 * two, not a contradiction.
 *
 * Every field defaults to sharp::MISSING.
 */
struct PrecipTypeProbabilities {
    /**
     * \brief Probability of rain (fraction)
     */
    float rain = MISSING;

    /**
     * \brief Probability of snow (fraction)
     */
    float snow = MISSING;

    /**
     * \brief Probability of freezing rain (fraction)
     */
    float freezing_rain = MISSING;

    /**
     * \brief Probability of ice pellets (fraction)
     */
    float ice_pellets = MISSING;
};

// ---------------------------------------------------------------------------
// Wet-bulb melting and refreezing energies from a sounding
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Precipitation generation layer from a sounding
// ---------------------------------------------------------------------------

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Finds the precipitation generation layer of the modified <!--
 * --> Bourgouin method.
 *
 * Splits the profile with sharp::for_each_threshold_layer into moist
 * layers, where relative humidity is above 75 %, and dry layers, where it
 * is below 75 %. Returns the highest moist layer deeper than 1000 m that
 * lies below the first dry layer deeper than 1500 m (Birk et al. 2021,
 * section 3e). The paper assumes that precipitation falling from above such
 * a dry layer sublimates, so the dry layer eliminates every layer above it,
 * even when it starts at the surface. Both depths are strict. A 1000 m
 * moist layer is not a generation layer, and a 1500 m dry layer eliminates
 * nothing. Layer boundaries are the linearly interpolated 75 % crossings.
 * The function walks up the profile once and stops at the first
 * eliminating dry layer.
 *
 * Relative humidity is over ice where the air temperature is below 0 C and
 * over liquid water otherwise, from sharp::relative_humidity_ice and
 * sharp::relative_humidity. The paper uses relative humidity over ice at
 * every temperature, so this deviates from it above 0 C. The two agree at
 * 0 C. For example, T = 283.15 K with Td = 280 K is 0.808 over liquid,
 * which is moist, but 0.733 over ice, which is dry.
 *
 * Levels at exactly 75 % continue the current layer and add their depth to
 * it, so a layer ends only where the relative humidity crosses 75 %. For
 * example, relative humidities of 0.8, 0.75, 0.75, and 0.8 at 0, 100, 2000,
 * and 2100 m form one 2100 m moist layer, which is a generation layer. The
 * paper defines moist as above 75 % and dry as below it, and doesn't say how
 * levels at exactly 75 % count. Read strictly, they belong to neither kind
 * of layer and end both. That reading turns the example into two 100 m
 * moist layers and no generation layer. The two readings differ only where
 * the relative humidity is exactly 0.75.
 *
 * min_depth is an opt-in filter for noisy, high-resolution profiles. Moist
 * and dry runs shallower than min_depth merge into the layers around them,
 * following the merge rule of sharp::for_each_threshold_layer, so a dry
 * layer's depth includes the shallow moist runs it absorbs. The default of
 * 0 merges nothing, which is the paper as written. With min_depth above 0,
 * a dry layer's depth is final once the moist run above it is min_depth
 * deep, so the walk stops there.
 *
 * A {sharp::MISSING, sharp::MISSING} result does not say why. The profile
 * may have no moist layer deeper than 1000 m, a deep dry layer under the
 * cloud may eliminate it (virga), or the moisture data may be missing.
 * Callers can't tell these cases apart from the result.
 *
 * Unless NO_QC is defined, the walk skips levels with a sharp::MISSING or
 * NaN height, temperature, or dewpoint, or a sharp::MISSING pressure, and
 * joins their valid neighbors with a straight line.
 *
 * Heights must be strictly increasing. This is not checked.
 *
 * References:
 * Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
 *
 * \param   pressure    (Pa)
 * \param   height      (meters)
 * \param   temperature (K)
 * \param   dewpoint    (K)
 * \param   N           (length of arrays)
 * \param   min_depth   (meters; 0 disables merging)
 *
 * \return  The precipitation generation layer (meters, AGL or MSL like
 *          height), or {sharp::MISSING, sharp::MISSING}
 */
[[nodiscard]] HeightLayer precipitation_generation_layer(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const std::ptrdiff_t N,
    const float min_depth = 0.0f);

// ---------------------------------------------------------------------------
// Probability of ice, and precipitation-type probabilities from energies
// ---------------------------------------------------------------------------

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Probability of ice in a precipitation generation layer
 *
 * Computes the probability of heterogeneous ice nucleation (ProbIce) from
 * the minimum air temperature T in the precipitation generation layer. It
 * uses the piecewise relation of Baumgardt et al. (2017), as reproduced by
 * Birk et al. (2021, Eq. 2 and Fig. 2):
 *
 * - 1 for T at or below -15 C
 * - 0 for T at or above -7 C
 * - otherwise (-0.065 T^4 - 3.1544 T^3 - 56.414 T^2 - 449.6 T - 1308) / 100,
 *   with T in C
 *
 * The polynomial does not meet the constant pieces. It gives 98.3 % at
 * -15 C and 0.81 % at -7 C, so ProbIce jumps by about 0.017 at -15 C and
 * about 0.008 at -7 C, as in the paper.
 *
 * \warning The input is in Kelvin, and Celsius input gives wrong results. A
 * value at or below 0 returns sharp::MISSING, and a positive value reads as
 * a Kelvin temperature far below -15 C and returns 1.
 *
 * sharp::MISSING or NaN input returns sharp::MISSING. The function checks
 * this in every build, including NO_QC builds.
 *
 * References:
 * Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
 *
 * Baumgardt et al. 2017:
 * https://ams.confex.com/ams/97Annual/webprogram/Paper313165.html
 *
 * \param   temperature     Minimum air temperature in the precipitation
 *                          generation layer (K)
 *
 * \return  Probability of ice (fraction)
 */
[[nodiscard]] float probability_of_ice(const float temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Precipitation-type probabilities from wet-bulb energies and the <!--
 * --> probability of ice (modified Bourgouin method).
 *
 * Computes the probabilities of rain, snow, freezing rain, and ice pellets
 * with the modified Bourgouin method, following the steps in the Appendix of
 * Birk et al. (2021). The function works in percent and returns fractions.
 * ProbIce is prob_ice in percent. Every clamp is to [0, 100] %.
 *
 * - Snow (Eq. 9): 1540 exp(-0.29 ME_total), clamped, then multiplied by
 *   ProbIce / 100 and clamped.
 * - Ice pellets (Eq. 8): when both ME_aloft and RE are above 0,
 *   2.3 RE - 42 ln(ME_aloft + 1) + 3, clamped; otherwise 0. Then multiplied
 *   by ProbIce / 100 and clamped.
 * - Freezing rain or rain (Eq. 7): -2.1 RE + 0.2 ME_total + 458, clamped
 *   first. When ME_total is below 5 J/kg, that clamped value is then
 *   multiplied by 0.2 ME_total. The result is
 *   (100 - ProbIce) + (ProbIce / 100) times that value, clamped.
 *
 * ME_total, ME_aloft, and RE are the fields of sharp::BourgouinEnergy.
 * Freezing rain and rain use the total melting energy, and ice pellets use
 * the melting energy above the near-surface cold layer.
 *
 * The four probabilities are independent and do not sum to 1. Freezing rain
 * and rain share one value. The function reports it as rain when
 * surface_wetbulb is above 0 C (273.15 K) and as freezing rain otherwise,
 * including at exactly 0 C, and sets the other to 0. A warm surface does not
 * suppress ice pellets.
 *
 * As in the paper, the ice pellet probability jumps from 0 to its Eq. 8
 * value as ME_aloft rises from 0. Just above 0, Eq. 8 gives 2.3 RE + 3 %
 * before the clamp and the ProbIce scaling.
 *
 * If any input is sharp::MISSING or NaN, all four probabilities are
 * sharp::MISSING. The function checks this in every build, including NO_QC
 * builds.
 *
 * References:
 * Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
 *
 * \param   energy          Wet-bulb melting and refreezing energies (J/kg)
 * \param   prob_ice        Probability of ice in the precipitation generation
 *                          layer (fraction; see sharp::probability_of_ice)
 * \param   surface_wetbulb Surface wet-bulb temperature (K)
 *
 * \return  {rain, snow, freezing_rain, ice_pellets} (fractions)
 */
[[nodiscard]] PrecipTypeProbabilities modified_bourgouin(
    const BourgouinEnergy& energy, const float prob_ice,
    const float surface_wetbulb);

// ---------------------------------------------------------------------------
// Precipitation type from a full sounding
// ---------------------------------------------------------------------------

}  // namespace sharp

#endif  // SHARP_PARAMS_WINTER_H
