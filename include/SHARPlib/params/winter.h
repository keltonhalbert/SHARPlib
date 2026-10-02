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

#include <array>
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

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Wet-bulb melting and refreezing energies from a sounding <!--
 * --> (modified Bourgouin method).
 *
 * Computes the energy areas of Birk et al. (2021, Eq. 1) with the wet-bulb
 * temperature Tw: g (Tw - T0) / T0 integrated over height, with
 * T0 = 273.15 K. Layers are bounded by linearly interpolated 0 C
 * crossings, and the areas are trapezoids split exactly at those crossings
 * (see sharp::for_each_threshold_layer). The function reports every
 * energy as a positive value in J/kg.
 *
 * - melting_energy_total (ME_total): all energy above 0 C in the column.
 * - refreezing_energy (RE): the energy of the near-surface cold layer, the
 *   lowest layer below 0 C that has a layer above 0 C over it. 0 if there
 *   is no such layer.
 * - melting_energy_aloft (ME_aloft): all energy above 0 C over the top of
 *   the near-surface cold layer. 0 if there is no such layer.
 *
 * sharp::modified_bourgouin uses ME_total for snow and for freezing rain or
 * rain, and ME_aloft for ice pellets. With the default min_energy = 0,
 * ME_total equals ME_aloft whenever the lowest layer is cold. When the only
 * warm layer is at the surface (Fig. 1b), ME_total holds all of the
 * melting energy and ME_aloft is 0.
 *
 * With several warm layers, RE comes from the near-surface cold layer only,
 * and ME_aloft adds up every warm layer above it. Other cold layers never
 * enter the equations. Example, from the surface up: cold 100, warm 30,
 * cold 80, and warm 20 J/kg give RE = 100, ME_aloft = 50, and
 * ME_total = 50.
 *
 * A warm layer at the surface, below the near-surface cold layer (Fig. 1d),
 * adds to ME_total only. As in the paper, it does not suppress ice pellets.
 * Example, from the surface up: warm 150, cold 51, and warm 10 J/kg give
 * ME_total = 160, ME_aloft = 10, and RE = 51. With ProbIce = 1,
 * sharp::modified_bourgouin then gives rain 100 % and ice pellets about
 * 20 % (19.6 %).
 *
 * min_energy sets the energy a layer needs to stand on its own. Weaker
 * layers merge into their neighbors, following the merge rule of
 * sharp::for_each_threshold_layer with min_area = min_energy T0 / g. The
 * paper sets no minimum melting energy for the modified method, so the
 * default of 0 is the paper as written and merges nothing. A positive value
 * is an opt-in for noisy, high-resolution data, and a deviation from the
 * paper. It changes more than the onset of ice pellets:
 *
 * - A warm or cold layer with less than min_energy no longer splits the
 *   layers around it, so ME_aloft is either 0 or at least min_energy.
 * - A weak warm layer between two cold layers merges them into one
 *   near-surface cold layer, and RE includes both. Example: cold 100, warm
 *   1.99, cold 80, and warm 20 J/kg give RE = 180 with min_energy = 2, but
 *   RE = 100 with min_energy = 0. RE jumps as the weak layer crosses the
 *   threshold.
 * - ME_total still counts every warm layer, including the merged ones, but
 *   ME_aloft leaves out a warm layer merged into the near-surface cold
 *   layer. So ME_total and ME_aloft can differ over a cold surface. Example:
 *   cold 100, warm 1.99, cold 80, and warm 2.01 J/kg give ME_total = 4,
 *   ME_aloft = 2.01, and RE = 180 with min_energy = 2. The weak-melting
 *   taper of sharp::modified_bourgouin acts on ME_total.
 *
 * Consider min_energy for 1 Hz soundings and other high-resolution profiles
 * whose wet-bulb temperature stays near 0 C over some depth, such as an
 * isothermal melting layer or a surface layer close to 0 C. There, noise makes
 * the profile cross 0 C many times. A warm sliver over a cold surface layer
 * makes that layer a near-surface cold layer, ME_aloft and RE become positive,
 * and the ice pellet probability jumps from 0 to about 2.3 RE + 3 %. Start with
 * 2 J/kg, the melting-layer minimum of the original Bourgouin method. On the
 * three 1 Hz soundings in the SHARPlib test data, noise alone on a layer at
 * exactly 0 C makes wet-bulb layers of which 99 % hold less than 1.2 J/kg, and
 * the largest held 2.8 J/kg. With that noise added to a saturated 800 m layer
 * at -0.15 C, a wet-bulb temperature of -0.28 C, over a cold surface, spurious
 * ice pellets appeared in 25 % of trials with min_energy = 0 and in none with
 * 0.5 J/kg or more. With the noise doubled, it took 2 J/kg. The cost is that
 * warm layers weaker than min_energy no longer count, and noise can split a
 * slightly stronger layer into pieces that each fall below it. A 3.1 J/kg
 * melting layer was lost in 7 % of trials at 2 J/kg and in 31 % at 3 J/kg.
 *
 * min_energy never changes ME_total, which counts every warm layer. On 1 Hz
 * data, noise moved ME_total by up to about 1 J/kg in 95 % of trials. With
 * ME_total below 5 J/kg, that alone moved rain or freezing rain by 0.1 or more
 * in 9 to 32 % of trials, whatever the options.
 *
 * Only levels at pressures at or above pressure_min (by default
 * sharp::BOURGOUIN_PRESSURE_MIN, 250 hPa) are used, which keeps
 * stratospheric temperatures out of the melting energy. The column ends at
 * the highest such level, with nothing interpolated to pressure_min itself.
 * A pressure_min of 0 uses every level. ME_total depends on how far up the
 * data reach, up to that limit.
 *
 * The function never reads a wet-bulb temperature at a pressure below
 * pressure_min. A caller can compute the wet-bulb temperature only up to
 * pressure_min and fill the rest of the array with sharp::MISSING.
 *
 * Unless NO_QC is defined, the function skips levels whose wet-bulb
 * temperature is sharp::MISSING or NaN, and joins the valid levels on either
 * side with a straight line. With N < 2, or with fewer than 2 valid
 * levels at pressures at or above pressure_min, every energy is
 * sharp::MISSING, never 0, since zero energy would read as certain snow.
 * N < 2 returns before any array element is read.
 *
 * Height must be strictly increasing, and pressure must be valid and
 * strictly decreasing. This is not checked.
 *
 * References:
 * Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
 *
 * Bourgouin 2000:
 * https://doi.org/10.1175/1520-0434(2000)015%3C0583:AMTDPT%3E2.0.CO;2
 *
 * \param   pressure        (Pa)
 * \param   height          (m)
 * \param   wetbulb         (K)
 * \param   N               (length of arrays)
 * \param   min_energy      (J/kg; 0 disables merging)
 * \param   pressure_min    (Pa; levels at lower pressures are ignored)
 *
 * \return  {melting_energy_total, melting_energy_aloft, refreezing_energy}
 */
[[nodiscard]] BourgouinEnergy bourgouin_energy(
    const float pressure[], const float height[], const float wetbulb[],
    const std::ptrdiff_t N, const float min_energy = 0.0f,
    const float pressure_min = BOURGOUIN_PRESSURE_MIN);

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
 * Consider min_depth for 1 Hz soundings and other high-resolution profiles
 * where a moist layer is close to 1000 m deep, a dry layer is close to 1500 m
 * deep, or the relative humidity stays near 75 %. Noise there splits and joins
 * layers, which can move the result between sharp::MISSING, a low cloud, and a
 * higher, colder one. Start with 100 m. On the three 1 Hz soundings in the
 * SHARPlib test data, noise alone on a layer at exactly 75 % makes humidity
 * layers of which 90 % are under 110 m deep and 99 % under 200 m. With that
 * noise added to profiles whose layers were 150 to 250 m from those depths, the
 * generation layer changed in 19 to 29 % of trials with min_depth = 0, in
 * 3 to 12 % with 50 m, and in at most 0.2 % with 100 m. Larger values merge
 * real layers. One test sounding has a 324 m dry layer at the surface under a
 * 684 m moist layer. With 300 m, noise often thinned the dry layer below 300 m,
 * and the moist layer absorbed it and became a generation layer in 17 % of
 * trials, against under 1 % with 200 m.
 *
 * A {sharp::MISSING, sharp::MISSING} result does not say why. The profile
 * may have no moist layer deeper than 1000 m, a deep dry layer under the
 * cloud may eliminate it (virga), or the moisture data may be missing.
 * Callers can't tell these cases apart from the result.
 *
 * Unless NO_QC is defined, the walk skips levels with a sharp::MISSING or
 * NaN temperature or dewpoint, and joins their valid neighbors with a
 * straight line.
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

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Precipitation-type probabilities from a sounding (modified <!--
 * --> Bourgouin method).
 *
 * Computes the probabilities of rain, snow, freezing rain, and ice pellets
 * of Birk et al. (2021) from pressure, height, temperature, dewpoint, and
 * wet-bulb temperature profiles. It runs these steps, and calling them
 * yourself with the same options gives the same result:
 *
 * 1. sharp::precipitation_generation_layer with min_depth finds the
 *    precipitation generation layer.
 * 2. sharp::layer_min finds the minimum air temperature in that layer, and
 *    sharp::probability_of_ice turns it into ProbIce.
 * 3. sharp::bourgouin_energy with min_energy and pressure_min computes the
 *    wet-bulb melting and refreezing energies of the whole column, not just
 *    the generation layer.
 * 4. The surface wet-bulb temperature is that of the lowest level. Above
 *    0 C, the liquid probability is rain, otherwise freezing rain.
 * 5. The sharp::modified_bourgouin overload that takes energies combines
 *    them into the four probabilities, which are independent and do not
 *    sum to 1.
 *
 * With the defaults (min_depth = 0, min_energy = 0, and pressure_min =
 * sharp::BOURGOUIN_PRESSURE_MIN), the function follows the paper except
 * in three ways:
 *
 * - The energies use only levels at pressures at or above 250 hPa. This
 *   has no effect on realistic tropospheric profiles.
 * - Relative humidity is over ice below 0 C and over liquid water
 *   otherwise, where the paper uses relative humidity over ice at every
 *   temperature. For example, T = 283.15 K with Td = 280 K is moist over
 *   liquid (0.808) but dry over ice (0.733).
 * - Levels at exactly 75 % relative humidity continue the current layer,
 *   where the paper puts them in neither the moist nor the dry class.
 *
 * Positive min_depth and min_energy are opt-ins for noisy, high-resolution
 * data and further deviations from the paper. min_depth applies only to
 * the generation layer, and min_energy and pressure_min only to the
 * energies. For 1 Hz soundings, start with min_depth = 100 m and
 * min_energy = 2 J/kg. sharp::precipitation_generation_layer and
 * sharp::bourgouin_energy describe their effects, the measurements behind
 * these values, and what the options cost.
 *
 * Every probability is sharp::MISSING when there is no generation layer.
 * The result does not say why. The profile may have no moist layer deeper
 * than 1000 m, a deep dry layer under the cloud may eliminate it (virga),
 * or the moisture data may be missing. In particular, a cloud 1 km deep or
 * less, such as a drizzle cloud, gives sharp::MISSING. A caller who wants
 * a result for such a cloud can call sharp::bourgouin_energy and then the
 * sharp::modified_bourgouin overload that takes energies, with
 * prob_ice = 0, which treats the cloud as having no ice.
 *
 * ME_total, and with it snow and freezing rain or rain, depends on how far
 * up the data reach, up to pressure_min. The wet-bulb temperature at
 * pressures below pressure_min does not affect the result, so a caller can
 * compute it only up to pressure_min and fill the rest of the array with
 * sharp::MISSING.
 *
 * Unless NO_QC is defined, the steps skip missing levels as their own
 * documentation describes, and the surface wet-bulb temperature is that of
 * the lowest level whose wet-bulb temperature is not sharp::MISSING or
 * NaN. With NO_QC, the surface wet-bulb temperature is wetbulb[0], and
 * keeping sharp::MISSING and NaN out of the profiles is the caller's job.
 * With N < 2, every probability is sharp::MISSING, and the function
 * returns before reading any array element.
 *
 * The profiles must start at the surface. Height must be strictly
 * increasing, and pressure must be valid and strictly decreasing. This is
 * not checked.
 *
 * References:
 * Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
 *
 * \param   pressure        (Pa)
 * \param   height          (meters)
 * \param   temperature     (K)
 * \param   dewpoint        (K)
 * \param   wetbulb         (K)
 * \param   N               (length of arrays)
 * \param   min_depth       (meters; generation layer only; 0 disables
 *                          merging)
 * \param   min_energy      (J/kg; energies only; 0 disables merging)
 * \param   pressure_min    (Pa; energies only; levels at lower pressures
 *                          are ignored)
 *
 * \return  {rain, snow, freezing_rain, ice_pellets} (fractions)
 */
[[nodiscard]] PrecipTypeProbabilities modified_bourgouin(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float wetbulb[], const std::ptrdiff_t N,
    const float min_depth = 0.0f, const float min_energy = 0.0f,
    const float pressure_min = BOURGOUIN_PRESSURE_MIN);

// ===========================================================================
// Precipitation type: the spectral bin classifier (Reeves et al. 2016)
// ===========================================================================

// ---------------------------------------------------------------------------
// Result types and the drop-size distribution
// ---------------------------------------------------------------------------

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Precipitation type from the spectral bin classifier.
 *
 * The seven categories of the 2023 version of the classifier, which adds
 * rain mixed with ice pellets to the six categories of Reeves et al. (2016).
 * The integer values are stable, for gridded output, and 0 is unused.
 * PrecipType::missing converts to sharp::MISSING as a float.
 *
 * References:
 * Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
 */
enum class PrecipType : int {
    /**
     * \brief No classification
     */
    missing = -9999,

    /**
     * \brief Rain (RA)
     */
    rain = 1,

    /**
     * \brief Snow (SN)
     */
    snow = 2,

    /**
     * \brief Rain and snow (RASN)
     */
    rain_snow = 3,

    /**
     * \brief Freezing rain (FZRA)
     */
    freezing_rain = 4,

    /**
     * \brief Ice pellets (PL)
     */
    ice_pellets = 5,

    /**
     * \brief Freezing rain and ice pellets (FZRAPL)
     */
    freezing_rain_ice_pellets = 6,

    /**
     * \brief Rain and ice pellets (RAPL)
     */
    rain_ice_pellets = 7,
};

/**
 * \brief Capacity of a sharp::SpectralBinDSD (bins)
 *
 * The classifier keeps its per-bin state in fixed-size blocks of this many
 * bins on the stack.
 */
static constexpr std::ptrdiff_t SBC_MAX_BINS = 64;

/**
 * \brief Default ice nucleation temperature of the spectral bin <!--
 * --> classifier (K)
 *
 * -6 C, as in Reeves et al. (2016) and the reference code.
 */
static constexpr float SBC_ICE_NUCLEATION_TEMPERATURE = 267.15f;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Result of the spectral bin classifier for one column.
 *
 * Every field defaults to missing: PrecipType::missing and sharp::MISSING.
 */
struct SpectralBinResult {
    /**
     * \brief Precipitation type at the surface
     */
    PrecipType precip_type = PrecipType::missing;

    /**
     * \brief Liquid share of the precipitation mass reaching the surface <!--
     * --> (fraction)
     *
     * It is not a probability. 0.5 means a mix of liquid and ice, not even
     * odds.
     */
    float liquid_fraction = MISSING;

    /**
     * \brief Height of the lowest supercooled liquid water (m AGL)
     *
     * 0 when the surface type is freezing rain, alone or with ice pellets.
     * sharp::MISSING when the classifier finds no supercooled liquid.
     */
    float supercooled_liquid_height = MISSING;
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A drop-size distribution and the per-bin constants that the <!--
 * --> spectral bin classifier reads.
 *
 * Build one with sharp::spectral_bin_dsd or sharp::spectral_bin_dsd_default,
 * and reuse it for every column. It holds up to sharp::SBC_MAX_BINS bins in
 * fixed-size arrays and allocates nothing. Only the first nbins() elements
 * of each array are bins, and the rest are 0. A default-constructed object
 * has nbins() == 0 and is invalid, like every distribution the builder
 * rejects.
 *
 * Each constant depends only on the diameters and the riming factor, and
 * comes from the 2023 Python reference. Each accessor names its variable in
 * the reference. Diameters are in mm, masses in g, densities in g cm^-3, and
 * fall speeds in m/s, the units of the reference's formulas.
 *
 * The snow_ constants describe the snow that a frozen cloud top produces,
 * and liquid_aa() and liquid_bb() belong to a liquid cloud top. The liquid_
 * aspect-ratio factors apply to a bin whose class at the level above takes
 * the raindrop aspect ratio in the reference. Those are its liquid_ar
 * classes, which include ice pellets. Other bins use an aspect ratio of 0.8.
 *
 * References:
 * Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
 *
 * Python reference (sbc_alg_2023Aug31.py): D. Tripp, 2023
 */
struct SpectralBinDSD {
    /**
     * \brief Number of bins, or 0 for an invalid distribution
     */
    [[nodiscard]] std::ptrdiff_t nbins() const { return m_nbins; }

    /**
     * \brief Degree of riming of the snow from a frozen cloud top (1 to 5)
     *
     * sharp::MISSING for an invalid distribution.
     */
    [[nodiscard]] float rime_factor() const { return m_rime_factor; }

    /**
     * \brief Melted diameter of each bin, exactly as passed (mm)
     *
     * particle_diam in the reference.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& diameter() const {
        return m_diameter;
    }

    /**
     * \brief Number concentration of each bin, exactly as passed
     *
     * particle_count in the reference.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& concentration()
        const {
        return m_concentration;
    }

    /**
     * \brief Mass of each hydrometeor, (pi / 6) 10^-3 D^3 (g)
     *
     * mass_hydro in the reference.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& mass() const {
        return m_mass;
    }

    /**
     * \brief Raindrop fall speed near the ground (m/s)
     *
     * velocity_RA_sfc in the reference: -0.1021 + 4.932 D - 0.9551 D^2 +
     * 0.07932 D^3 - 0.002362 D^4, raised to at least 0.01 m/s and to the
     * fall speed of the bin below it, so it never decreases with diameter.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& rain_fall_speed()
        const {
        return m_rain_fall_speed;
    }

    /**
     * \brief Ice pellet fall speed near the ground (m/s)
     *
     * velocity_PL_sfc in the reference: 0.2259 + 1.5954 D_i - 0.0405 D_i^2,
     * with the frozen-drop diameter D_i = D (1 / 0.917)^(1/3).
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& pellet_fall_speed()
        const {
        return m_pellet_fall_speed;
    }

    /**
     * \brief Raindrop fall speed of Foote and du Toit (1969) (m/s)
     *
     * -0.193 + 4.96 D - 0.904 D^2 + 0.0566 D^3, which the reference scales
     * by exp(z / 20 km) below a liquid cloud top. It is negative for
     * diameters below about 0.039 mm, as in the reference.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>&
    foote_du_toit_fall_speed() const {
        return m_foote_du_toit_fall_speed;
    }

    /**
     * \brief Capacitance factor for the raindrop aspect ratio
     *
     * The reference's capac is this factor times the diameter of the level
     * above, plus 0.2 fw: 0.5 a^(-1/3) e / asin(e) 0.8, with
     * e = 0.198997487421 and the raindrop aspect ratio
     * a = min(0.9951 + 0.02510 D - 0.03644 D^2 + 0.005303 D^3 -
     * 0.0002492 D^4, 1).
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>&
    liquid_capacitance_factor() const {
        return m_liquid_capacitance_factor;
    }

    /**
     * \brief Length factor for the raindrop aspect ratio
     *
     * The reference's unknwn_leng is this factor times the diameter of the
     * level above: (2 + a^2 / e ln((1 + e) / (1 - e))) / (4 a^(1/3)), with
     * e and a as for liquid_capacitance_factor().
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& liquid_length_factor()
        const {
        return m_liquid_length_factor;
    }

    /**
     * \brief Diameter of the snow at a frozen cloud top (mm)
     *
     * diam_melt_snow at a frozen cloud top in the reference:
     * 2.29 f_rim^-0.48 D^1.443 for D at or above
     * 0.154 f_rim^1.08 0.5^-0.75 mm, and 0.5^(-1/3) D below it.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& snow_diameter()
        const {
        return m_snow_diameter;
    }

    /**
     * \brief Density of the dry snow at a frozen cloud top (g cm^-3)
     *
     * density_drySnow at a frozen cloud top in the reference:
     * 0.178 f_rim D_s^-0.922 for the snow diameter D_s, capped at 0.5.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& snow_density()
        const {
        return m_snow_density;
    }

    /**
     * \brief Fall-speed coefficient aa of the snow at a frozen cloud top
     *
     * uknwn_aa in the reference: 1.26 rho^(-1/3), from the snow density
     * before the cap at 0.5. Melting snow falls at the raindrop fall speed,
     * corrected for air density, divided by aa - bb fw (1 + fw).
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& snow_aa() const {
        return m_snow_aa;
    }

    /**
     * \brief Fall-speed coefficient bb of the snow at a frozen cloud top
     *
     * uknwn_bb in the reference: (aa - 1) / 2.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& snow_bb() const {
        return m_snow_bb;
    }

    /**
     * \brief Fall-speed coefficient aa at a liquid cloud top
     *
     * uknwn_aa at a liquid cloud top in the reference: 1.26 rho^(-1/3), with
     * the dry-snow density of the melted diameter, 0.178 f_rim D^-0.922,
     * uncapped.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& liquid_aa() const {
        return m_liquid_aa;
    }

    /**
     * \brief Fall-speed coefficient bb at a liquid cloud top
     *
     * uknwn_bb at a liquid cloud top in the reference: (aa - 1) / 2.
     */
    [[nodiscard]] const std::array<float, SBC_MAX_BINS>& liquid_bb() const {
        return m_liquid_bb;
    }

   private:
    std::ptrdiff_t m_nbins = 0;
    float m_rime_factor = MISSING;
    std::array<float, SBC_MAX_BINS> m_diameter{};
    std::array<float, SBC_MAX_BINS> m_concentration{};
    std::array<float, SBC_MAX_BINS> m_mass{};
    std::array<float, SBC_MAX_BINS> m_rain_fall_speed{};
    std::array<float, SBC_MAX_BINS> m_pellet_fall_speed{};
    std::array<float, SBC_MAX_BINS> m_foote_du_toit_fall_speed{};
    std::array<float, SBC_MAX_BINS> m_liquid_capacitance_factor{};
    std::array<float, SBC_MAX_BINS> m_liquid_length_factor{};
    std::array<float, SBC_MAX_BINS> m_snow_diameter{};
    std::array<float, SBC_MAX_BINS> m_snow_density{};
    std::array<float, SBC_MAX_BINS> m_snow_aa{};
    std::array<float, SBC_MAX_BINS> m_snow_bb{};
    std::array<float, SBC_MAX_BINS> m_liquid_aa{};
    std::array<float, SBC_MAX_BINS> m_liquid_bb{};

    friend SpectralBinDSD spectral_bin_dsd(const float diameter[],
                                           const float concentration[],
                                           const std::ptrdiff_t nbins,
                                           const float rime_factor);
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Builds a drop-size distribution for the spectral bin classifier.
 *
 * Validates the bins and precomputes the per-bin constants of
 * sharp::SpectralBinDSD. The diameters are melted diameters in mm, as in the
 * paper and the reference code. Only the ratios of the concentrations
 * matter, so any unit works, and a bin may hold 0.
 *
 * The result is invalid, with nbins() == 0, unless:
 *
 * - 1 <= nbins <= sharp::SBC_MAX_BINS
 * - the diameters are finite, positive, and strictly increasing
 * - the raindrop aspect-ratio fit of the reference is positive at every
 *   diameter, which holds below about 12.16 mm
 * - the concentrations are finite and non-negative, and at least one is
 *   positive
 * - rime_factor is in [1, 5]
 *
 * The function checks this in every build, including NO_QC builds. With
 * nbins outside [1, sharp::SBC_MAX_BINS], it returns before reading any
 * array element.
 *
 * rime_factor is the degree of riming of the snow that a frozen cloud top
 * produces, from 1 (none) to 5 (graupel). It replaces the reference's fixed
 * value of 1. Melting layers below a refreezing layer use 5 instead, as in
 * the reference.
 *
 * References:
 * Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
 *
 * Python reference (sbc_alg_2023Aug31.py, run_sbc.py): D. Tripp, 2023
 *
 * C++ MRMS code (sbcmodel_core.cc): A. Rosenow and D. Tripp
 *
 * \param   diameter        Melted diameter of each bin (mm)
 * \param   concentration   Number concentration of each bin (any unit)
 * \param   nbins           (length of arrays)
 * \param   rime_factor     Degree of riming (1 to 5, unitless)
 *
 * \return  The drop-size distribution, or one with nbins() == 0
 */
[[nodiscard]] SpectralBinDSD spectral_bin_dsd(const float diameter[],
                                              const float concentration[],
                                              const std::ptrdiff_t nbins,
                                              const float rime_factor = 1.0f);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief The default drop-size distribution of the spectral bin classifier.
 *
 * The 4-bin distribution of the 2023 Python reference (run_sbc.py), with
 * rime_factor = 1:
 *
 * | Diameter (mm) | 0.05    | 0.75    | 1.45    | 2.15    |
 * |---------------|---------|---------|---------|---------|
 * | Concentration | 55.1843 | 146.647 | 11.6891 | 3.60886 |
 *
 * The reference interpolates its table of the DSD25 distribution of Reeves
 * et al. (2016), 0.05 to 1.65 mm, every 0.7 mm. The 2.15 mm bin lies past
 * the end of the table and takes its last value. The paper uses DSD25 as
 * measured, with 18 bins 0.1 mm apart and a largest diameter of 1.84 mm.
 * The C++ MRMS code (version 2.0.3) caps the bins at 1.85 mm: 0.05, 0.65,
 * 1.25, and 1.85 mm, with concentrations 55.1843, 206.606, 25.4924, and
 * 3.60886. Either can be passed to sharp::spectral_bin_dsd.
 *
 * References:
 * Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
 *
 * Python reference (run_sbc.py): D. Tripp, 2023
 *
 * \return  The default drop-size distribution
 */
[[nodiscard]] SpectralBinDSD spectral_bin_dsd_default();

// ---------------------------------------------------------------------------
// Cloud top from a sounding
// ---------------------------------------------------------------------------

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Cloud-top height of the spectral bin classifier from a sounding.
 *
 * Finds the cloud top of the 2023 version of the classifier from the
 * dewpoint depression T - Td and the relative humidity, searching from the
 * highest level down. A cloud level has T - Td of at most 6 K and relative
 * humidity above 0.60.
 *
 * 1. The cloud top is the highest cloud level.
 * 2. If a level below it has T - Td above 10 K or relative humidity below
 *    0.40, the cloud top moves down to the highest cloud level at or below
 *    the driest level. The driest level has the largest T - Td at or below
 *    the cloud top, and the highest of tied levels wins. With no cloud
 *    level there, the cloud top stays.
 * 3. With no cloud level at all, the cloud top is the highest level, other
 *    than the highest level of the profile, with relative humidity of at
 *    least 0.80.
 * 4. Otherwise there is no cloud, and the result is sharp::MISSING.
 *
 * The thresholds compare as written. T - Td of exactly 6 K makes a cloud
 * level, and exactly 10 K is not dry. Relative humidity of exactly 0.60
 * does not make a cloud level, exactly 0.40 is not dry, and exactly 0.80
 * passes step 3.
 *
 * The Python reference and the C++ MRMS code use this rule, and it departs
 * from the paper. Reeves et al. (2016) put the cloud top at the level of
 * highest relative humidity when the column maximum is above 80 %, and
 * otherwise classify rain or snow from the surface wet-bulb temperature.
 * Here, no cloud gives sharp::MISSING, as in the Python reference. This port
 * leaves out the rain and snow fallback of the C++ MRMS code.
 *
 * The C++ MRMS code differs from the Python reference in two ways. This
 * function follows the Python, which the authors consider authoritative:
 *
 * - A negative T - Td, from supersaturated data, is used as it is. The C++
 *   code raises it to 0, which can change the driest level.
 * - The highest level can be the cloud top in steps 1 and 2. The C++ code
 *   treats a cloud top there as no cloud.
 *
 * Relative humidity is an input, as in the reference, which reads it from
 * the model. That is why step 3 can fire. A relative humidity of 0.80 or
 * more computed from T and Td would mean T - Td under 6 K, which already
 * makes a cloud level.
 *
 * The rule does not read pressure. The parameter keeps the argument list of
 * the classifier.
 *
 * Unless NO_QC is defined, the search skips levels whose temperature,
 * dewpoint, or relative humidity is sharp::MISSING or NaN, and the highest
 * level of the profile is the highest valid one. With NO_QC, keeping
 * sharp::MISSING and NaN out of the profiles is the caller's job. With
 * N < 1, the function returns sharp::MISSING before reading any array
 * element.
 *
 * The profiles must start at the surface. This is not checked.
 *
 * References:
 * Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
 *
 * Python reference (run_sbc.py): D. Tripp, 2023
 *
 * C++ MRMS code (topCalc.cc): A. Rosenow and D. Tripp
 *
 * \param   pressure    (Pa; not read)
 * \param   height      (m)
 * \param   temperature (K)
 * \param   dewpoint    (K)
 * \param   relh        Relative humidity (fraction)
 * \param   N           (length of arrays)
 *
 * \return  The cloud-top height (m, AGL or MSL like height), or
 *          sharp::MISSING
 */
[[nodiscard]] float spectral_bin_cloud_top(const float pressure[],
                                           const float height[],
                                           const float temperature[],
                                           const float dewpoint[],
                                           const float relh[],
                                           const std::ptrdiff_t N);

// ---------------------------------------------------------------------------
// Precipitation type from a given cloud top: pre-classifier
// ---------------------------------------------------------------------------

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Precipitation type from the spectral bin classifier, given a <!--
 * --> cloud top.
 *
 * The column runs from the surface, the lowest valid level, up to the
 * highest valid level at or below cloud_top. The function never reads the
 * temperature, dewpoint, relative humidity, or wet-bulb temperature of a
 * level above cloud_top, so a caller can compute the wet-bulb temperature
 * only up to the cloud top and fill the rest of the array with
 * sharp::MISSING. Heights above ground level (AGL) are height minus the
 * height of the surface. The pre-classifier of the 2023 version of the
 * algorithm (run_sbc.py) then applies these rules in order, with Tw the
 * wet-bulb temperature, Tice the ice nucleation temperature, and
 * 0 C = 273.15 K:
 *
 * 1. If the maximum Tw of the column is at or below 0 C: snow (SN) if Tw at
 *    the cloud top is below Tice and the minimum Tw from the level nearest
 *    3000 m AGL down to the surface is below Tice, otherwise freezing rain
 *    (FZRA). Of two levels equally near 3000 m AGL, the higher one counts,
 *    as in the reference. With a cloud top below 3000 m AGL, the level
 *    nearest 3000 m AGL is the cloud top.
 * 2. Otherwise, if Tw at the cloud top is above Tice and Tw at the surface
 *    is below 0 C: FZRA.
 * 3. Otherwise, if the minimum Tw of the column is above 0 C: rain (RA).
 *
 * Each comparison uses the operator of the reference against a float
 * constant, so a Tw of exactly 273.15f is at 0 C: it counts as subfreezing
 * in rule 1 and does not fire rule 2.
 *
 * The rules are unpublished. The comments of the C++ MRMS code credit
 * H. Reeves, who developed them from tests on a large dataset in the study
 * published as Reeves et al. (2023, Weather and Forecasting). They depart
 * from Fig. 2 of Reeves et al. (2016):
 *
 * - For a column at or below 0 C, the paper gives SN when Tw at the cloud
 *   top is at or below Tice, and FZRA otherwise. Rule 1 also requires the
 *   minimum Tw from 3000 m AGL down to be below Tice, and gives FZRA,
 *   "non-classical freezing rain", for a cloud top at exactly Tice.
 * - Every cloud top warmer than Tice over a subfreezing surface gives FZRA,
 *   even with a deep layer colder than Tice below a warm layer, where the
 *   paper integrates the microphysics and can give ice pellets.
 * - A column above 0 C everywhere gives RA without the microphysics.
 *
 * The reference reports no liquid fraction or supercooled-liquid height
 * for a column the pre-classifier decides. This function returns:
 *
 * | Category | liquid_fraction | supercooled_liquid_height |
 * |----------|-----------------|---------------------------|
 * | SN       | 0               | sharp::MISSING            |
 * | FZRA     | 1               | 0 m                       |
 * | RA       | 1               | sharp::MISSING            |
 *
 * The result is missing (PrecipType::missing, with sharp::MISSING fields)
 * when:
 *
 * - N < 2. The function then returns before reading any input array
 *   element.
 * - dsd is invalid (nbins() == 0).
 * - ice_nucleation_temperature is sharp::MISSING, NaN, or not above 0.
 * - cloud_top is sharp::MISSING or NaN, as when there is no cloud.
 * - The column has fewer than 2 valid levels, including a cloud top below
 *   the surface. This departs from the reference, which pre-classifies a
 *   column whose cloud top is the surface level.
 *
 * Unless NO_QC is defined, a level is valid when its temperature, dewpoint,
 * relative humidity, and wet-bulb temperature are neither sharp::MISSING
 * nor NaN, and the function skips the other levels. With NO_QC, every level
 * is valid, and keeping sharp::MISSING and NaN out of the profiles is the
 * caller's job. Pressure and height are never checked.
 *
 * When liquid_fraction_profile is not null, the function writes all of its
 * N x nbins elements, row-major by level in input order: [k * nbins + j] is
 * level k and bin j. It holds sharp::MISSING outside the integrated column,
 * at skipped levels, and everywhere for a column that the pre-classifier
 * decides or that gives missing.
 *
 * The profiles must start at the surface, and height must be strictly
 * increasing. This is not checked.
 *
 * References:
 * Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
 *
 * Python reference (run_sbc.py): D. Tripp, 2023
 *
 * C++ MRMS code (sbcmodel_core.cc): A. Rosenow and D. Tripp
 *
 * \param   pressure                    (Pa)
 * \param   height                      (m)
 * \param   temperature                 (K)
 * \param   dewpoint                    (K)
 * \param   relh                        Relative humidity over liquid water
 *                                      (fraction)
 * \param   wetbulb                     Wet-bulb temperature (K)
 * \param   N                           (length of arrays)
 * \param   cloud_top                   Cloud-top height, AGL or MSL like
 *                                      height (m)
 * \param   dsd                         Drop-size distribution, with
 *                                      diameters in mm (see
 *                                      sharp::spectral_bin_dsd)
 * \param   ice_nucleation_temperature  Tice (K)
 * \param   liquid_fraction_profile     Optional output, N x nbins liquid
 *                                      fractions, or nullptr (fraction)
 *
 * \return  {precip_type, liquid_fraction, supercooled_liquid_height}
 */
[[nodiscard]] SpectralBinResult spectral_bin_classifier(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float relh[], const float wetbulb[],
    const std::ptrdiff_t N, const float cloud_top, const SpectralBinDSD& dsd,
    const float ice_nucleation_temperature = SBC_ICE_NUCLEATION_TEMPERATURE,
    float liquid_fraction_profile[] = nullptr);

// ---------------------------------------------------------------------------
// Microphysics: frozen cloud tops and melting
// ---------------------------------------------------------------------------

/**
 * \fn SpectralBinResult spectral_bin_classifier( \
 *     const float[], const float[], const float[], const float[], \
 *     const float[], const float[], const std::ptrdiff_t, const float, \
 *     const SpectralBinDSD&, const float, float[])
 *
 * <b>Microphysics: frozen cloud tops and melting.</b> A column that the
 * pre-classifier does not decide runs the microphysics of the reference
 * (classify in sbc_alg_2023Aug31.py) for every bin of dsd, level by level
 * from the cloud top down to the surface. A layer spans two adjacent valid
 * levels. The function first counts the 0 C crossings Nc between adjacent
 * levels, where a level at exactly 0 C counts as subfreezing, and finds the
 * first crossing.
 *
 * With Tw at or below Tice at the cloud top, every bin starts as snow, with
 * the diameter, density, and fall-speed coefficients aa and bb of
 * sharp::SpectralBinDSD. Above the first crossing the snow falls unchanged.
 * Its fall speed is the raindrop fall speed near the ground times
 * sqrt(rho_0 / rho) / aa, with the air density rho and
 * rho_0 = 1.292e-3 g cm^-3. At each level with Tw at or above 0 C, a bin
 * that is not all liquid melts. Its liquid fraction grows with the heat
 * flux from the air across the layer above, which depends on Tw and the
 * relative humidity. A bin melts with the aspect ratio of a raindrop when
 * its class at the level above is liquid-like, ice pellets included.
 * Otherwise its aspect ratio is 0.8. The class of a bin is rain, snow, or
 * both, from the liquid share of its mass flux and a threshold of 0.15.
 *
 * At the surface, the liquid and ice mass fluxes are Pw = sum(m0 fw v N)
 * and Pi = sum(m0 (1 - fw) v N) / 0.917 over the bins, with the mass m0,
 * liquid fraction fw, fall speed v, and concentration N of each bin.
 * liquid_fraction is Pw / (Pw + Pi). The reference rounds it to 0.1 %, and
 * this function does not. A surface with Tw above 0 C is warm:
 *
 * | Nc  | Surface | Category by liquid_fraction                      |
 * |-----|---------|--------------------------------------------------|
 * | 1   | warm    | RA above 0.85, SN below 0.60, otherwise RASN     |
 * | > 1 | warm    | PL below 0.15, RA above 0.85, otherwise RAPL     |
 * | any | cold    | PL below 0.15, FZRA above 0.85, otherwise FZRAPL |
 *
 * FZRA and FZRAPL set supercooled_liquid_height to 0 m. This decision
 * departs from Fig. 2 of Reeves et al. (2016), which compares Pw and Pi
 * with a ratio of 0.15, has no RAPL, and gives RA, RASN, or PL, never SN,
 * over a warm surface.
 *
 * The function keeps these quirks of the reference:
 *
 * - The air density p / (R_d Tw), with R_d = 287 J kg^-1 K^-1, uses the
 *   wet-bulb temperature.
 * - aa and bb keep their cloud-top values all the way down. The reference
 *   recomputes them at melting levels but never reads the new values.
 * - The surface ice flux uses an ice density of 0.917 g cm^-3, and the
 *   class of each bin 0.918.
 * - A quantity that the reference does not set at a level reads back as 0
 *   at the next level, because its arrays start at 0. The snow of a frozen
 *   cloud top has a liquid fraction of 0, and a bin that has melted
 *   completely has no snow mass and no ice or snow volume.
 * - The reference weights the mass flux of each bin by the bin width and by
 *   the ratio of its fall speed at the surface to its fall speed at the
 *   level. That ratio is always 1 where it is used, and both factors cancel
 *   in every ratio of fluxes, so the function leaves them out.
 *
 * liquid_fraction_profile receives the liquid fraction of every bin at
 * every level of the integration.
 */

// ---------------------------------------------------------------------------
// Microphysics: refreezing
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Microphysics: liquid cloud tops
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Precipitation type from a full sounding
// ---------------------------------------------------------------------------

}  // namespace sharp

#endif  // SHARP_PARAMS_WINTER_H
