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

}  // namespace sharp

#endif  // SHARP_PARAMS_WINTER_H
