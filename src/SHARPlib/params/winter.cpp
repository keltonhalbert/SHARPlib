/**
 * \file
 * \brief Winter weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Prediction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/params/winter.h>
#include <SHARPlib/thermo.h>

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sharp {

PressureLayer dendritic_layer(const float pressure[], const float temperature[],
                              const std::ptrdiff_t N) {
    return temperature_layer(pressure, temperature, -17.0f + ZEROCNK,
                             -12.0f + ZEROCNK, N);
}

float snow_squall_parameter(const float wetbulb_2m, const float mean_relh_0_2km,
                            const float delta_thetae_0_2km,
                            const float mean_wind_0_2km) {
    if (wetbulb_2m > ZEROCNK + 1) return 0.0f;
    const float relh_term = (mean_relh_0_2km - 0.60f) / 0.15f;
    const float thetae_term = (4.0f - delta_thetae_0_2km) / 4.0f;
    if ((relh_term < 0) || (thetae_term < 0)) return 0.0f;

    return relh_term * thetae_term * (mean_wind_0_2km / 9.0f);
}

// ===========================================================================
// Precipitation type: the modified Bourgouin method (Birk et al. 2021)
// ===========================================================================

// ---------------------------------------------------------------------------
// Wet-bulb melting and refreezing energies from a sounding
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Precipitation generation layer from a sounding
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Probability of ice, and precipitation-type probabilities from energies
// ---------------------------------------------------------------------------

float probability_of_ice(const float temperature) {
    // Fails for both MISSING and NaN.
    if (!(temperature > 0.0f)) return MISSING;

    // Birk et al. (2021) Eq. 2, from Baumgardt et al. (2017). T in Celsius,
    // result in percent.
    constexpr float ALWAYS_ICE_TEMP = -15.0f;
    constexpr float NEVER_ICE_TEMP = -7.0f;
    constexpr float C4 = -0.065f;
    constexpr float C3 = -3.1544f;
    constexpr float C2 = -56.414f;
    constexpr float C1 = -449.6f;
    constexpr float C0 = -1308.0f;

    const float tmpc = temperature - ZEROCNK;
    if (tmpc <= ALWAYS_ICE_TEMP) return 1.0f;
    if (tmpc >= NEVER_ICE_TEMP) return 0.0f;
    const float pct = (((C4 * tmpc + C3) * tmpc + C2) * tmpc + C1) * tmpc + C0;
    return std::clamp(pct, 0.0f, 100.0f) / 100.0f;
}

PrecipTypeProbabilities modified_bourgouin(const BourgouinEnergy& energy,
                                           const float prob_ice,
                                           const float surface_wetbulb) {
    const float me_total = energy.melting_energy_total;
    const float me_aloft = energy.melting_energy_aloft;
    const float re = energy.refreezing_energy;

    // Each comparison fails for both MISSING and NaN.
    if (!(prob_ice >= 0.0f) || !(me_total >= 0.0f) || !(me_aloft >= 0.0f) ||
        !(re >= 0.0f) || !(surface_wetbulb > 0.0f)) {
        return PrecipTypeProbabilities{};
    }

    // Birk et al. (2021) Appendix. Energies in J/kg, results in percent.
    // Snow, Eq. 9
    constexpr float SN_SCALE = 1540.0f;
    constexpr float SN_ME = -0.29f;
    // Ice pellets, Eq. 8
    constexpr float PL_RE = 2.3f;
    constexpr float PL_ME = -42.0f;
    constexpr float PL_OFFSET = 3.0f;
    // Freezing rain or rain, Eq. 7 and the taper for weak melting
    constexpr float FZ_RE = -2.1f;
    constexpr float FZ_ME = 0.2f;
    constexpr float FZ_OFFSET = 458.0f;
    constexpr float FZ_TAPER_BELOW_ME = 5.0f;
    constexpr float FZ_TAPER_PER_ME = 0.2f;

    const auto clamp_pct = [](const float pct) {
        return std::clamp(pct, 0.0f, 100.0f);
    };

    const float snow_i = clamp_pct(SN_SCALE * std::exp(SN_ME * me_total));
    const float snow = clamp_pct(prob_ice * snow_i);

    // The paper computes ice pellets only for melting aloft above a
    // refreezing layer, so the probability jumps as me_aloft leaves 0.
    const float pellets_i =
        (me_aloft > 0.0f && re > 0.0f)
            ? clamp_pct(PL_RE * re + PL_ME * std::log(me_aloft + 1.0f) +
                        PL_OFFSET)
            : 0.0f;
    const float pellets = clamp_pct(prob_ice * pellets_i);

    // Clamp first, then taper (Appendix c, steps 1 and 2).
    float liquid_i = clamp_pct(FZ_RE * re + FZ_ME * me_total + FZ_OFFSET);
    if (me_total < FZ_TAPER_BELOW_ME) liquid_i *= FZ_TAPER_PER_ME * me_total;
    const float liquid =
        clamp_pct(100.0f * (1.0f - prob_ice) + prob_ice * liquid_i);

    const bool warm_surface = surface_wetbulb > ZEROCNK;
    PrecipTypeProbabilities probs;
    probs.rain = warm_surface ? liquid / 100.0f : 0.0f;
    probs.snow = snow / 100.0f;
    probs.freezing_rain = warm_surface ? 0.0f : liquid / 100.0f;
    probs.ice_pellets = pellets / 100.0f;
    return probs;
}

// ---------------------------------------------------------------------------
// Precipitation type from a full sounding
// ---------------------------------------------------------------------------

}  // namespace sharp
