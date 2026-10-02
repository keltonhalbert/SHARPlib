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
#include <SHARPlib/qc.h>
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

BourgouinEnergy bourgouin_energy(const float pressure[], const float height[],
                                 const float wetbulb[], const std::ptrdiff_t N,
                                 const float min_energy,
                                 const float pressure_min) {
    // Before the cap below reads pressure[N - 1].
    if (N < 2) return BourgouinEnergy{};

    // Levels [0, cut) are at pressures at or above pressure_min.
    // sharp::upper_bound never returns N, so without the guard it would
    // drop a top level at or above pressure_min.
    const std::ptrdiff_t cut =
        (pressure[N - 1] >= pressure_min)
            ? N
            : upper_bound(pressure, N, pressure_min, std::greater<float>());

    // Areas between Tw and T0 (K m). The lowest cold layer becomes the
    // near-surface cold layer once a warm layer turns up above it.
    float melting_total = 0.0f;
    float melting_above_cold = 0.0f;
    float refreezing_lowest_cold = 0.0f;
    bool found_cold = false;
    bool warm_above_cold = false;
    // The walker reports a lone valid level as one layer with no depth, so
    // the column has depth only with at least 2 valid levels.
    bool has_depth = false;
    for_each_threshold_layer(
        height, [wetbulb](const std::ptrdiff_t k) { return wetbulb[k]; }, cut,
        ZEROCNK, 0.0f, min_energy * ZEROCNK / GRAVITY,
        [&](const HeightLayer& layer, const bool above, const float pos_area,
            const float neg_area) {
            has_depth |= (layer.top > layer.bottom);
            melting_total += pos_area;
            if (found_cold) {
                melting_above_cold += pos_area;
                warm_above_cold |= above;
            } else if (!above) {
                found_cold = true;
                refreezing_lowest_cold = -neg_area;
            }
            return true;
        });

    // Zero energies would read as certain snow.
    if (!has_depth) return BourgouinEnergy{};

    // Eq. 1: energy is g / T0 times the area.
    constexpr float ENERGY_PER_AREA = GRAVITY / ZEROCNK;
    BourgouinEnergy energy;
    energy.melting_energy_total = ENERGY_PER_AREA * melting_total;
    energy.melting_energy_aloft =
        (warm_above_cold) ? ENERGY_PER_AREA * melting_above_cold : 0.0f;
    energy.refreezing_energy =
        (warm_above_cold) ? ENERGY_PER_AREA * refreezing_lowest_cold : 0.0f;
    return energy;
}

// ---------------------------------------------------------------------------
// Precipitation generation layer from a sounding
// ---------------------------------------------------------------------------

HeightLayer precipitation_generation_layer(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const std::ptrdiff_t N, const float min_depth) {
    // Birk et al. (2021) section 3e. Relative humidity as a fraction, depths
    // in meters.
    constexpr float MOIST_RELH = 0.75f;
    constexpr float MIN_GENERATION_DEPTH = 1000.0f;
    constexpr float MAX_DRY_DEPTH = 1500.0f;

    // Over ice below 0 C, over liquid otherwise. MISSING and NaN inputs
    // give MISSING or NaN, which the walker skips.
    const auto relh_at = [&](const std::ptrdiff_t k) {
        return (temperature[k] < ZEROCNK)
                   ? relative_humidity_ice(pressure[k], temperature[k],
                                           dewpoint[k])
                   : relative_humidity(pressure[k], temperature[k],
                                       dewpoint[k]);
    };

    // Layers arrive bottom up, so the last eligible one is the highest.
    HeightLayer generation_layer;
    for_each_threshold_layer(
        height, relh_at, N, MOIST_RELH, min_depth, 0.0f,
        [&](const HeightLayer& layer, const bool moist, float, float) {
            const float depth = layer.top - layer.bottom;
            if (moist) {
                if (depth > MIN_GENERATION_DEPTH) generation_layer = layer;
                return true;
            }
            // Sublimation in a deep dry layer eliminates everything above.
            return depth <= MAX_DRY_DEPTH;
        });
    return generation_layer;
}

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

PrecipTypeProbabilities modified_bourgouin(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float wetbulb[], const std::ptrdiff_t N,
    const float min_depth, const float min_energy, const float pressure_min) {
    // One level has no layers. Return before reading any element.
    if (N < 2) return PrecipTypeProbabilities{};

    const HeightLayer generation_layer = precipitation_generation_layer(
        pressure, height, temperature, dewpoint, N, min_depth);
    // layer_min checks for a MISSING layer only in QC builds.
    if (generation_layer.bottom == MISSING) return PrecipTypeProbabilities{};

    const float prob_ice =
        probability_of_ice(layer_min(generation_layer, height, temperature, N));
    const BourgouinEnergy energy = bourgouin_energy(
        pressure, height, wetbulb, N, min_energy, pressure_min);

    // The surface is the lowest level with a wet-bulb temperature. If there
    // is none, wetbulb[N - 1] is missing too, and every result is MISSING.
    std::ptrdiff_t surface = 0;
#ifndef NO_QC
    while ((surface < N - 1) && is_missing(wetbulb[surface])) {
        ++surface;
    }
#endif
    return modified_bourgouin(energy, prob_ice, wetbulb[surface]);
}

// ===========================================================================
// Precipitation type: the spectral bin classifier (Reeves et al. 2016)
// ===========================================================================

// A port of the spectral bin classifier (SBC) of Reeves, H. D., A. V.
// Ryzhkov, and J. Krause, 2016: Discrimination between winter precipitation
// types based on spectral-bin microphysical modeling. J. Appl. Meteor.
// Climatol., 55, 1747-1761, https://doi.org/10.1175/JAMC-D-16-0044.1
//
// It follows the 2023 version of the algorithm in the Python reference by
// D. Tripp (sbc_alg_2023Aug31.py and run_sbc.py, 2023-08-31), which the
// authors consider authoritative, and in the C++ MRMS code by A. Rosenow and
// D. Tripp (sbcmodel_core.cc and topCalc.cc, versions 2.0.0 to 2.0.3). The
// authors gave permission for this port. The work is NOAA-funded, and the
// port is in the public domain.

// ---------------------------------------------------------------------------
// Result types and the drop-size distribution
// ---------------------------------------------------------------------------

namespace {
// Values of the reference, which differ slightly from SHARPlib's constants.
// Densities in g cm^-3.
constexpr float SBC_ICE_DENSITY = 0.917f;
constexpr float SBC_MAX_SNOW_DENSITY = 0.5f;
// unknwn_e of the reference, sqrt(1 - 0.98^2)
constexpr float SBC_E = 0.198997487f;
// e / asin(e)
constexpr float SBC_E_OVER_ASIN_E = 0.993324402f;
// log((1 + e) / (1 - e))
constexpr float SBC_LOG_E_RATIO = 0.403376976f;

// The reference's capac is this factor times the diameter, plus 0.2 fw.
inline float sbc_capacitance_factor(const float aspect_ratio) {
    return 0.5f * SBC_E_OVER_ASIN_E * 0.8f / std::cbrt(aspect_ratio);
}

// The reference's unknwn_leng is this factor times the diameter.
inline float sbc_length_factor(const float aspect_ratio) {
    return (2.0f + aspect_ratio * aspect_ratio / SBC_E * SBC_LOG_E_RATIO) /
           (4.0f * std::cbrt(aspect_ratio));
}

// Dry-snow density (g cm^-3) at a diameter (mm).
inline float sbc_dry_snow_density(const float rime_factor,
                                  const float diameter) {
    return 0.178f * rime_factor * std::pow(diameter, -0.922f);
}

// uknwn_aa of the reference, from a dry-snow density (g cm^-3).
inline float sbc_fall_speed_aa(const float density) {
    return 1.26f / std::cbrt(density);
}
}  // namespace

SpectralBinDSD spectral_bin_dsd(const float diameter[],
                                const float concentration[],
                                const std::ptrdiff_t nbins,
                                const float rime_factor) {
    // Before any array read.
    if ((nbins < 1) || (nbins > SBC_MAX_BINS)) return SpectralBinDSD{};
    // Fails for NaN too.
    if (!((rime_factor >= 1.0f) && (rime_factor <= 5.0f))) {
        return SpectralBinDSD{};
    }

    // maxD_melt_snow: snow with a smaller melted diameter (mm) has the
    // maximum density.
    const float max_dense_snow_diameter =
        0.154f * std::pow(rime_factor, 1.08f) *
        std::pow(SBC_MAX_SNOW_DENSITY, -0.75f);
    const float large_snow_coeff = 2.29f * std::pow(rime_factor, -0.48f);
    const float dense_snow_coeff = 1.0f / std::cbrt(SBC_MAX_SNOW_DENSITY);
    const float frozen_drop_coeff = std::cbrt(1.0f / SBC_ICE_DENSITY);
    constexpr float MASS_COEFF = PI / 6.0f * 1.0e-3f;

    SpectralBinDSD dsd;
    bool has_positive = false;
    float smaller = 0.0f;
    float rain_fall_speed = 0.01f;
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        const float D = diameter[j];
        const float N = concentration[j];
        if (!std::isfinite(D) || (D <= smaller)) return SpectralBinDSD{};
        if (!std::isfinite(N) || (N < 0.0f)) return SpectralBinDSD{};
        const float D2 = D * D;
        const float D3 = D2 * D;
        const float D4 = D3 * D;
        const float aspect_ratio =
            std::min(0.9951f + 0.02510f * D - 0.03644f * D2 + 0.005303f * D3 -
                         0.0002492f * D4,
                     1.0f);
        if (!(aspect_ratio > 0.0f)) return SpectralBinDSD{};
        has_positive |= (N > 0.0f);
        smaller = D;

        dsd.m_diameter[j] = D;
        dsd.m_concentration[j] = N;
        dsd.m_mass[j] = MASS_COEFF * D3;

        rain_fall_speed =
            std::max(rain_fall_speed, -0.1021f + 4.932f * D - 0.9551f * D2 +
                                          0.07932f * D3 - 0.002362f * D4);
        dsd.m_rain_fall_speed[j] = rain_fall_speed;

        const float Di = frozen_drop_coeff * D;
        dsd.m_pellet_fall_speed[j] = 0.2259f + 1.5954f * Di - 0.0405f * Di * Di;
        dsd.m_foote_du_toit_fall_speed[j] =
            -0.193f + 4.96f * D - 0.904f * D2 + 0.0566f * D3;

        dsd.m_liquid_capacitance_factor[j] =
            sbc_capacitance_factor(aspect_ratio);
        dsd.m_liquid_length_factor[j] = sbc_length_factor(aspect_ratio);

        const float snow_diameter =
            (D >= max_dense_snow_diameter)
                ? large_snow_coeff * std::pow(D, 1.443f)
                : dense_snow_coeff * D;
        const float snow_density =
            sbc_dry_snow_density(rime_factor, snow_diameter);
        dsd.m_snow_diameter[j] = snow_diameter;
        dsd.m_snow_density[j] = std::min(snow_density, SBC_MAX_SNOW_DENSITY);
        dsd.m_snow_aa[j] = sbc_fall_speed_aa(snow_density);
        dsd.m_snow_bb[j] = (dsd.m_snow_aa[j] - 1.0f) / 2.0f;

        dsd.m_liquid_aa[j] =
            sbc_fall_speed_aa(sbc_dry_snow_density(rime_factor, D));
        dsd.m_liquid_bb[j] = (dsd.m_liquid_aa[j] - 1.0f) / 2.0f;
    }
    if (!has_positive) return SpectralBinDSD{};

    dsd.m_nbins = nbins;
    dsd.m_rime_factor = rime_factor;
    return dsd;
}

SpectralBinDSD spectral_bin_dsd_default() {
    constexpr std::ptrdiff_t NBINS = 4;
    constexpr float DIAMETER[NBINS] = {0.05f, 0.75f, 1.45f, 2.15f};
    constexpr float CONCENTRATION[NBINS] = {55.1843f, 146.647f, 11.6891f,
                                            3.60886f};
    return spectral_bin_dsd(DIAMETER, CONCENTRATION, NBINS);
}

// ---------------------------------------------------------------------------
// Cloud top from a sounding
// ---------------------------------------------------------------------------

namespace {
// Thresholds of the reference's cloud-top rule. Dewpoint depressions in K,
// relative humidities as fractions.
constexpr float SBC_CLOUD_DEPRESSION = 6.0f;
constexpr float SBC_CLOUD_RELH = 0.60f;
constexpr float SBC_DRY_DEPRESSION = 10.0f;
constexpr float SBC_DRY_RELH = 0.40f;
constexpr float SBC_FALLBACK_RELH = 0.80f;

inline bool sbc_is_cloud(const float depression, const float relh) {
    return (depression <= SBC_CLOUD_DEPRESSION) && (relh > SBC_CLOUD_RELH);
}

// Level of the cloud top, or -1 for no cloud. One pass from the top: the
// first loop finds the highest cloud level, and the second tracks the dry
// test and the highest cloud level at or below the driest level so far.
std::ptrdiff_t sbc_cloud_top_level(const float temperature[],
                                   const float dewpoint[], const float relh[],
                                   const std::ptrdiff_t N) {
    std::ptrdiff_t k = N - 1;
    std::ptrdiff_t fallback = -1;
    bool below_highest = false;
    for (; k >= 0; --k) {
#ifndef NO_QC
        if (is_missing(temperature[k]) || is_missing(dewpoint[k]) ||
            is_missing(relh[k])) {
            continue;
        }
#endif
        if (sbc_is_cloud(temperature[k] - dewpoint[k], relh[k])) break;
        if ((fallback < 0) && below_highest &&
            (relh[k] >= SBC_FALLBACK_RELH)) {
            fallback = k;
        }
        below_highest = true;
    }
    if (k < 0) return fallback;

    const std::ptrdiff_t first_top = k;
    float max_depression = temperature[k] - dewpoint[k];
    std::ptrdiff_t below_driest = k;
    bool dry = false;
    for (--k; k >= 0; --k) {
#ifndef NO_QC
        if (is_missing(temperature[k]) || is_missing(dewpoint[k]) ||
            is_missing(relh[k])) {
            continue;
        }
#endif
        const float depression = temperature[k] - dewpoint[k];
        dry |= (depression > SBC_DRY_DEPRESSION) || (relh[k] < SBC_DRY_RELH);
        if (depression > max_depression) {
            max_depression = depression;
            below_driest = -1;
        }
        if ((below_driest < 0) && sbc_is_cloud(depression, relh[k])) {
            below_driest = k;
        }
    }
    return (dry && (below_driest >= 0)) ? below_driest : first_top;
}
}  // namespace

float spectral_bin_cloud_top([[maybe_unused]] const float pressure[],
                             const float height[], const float temperature[],
                             const float dewpoint[], const float relh[],
                             const std::ptrdiff_t N) {
    if (N < 1) return MISSING;
    const std::ptrdiff_t top =
        sbc_cloud_top_level(temperature, dewpoint, relh, N);
    return (top < 0) ? MISSING : height[top];
}

// ---------------------------------------------------------------------------
// Precipitation type from a given cloud top: pre-classifier
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Microphysics: frozen cloud tops and melting
// ---------------------------------------------------------------------------

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
