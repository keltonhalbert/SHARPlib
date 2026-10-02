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
    // Every return returns dsd, so that it is built in place. A rejected
    // distribution is a default-constructed one.
    SpectralBinDSD dsd;
    // Before any array read.
    if ((nbins < 1) || (nbins > SBC_MAX_BINS)) return dsd;
    // Fails for NaN too.
    if (!((rime_factor >= 1.0f) && (rime_factor <= 5.0f))) return dsd;

    // pow(0.5, -0.75), with the maximum snow density 0.5
    constexpr float MAX_DENSITY_POW = 1.68179286f;
    // 1 / cbrt(0.5)
    constexpr float DENSE_SNOW_COEFF = 1.25992107f;
    // cbrt(1 / 0.917), with the ice density 0.917
    constexpr float FROZEN_DROP_COEFF = 1.02930379f;
    // maxD_melt_snow: snow with a smaller melted diameter (mm) has the
    // maximum density.
    const float max_dense_snow_diameter =
        0.154f * std::pow(rime_factor, 1.08f) * MAX_DENSITY_POW;
    const float large_snow_coeff = 2.29f * std::pow(rime_factor, -0.48f);
    constexpr float MASS_COEFF = PI / 6.0f * 1.0e-3f;

    bool has_positive = false;
    float smaller = 0.0f;
    float rain_fall_speed = 0.01f;
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        const float D = diameter[j];
        const float N = concentration[j];
        if (!std::isfinite(D) || (D <= smaller) || !std::isfinite(N) ||
            (N < 0.0f)) {
            dsd = SpectralBinDSD{};
            return dsd;
        }
        const float D2 = D * D;
        const float D3 = D2 * D;
        const float D4 = D3 * D;
        const float aspect_ratio =
            std::min(0.9951f + 0.02510f * D - 0.03644f * D2 + 0.005303f * D3 -
                         0.0002492f * D4,
                     1.0f);
        if (!(aspect_ratio > 0.0f)) {
            dsd = SpectralBinDSD{};
            return dsd;
        }
        has_positive |= (N > 0.0f);
        smaller = D;

        dsd.m_diameter[j] = D;
        dsd.m_concentration[j] = N;
        dsd.m_mass[j] = MASS_COEFF * D3;

        rain_fall_speed =
            std::max(rain_fall_speed, -0.1021f + 4.932f * D - 0.9551f * D2 +
                                          0.07932f * D3 - 0.002362f * D4);
        dsd.m_rain_fall_speed[j] = rain_fall_speed;

        const float Di = FROZEN_DROP_COEFF * D;
        dsd.m_pellet_fall_speed[j] = 0.2259f + 1.5954f * Di - 0.0405f * Di * Di;
        dsd.m_foote_du_toit_fall_speed[j] =
            -0.193f + 4.96f * D - 0.904f * D2 + 0.0566f * D3;

        dsd.m_liquid_capacitance_factor[j] =
            sbc_capacitance_factor(aspect_ratio);
        dsd.m_liquid_length_factor[j] = sbc_length_factor(aspect_ratio);

        const float snow_diameter =
            (D >= max_dense_snow_diameter)
                ? large_snow_coeff * std::pow(D, 1.443f)
                : DENSE_SNOW_COEFF * D;
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
    if (!has_positive) {
        dsd = SpectralBinDSD{};
        return dsd;
    }

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
// test and the highest cloud level at or below the driest level so far. The
// missing-data tests combine with | rather than ||, for one branch per level.
std::ptrdiff_t sbc_cloud_top_level(const float temperature[],
                                   const float dewpoint[], const float relh[],
                                   const std::ptrdiff_t N) {
    std::ptrdiff_t k = N - 1;
    std::ptrdiff_t fallback = -1;
    bool below_highest = false;
    for (; k >= 0; --k) {
#ifndef NO_QC
        if (is_missing(temperature[k]) | is_missing(dewpoint[k]) |
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
        if (is_missing(temperature[k]) | is_missing(dewpoint[k]) |
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

namespace {
// The levels the classifier reads, from the surface (the lowest valid level)
// up to the top (the highest valid level at or below the cloud top). Levels
// between them that are not valid are skipped.
struct SpectralBinColumn {
    const float* pressure;
    const float* height;
    const float* temperature;
    const float* dewpoint;
    const float* relh;
    const float* wetbulb;
    std::ptrdiff_t surface;
    std::ptrdiff_t top;

    [[nodiscard]] bool valid([[maybe_unused]] const std::ptrdiff_t k) const {
#ifdef NO_QC
        return true;
#else
        // | rather than ||, so that a level costs one branch, not four
        return !(is_missing(temperature[k]) | is_missing(dewpoint[k]) |
                 is_missing(relh[k]) | is_missing(wetbulb[k]));
#endif
    }

    // Height above the surface (m)
    [[nodiscard]] float height_agl(const std::ptrdiff_t k) const {
        return height[k] - height[surface];
    }
};

constexpr SpectralBinResult SBC_SNOW{PrecipType::snow, 0.0f, MISSING};
constexpr SpectralBinResult SBC_FREEZING_RAIN{PrecipType::freezing_rain, 1.0f,
                                              0.0f};
constexpr SpectralBinResult SBC_RAIN{PrecipType::rain, 1.0f, MISSING};

SpectralBinResult sbc_microphysics(const SpectralBinColumn& column,
                                   const SpectralBinDSD& dsd,
                                   float ice_nucleation_temperature,
                                   float liquid_fraction_profile[]);
}  // namespace

SpectralBinResult spectral_bin_classifier(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float relh[], const float wetbulb[],
    const std::ptrdiff_t N, const float cloud_top, const SpectralBinDSD& dsd,
    const float ice_nucleation_temperature, float liquid_fraction_profile[]) {
    const std::ptrdiff_t nbins = dsd.nbins();
    if ((liquid_fraction_profile != nullptr) && (N > 0)) {
        std::fill_n(liquid_fraction_profile, N * nbins, MISSING);
    }
    // Before any array read.
    if (N < 2) return SpectralBinResult{};
    // Fails for NaN too.
    if ((nbins == 0) || !(ice_nucleation_temperature > 0.0f) ||
        is_missing(cloud_top)) {
        return SpectralBinResult{};
    }

    // sharp::upper_bound never returns N, so a cloud top at or above the
    // highest level needs the guard.
    const std::ptrdiff_t top = (height[N - 1] <= cloud_top)
                                   ? N - 1
                                   : upper_bound(height, N, cloud_top) - 1;
    SpectralBinColumn column{pressure, height, temperature, dewpoint,
                             relh,     wetbulb, 0,          top};
    while ((column.surface < column.top) && !column.valid(column.surface)) {
        ++column.surface;
    }
    while ((column.top > column.surface) && !column.valid(column.top)) {
        --column.top;
    }
    if (column.top <= column.surface) return SpectralBinResult{};

    const float tice = ice_nucleation_temperature;
    const float tw_top = wetbulb[column.top];
    const float tw_surface = wetbulb[column.surface];

    // Rule 2 needs no pass over the column, and where it holds, rule 1 gives
    // FZRA too.
    if ((tw_top > tice) && (tw_surface < ZEROCNK)) return SBC_FREEZING_RAIN;

    // Bottom up, so of two levels equally near 3 km, the higher one wins, as
    // the reference takes the first from the top. The surface is 3 km away.
    constexpr float AGL_3KM = 3000.0f;
    float tw_min = tw_surface;
    float tw_max = tw_surface;
    float nearest_3km = AGL_3KM;
    float tw_min_below_3km = tw_surface;
    for (std::ptrdiff_t k = column.surface + 1; k <= column.top; ++k) {
        if (!column.valid(k)) continue;
        tw_min = std::min(tw_min, wetbulb[k]);
        tw_max = std::max(tw_max, wetbulb[k]);
        const float distance = std::abs(column.height_agl(k) - AGL_3KM);
        if (distance <= nearest_3km) {
            nearest_3km = distance;
            tw_min_below_3km = tw_min;
        }
    }

    if (tw_max <= ZEROCNK) {
        return ((tw_top < tice) && (tw_min_below_3km < tice))
                   ? SBC_SNOW
                   : SBC_FREEZING_RAIN;
    }
    if (tw_min > ZEROCNK) return SBC_RAIN;

    return sbc_microphysics(column, dsd, ice_nucleation_temperature,
                            liquid_fraction_profile);
}

// ---------------------------------------------------------------------------
// Microphysics: frozen cloud tops and melting
// ---------------------------------------------------------------------------

namespace {
// Values of the reference, in its units. Densities in g cm^-3.
constexpr float SBC_SEA_LEVEL_AIR_DENSITY = 1.292e-3f;
constexpr float SBC_DRY_AIR_GAS_CONSTANT = 287.0f;
// The per-bin class uses this ice density, and the surface SBC_ICE_DENSITY.
constexpr float SBC_CLASS_ICE_DENSITY = 0.918f;
constexpr float SBC_CLASS_FRACTION = 0.15f;
constexpr float SBC_GRAUPEL_RIME_FACTOR = 5.0f;
// 1.81e-5 * 5^3.26 (cm^3)
constexpr float SBC_GRAUPEL_SNOW_VOLUME_MIN = 0.00343811815f;
// 1e-3 / 1.8e-5, mm to m over the kinematic viscosity of air (s m^-2)
constexpr float SBC_REYNOLDS_COEFF = 55.5555573f;
// (1.8e-5 / 2.0e-5)^(1/3), the cube root of the Schmidt number
constexpr float SBC_SCHMIDT_CBRT = 0.965489388f;
constexpr float SBC_AIR_CONDUCTIVITY = 0.023f;
// 2.0e-5 * 2.5e6 / 461.5, vapor diffusivity times the latent heat of
// vaporization over the gas constant of water vapor
constexpr float SBC_VAPOR_HEAT_COEFF = 0.108342364f;
// 611 / 273.15
constexpr float SBC_SVP_OVER_T0 = 2.23686624f;
constexpr float SBC_HEAT_CLAMP_TW = 273.65f;
constexpr float SBC_MIN_HEAT = 0.02f;
// 4 pi / 3.35e5, over the latent heat of melting
constexpr float SBC_MELT_COEFF = 3.75115522e-05f;
// 6 / pi
constexpr float SBC_SIX_OVER_PI = 1.90985930f;

// sbc_capacitance_factor(0.8) and sbc_length_factor(0.8), for the snow aspect
// ratio of 0.8
constexpr float SBC_SNOW_CAPACITANCE_FACTOR = 0.428010523f;
constexpr float SBC_SNOW_LENGTH_FACTOR = 0.887979627f;

// t - 273.15 (C): t - 273.15f is exact, and 6.103515625e-6 = 273.15 - 273.15f
inline float sbc_celsius(const float t) {
    return (t - ZEROCNK) - 6.103515625e-6f;
}

// psd_ptype of the reference. Its codes are these values plus 0.5, and -999
// for unset.
enum class SBCClass : unsigned char {
    unset = 0,
    rain = 1,                          // RA
    freezing_drizzle = 2,              // FZDZ
    freezing_rain = 3,                 // FZRA
    ice_pellets = 4,                   // PL
    freezing_drizzle_ice_pellets = 5,  // FZDZPL
    freezing_rain_ice_pellets = 6,     // FZRAPL
    snow = 7,                          // SN
    rain_snow = 8,                     // RASN
    rain_ice_pellets = 9,              // RAPL
};

// liquid_ar of the reference: the classes that melt with the raindrop
// aspect ratio. The others melt with the snow aspect ratio of 0.8.
inline bool sbc_liquid_ar(const SBCClass ptype) {
    switch (ptype) {
        case SBCClass::rain:
        case SBCClass::freezing_drizzle:
        case SBCClass::freezing_rain:
        case SBCClass::ice_pellets:
        case SBCClass::freezing_drizzle_ice_pellets:
        case SBCClass::freezing_rain_ice_pellets:
        case SBCClass::rain_ice_pellets:
            return true;
        case SBCClass::unset:
        case SBCClass::snow:
        case SBCClass::rain_snow:
            return false;
    }
    return false;
}

// The class of a bin with supercooled liquid, alone or with ice pellets:
// drizzle below 0.6 mm (diameter in mm), rain otherwise
inline SBCClass sbc_freezing_class(const float diameter,
                                   const bool ice_pellets) {
    if (diameter < 0.6f) {
        return ice_pellets ? SBCClass::freezing_drizzle_ice_pellets
                           : SBCClass::freezing_drizzle;
    }
    return ice_pellets ? SBCClass::freezing_rain_ice_pellets
                       : SBCClass::freezing_rain;
}

// Every field of every bin at one level, named as in the reference.
struct SBCLevelState {
    std::array<float, SBC_MAX_BINS> water_fraction;      // fw
    std::array<float, SBC_MAX_BINS> velocity_melt_snow;  // v (m/s)
    std::array<float, SBC_MAX_BINS> mass_water;          // mw (g)
    std::array<float, SBC_MAX_BINS> mass_snow;           // ms (g)
    std::array<float, SBC_MAX_BINS> volume_liq;          // vl (cm^3)
    std::array<float, SBC_MAX_BINS> volume_ice;          // vi (cm^3)
    std::array<float, SBC_MAX_BINS> volume_snow;         // vs (cm^3)
    std::array<float, SBC_MAX_BINS> diam_melt_snow;      // dm (mm)
    // in_diam_sn_core (dc, mm) is sn_core_scale cbrt(sn_core_cube), or 0
    // where sn_core_scale is 0, at a level with no snow core. Only a refreeze
    // level reads it, so the cube root waits until then.
    std::array<float, SBC_MAX_BINS> sn_core_cube;
    std::array<SBCClass, SBC_MAX_BINS> psd_ptype;
    float sn_core_scale;

    // The reference's arrays start at 0 (psd_ptype at unset), so a field
    // that a branch does not write reads back as 0 at the next level. The
    // core clears a level before a branch writes it, except before branch
    // B, which writes the zeros itself.
    void clear(const std::ptrdiff_t nbins) {
        for (std::ptrdiff_t j = 0; j < nbins; ++j) {
            water_fraction[j] = 0.0f;
            velocity_melt_snow[j] = 0.0f;
            mass_water[j] = 0.0f;
            mass_snow[j] = 0.0f;
            volume_liq[j] = 0.0f;
            volume_ice[j] = 0.0f;
            volume_snow[j] = 0.0f;
            diam_melt_snow[j] = 0.0f;
            sn_core_cube[j] = 0.0f;
            psd_ptype[j] = SBCClass::unset;
        }
        sn_core_scale = 0.0f;
    }
};

// One level of the integration, which runs from the cloud top down.
struct SBCLevel {
    // Index in the input arrays
    std::ptrdiff_t k;
    // The valid level above, the reference's K - 1; -1 at the cloud top
    std::ptrdiff_t k_above;
    // Valid levels below the cloud top, the reference's K
    std::ptrdiff_t K;
    // (K)
    float wetbulb;
    // density_air_g, from the wet-bulb temperature (g cm^-3)
    float air_density;
    // layer_depth: height[k_above] - height[k], 0 at the cloud top (m)
    float layer_depth;

    // Scales the fall speeds near the ground to this level.
    [[nodiscard]] float density_correction() const {
        return std::sqrt(SBC_SEA_LEVEL_AIR_DENSITY / air_density);
    }
};

// What a column carries from level to level, besides the records of the
// levels K - 1 and K.
struct SBCColumnState {
    // Level 0, the cloud top
    SBCLevelState top;
    // uknwn_aa and uknwn_bb, set by the cloud-top branch
    const float* uknwn_aa = nullptr;
    const float* uknwn_bb = nullptr;

    // The pre-pass, with the original Tice. A level index of the number of
    // levels in the column means none, like the reference's len + 1000.
    std::ptrdiff_t crossings = 0;
    // cross_index[0]: the first 0 C crossing
    std::ptrdiff_t cross_index = 0;
    // ice_nuc_tops[0]: the first level that warms through Tice
    std::ptrdiff_t ice_nuc_top = 0;

    // temp_nuc: the Tice in use (K)
    float temp_nuc = MISSING;
    float rime_factor = MISSING;
    // 1.81e-5 * rime_factor^3.26 (cm^3)
    float snow_volume_min = MISSING;
    bool refrz_detected_flag = false;
    std::ptrdiff_t refrz_lvl = 0;
    std::array<bool, SBC_MAX_BINS> refrz_lvl_flag;
    // in_diam_sn_core, water_fraction, and velocity_melt_snow at refrz_lvl
    std::array<float, SBC_MAX_BINS> refrz_in_diam_sn_core;
    std::array<float, SBC_MAX_BINS> refrz_water_fraction;
    std::array<float, SBC_MAX_BINS> refrz_velocity_melt_snow;

    // slw_hgt (m AGL)
    float slw_hgt = MISSING;

    // The refreeze level becomes K, and level holds its values.
    void set_refrz_lvl(const std::ptrdiff_t K, const SBCLevelState& level,
                       const std::ptrdiff_t nbins) {
        refrz_lvl = K;
        if (level.sn_core_scale == 0.0f) {
            std::fill_n(refrz_in_diam_sn_core.begin(), nbins, 0.0f);
        } else {
            for (std::ptrdiff_t j = 0; j < nbins; ++j) {
                refrz_in_diam_sn_core[j] =
                    level.sn_core_scale * std::cbrt(level.sn_core_cube[j]);
            }
        }
        std::copy_n(level.water_fraction.begin(), nbins,
                    refrz_water_fraction.begin());
        std::copy_n(level.velocity_melt_snow.begin(), nbins,
                    refrz_velocity_melt_snow.begin());
    }
};

// The branches of classify, in its order. A cloud-top branch writes level 0
// to next. The others write level K to next from level K - 1 in prev. next
// arrives cleared, so a branch writes only the fields the reference writes.
// Branch B is the exception: it writes every field.

// Branch A (classify:138): a frozen cloud top
inline void sbc_frozen_cloud_top(
    SBCColumnState& state, SBCLevelState& next, const SBCLevel& level,
    [[maybe_unused]] const SpectralBinColumn& column,
    const SpectralBinDSD& dsd) {
    state.uknwn_aa = dsd.snow_aa().data();
    state.uknwn_bb = dsd.snow_bb().data();
    const float rho_air = level.air_density;
    const float ice_air = SBC_ICE_DENSITY - rho_air;
    const float density_correction = level.density_correction();
    const float* m0 = dsd.mass().data();
    const float* diameter = dsd.snow_diameter().data();
    const float* rho_snow = dsd.snow_density().data();
    const float* v_rain = dsd.rain_fall_speed().data();
    for (std::ptrdiff_t j = 0; j < dsd.nbins(); ++j) {
        next.diam_melt_snow[j] = diameter[j];
        next.mass_snow[j] = m0[j];
        next.volume_ice[j] =
            m0[j] / rho_snow[j] * (rho_snow[j] - rho_air) / ice_air;
        next.volume_snow[j] = m0[j] / rho_snow[j];
        next.velocity_melt_snow[j] =
            v_rain[j] * density_correction / state.uknwn_aa[j];
    }
}

// Branch B (classify:164): the snow of a frozen cloud top, above the first
// crossing. Most levels of such a column take it, so it writes the zeros of
// a cleared level itself, which saves a pass over the level.
inline void sbc_frozen_above_crossing(
    SBCColumnState& state, [[maybe_unused]] const SBCLevelState& prev,
    SBCLevelState& next, const SBCLevel& level,
    [[maybe_unused]] const SpectralBinColumn& column,
    const SpectralBinDSD& dsd) {
    const SBCLevelState& top = state.top;
    const float density_correction = level.density_correction();
    const float* v_rain = dsd.rain_fall_speed().data();
    for (std::ptrdiff_t j = 0; j < dsd.nbins(); ++j) {
        next.water_fraction[j] = 0.0f;
        next.mass_water[j] = 0.0f;
        next.volume_liq[j] = 0.0f;
        next.sn_core_cube[j] = 0.0f;
        next.mass_snow[j] = top.mass_snow[j];
        next.volume_ice[j] = top.volume_ice[j];
        next.volume_snow[j] = top.volume_snow[j];
        next.diam_melt_snow[j] = top.diam_melt_snow[j];
        next.velocity_melt_snow[j] =
            v_rain[j] * density_correction / state.uknwn_aa[j];
        next.psd_ptype[j] = SBCClass::snow;
    }
    next.sn_core_scale = 0.0f;
}

// Branches C (classify:186), D (:217), and E (:243), defined under
// "Microphysics: liquid cloud tops"
inline void sbc_liquid_cloud_top(SBCColumnState& state, SBCLevelState& next,
                                 const SBCLevel& level,
                                 const SpectralBinColumn& column,
                                 const SpectralBinDSD& dsd);
inline void sbc_supercooled_above_crossing(SBCColumnState& state,
                                           const SBCLevelState& prev,
                                           SBCLevelState& next,
                                           const SBCLevel& level,
                                           const SpectralBinColumn& column,
                                           const SpectralBinDSD& dsd);
inline void sbc_warm_above_crossing(SBCColumnState& state,
                                    const SBCLevelState& prev,
                                    SBCLevelState& next, const SBCLevel& level,
                                    const SpectralBinColumn& column,
                                    const SpectralBinDSD& dsd);

// Branch F (classify:261): a melting level
inline void sbc_melting(SBCColumnState& state, const SBCLevelState& prev,
                        SBCLevelState& next, const SBCLevel& level,
                        const SpectralBinColumn& column,
                        const SpectralBinDSD& dsd) {
    const std::ptrdiff_t nbins = dsd.nbins();
    if (state.refrz_detected_flag) {
        state.rime_factor = SBC_GRAUPEL_RIME_FACTOR;
        state.snow_volume_min = SBC_GRAUPEL_SNOW_VOLUME_MIN;
    }
    std::fill_n(state.refrz_lvl_flag.begin(), nbins, false);

    bool melted = true;
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        melted &= (prev.water_fraction[j] == 1.0f);
    }
    if (melted) {
        for (std::ptrdiff_t j = 0; j < nbins; ++j) {
            next.water_fraction[j] = 1.0f;
            next.velocity_melt_snow[j] = prev.velocity_melt_snow[j];
            next.mass_water[j] = prev.mass_water[j];
            next.volume_liq[j] = prev.volume_liq[j];
            next.diam_melt_snow[j] = prev.diam_melt_snow[j];
            next.psd_ptype[j] = SBCClass::rain;
        }
        return;
    }

    // cm to mm
    next.sn_core_scale = 10.0f;
    const float tw = level.wetbulb;
    const float tw_c = sbc_celsius(tw);
    const float rho_air = level.air_density;
    const float density_correction = level.density_correction();
    // svp_wrt_water (Pa)
    const float svp = 611.0f * std::exp(17.269f * tw_c / (tw - 35.86f));
    float heat = SBC_AIR_CONDUCTIVITY * tw_c +
                 SBC_VAPOR_HEAT_COEFF *
                     (column.relh[level.k] * svp / tw - SBC_SVP_OVER_T0);
    if ((heat < 0.0f) && (tw < SBC_HEAT_CLAMP_TW)) heat = 0.0f;
    if ((heat < SBC_MIN_HEAT) && (tw > SBC_HEAT_CLAMP_TW)) heat = SBC_MIN_HEAT;
    const float melt_coeff = SBC_MELT_COEFF * heat * level.layer_depth;

    const float ice_air = SBC_ICE_DENSITY - rho_air;
    const float dense_snow_ratio = ice_air / (SBC_MAX_SNOW_DENSITY - rho_air);
    const float rime = state.rime_factor;
    const float snow_volume_min = state.snow_volume_min;
    const float snow_density_coeff = 1.75e-2f * rime;
    const bool refrozen = (state.refrz_lvl != 0);
    const SBCClass ice_class =
        refrozen ? SBCClass::ice_pellets : SBCClass::snow;
    const SBCClass mix_class =
        refrozen ? SBCClass::rain_ice_pellets : SBCClass::rain_snow;

    const float* m0 = dsd.mass().data();
    const float* N = dsd.concentration().data();
    const float* v_rain = dsd.rain_fall_speed().data();
    const float* v_pellet = dsd.pellet_fall_speed().data();
    const float* liquid_capacitance = dsd.liquid_capacitance_factor().data();
    const float* liquid_length = dsd.liquid_length_factor().data();
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        const float fw_prev = prev.water_fraction[j];
        if (fw_prev == 1.0f) {
            next.water_fraction[j] = 1.0f;
            next.velocity_melt_snow[j] = prev.velocity_melt_snow[j];
            next.mass_water[j] = prev.mass_water[j];
            next.volume_liq[j] = prev.volume_liq[j];
            next.diam_melt_snow[j] = prev.diam_melt_snow[j];
            next.psd_ptype[j] = SBCClass::rain;
            continue;
        }

        const bool liquid_ar = sbc_liquid_ar(prev.psd_ptype[j]);
        const float dm_prev = prev.diam_melt_snow[j];
        const float capac =
            (liquid_ar ? liquid_capacitance[j] : SBC_SNOW_CAPACITANCE_FACTOR) *
                dm_prev +
            0.2f * fw_prev;
        const float length =
            (liquid_ar ? liquid_length[j] : SBC_SNOW_LENGTH_FACTOR) * dm_prev;
        const float v =
            state.refrz_detected_flag
                ? v_pellet[j] * density_correction
                : v_rain[j] * density_correction /
                      (state.uknwn_aa[j] -
                       state.uknwn_bb[j] * fw_prev * (1.0f + fw_prev));
        const float melting_schmidt =
            SBC_SCHMIDT_CBRT * std::sqrt(length * SBC_REYNOLDS_COEFF * v);
        const float vent =
            (melting_schmidt <= 1.0f)
                ? 1.0f + 0.14f * melting_schmidt * melting_schmidt
                : 0.86f + 0.28f * melting_schmidt;

        const float ms_prev = prev.mass_snow[j];
        float fw = fw_prev;
        if (ms_prev == 0.0f) {
            fw = 1.0f;
        } else if (heat > 0.0f) {
            fw = std::min(fw_prev + vent * capac / (ms_prev * v) * melt_coeff,
                          1.0f);
        }
        fw = std::max(fw, 0.0f);

        const float change_ice = ms_prev * fw / SBC_ICE_DENSITY;
        const float vi = std::max(state.top.volume_ice[j] - change_ice, 0.0f);
        const float rhs = ice_air * vi + rho_air * prev.volume_snow[j];
        const float vl = change_ice * SBC_ICE_DENSITY;

        next.water_fraction[j] = fw;
        next.velocity_melt_snow[j] = v;
        next.mass_water[j] = m0[j] * fw;
        next.volume_liq[j] = vl;
        next.volume_ice[j] = vi;
        // The base of the power law for vs, which the next loop replaces
        next.volume_snow[j] = rhs / rime;
    }
    // At the surface, the surface decision and the profile read only fw and
    // v, so the snow and the diameters are left undone.
    if (level.k == column.surface) return;
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        if (prev.water_fraction[j] == 1.0f) continue;
        const float vi = next.volume_ice[j];
        const float vl = next.volume_liq[j];
        float vs = 343.0f * std::pow(next.volume_snow[j], 1.443f);
        if (vs < snow_volume_min) vs = dense_snow_ratio * vi;
        const float snow_density =
            (vs >= snow_volume_min)
                ? snow_density_coeff * std::pow(vs, -0.307f)
                : SBC_MAX_SNOW_DENSITY;
        next.mass_snow[j] = snow_density * vs + vl;
        next.volume_snow[j] = vs;
        // Its cube root follows, in a loop of its own.
        next.diam_melt_snow[j] = SBC_SIX_OVER_PI * (vs + vl + vi);
        next.sn_core_cube[j] = SBC_SIX_OVER_PI * vs;
    }
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        if (prev.water_fraction[j] == 1.0f) continue;
        next.diam_melt_snow[j] = 10.0f * std::cbrt(next.diam_melt_snow[j]);
        const float fw = next.water_fraction[j];
        const float v = next.velocity_melt_snow[j];
        const float rw = m0[j] * fw * v * N[j];
        const float ri = m0[j] * (1.0f - fw) * v * N[j] / SBC_CLASS_ICE_DENSITY;
        const float total = ri + rw;
        next.psd_ptype[j] =
            ((ri == 0.0f) || (ri / total < SBC_CLASS_FRACTION)) ? SBCClass::rain
            : ((rw == 0.0f) || (rw / total < SBC_CLASS_FRACTION)) ? ice_class
                                                                  : mix_class;
    }
}

// Branch G (classify:402), defined under "Microphysics: refreezing"
inline void sbc_subfreezing(SBCColumnState& state, const SBCLevelState& prev,
                            SBCLevelState& next, const SBCLevel& level,
                            const SpectralBinColumn& column,
                            const SpectralBinDSD& dsd);

// classify:50-71 over the valid levels, from the cloud top down. A level at
// exactly 0 C is cold.
void sbc_crossings(const SpectralBinColumn& column, const float tice,
                   SBCColumnState& state) {
    const std::ptrdiff_t none = column.top - column.surface + 1;
    state.crossings = 0;
    state.cross_index = none;
    state.ice_nuc_top = none;
    float tw_above = column.wetbulb[column.top];
    std::ptrdiff_t K = 0;
    for (std::ptrdiff_t k = column.top - 1; k >= column.surface; --k) {
        if (!column.valid(k)) continue;
        ++K;
        const float tw = column.wetbulb[k];
        if ((tw <= ZEROCNK) != (tw_above <= ZEROCNK)) {
            if (state.crossings == 0) state.cross_index = K;
            ++state.crossings;
        }
        if ((tw > tice) && (tw_above <= tice) && (state.ice_nuc_top == none)) {
            state.ice_nuc_top = K;
        }
        tw_above = tw;
    }
}

// classify:550-616, without the 3600 s and the bin width of the reference's
// masses, which cancel.
SpectralBinResult sbc_surface(const SBCColumnState& state,
                              const SBCLevelState& surface,
                              const float tw_surface,
                              const SpectralBinDSD& dsd) {
    const float* m0 = dsd.mass().data();
    const float* N = dsd.concentration().data();
    float rainw = 0.0f;
    float raini = 0.0f;
    for (std::ptrdiff_t j = 0; j < dsd.nbins(); ++j) {
        const float fw = surface.water_fraction[j];
        const float v = surface.velocity_melt_snow[j];
        rainw += m0[j] * fw * v * N[j];
        raini += m0[j] * (1.0f - fw) * v * N[j];
    }
    raini /= SBC_ICE_DENSITY;

    const float liquid = rainw / (raini + rainw);
    SpectralBinResult result{PrecipType::missing, liquid, state.slw_hgt};
    const bool warm = (tw_surface > ZEROCNK);
    if (warm && (state.crossings == 1)) {
        if ((raini == 0.0f) || (liquid > 0.85f)) {
            result.precip_type = PrecipType::rain;
        } else if ((rainw == 0.0f) || (liquid < 0.60f)) {
            result.precip_type = PrecipType::snow;
        } else {
            result.precip_type = PrecipType::rain_snow;
        }
    } else if (warm) {
        if ((rainw == 0.0f) || (liquid < 0.15f)) {
            result.precip_type = PrecipType::ice_pellets;
        } else if ((raini == 0.0f) || (raini / (raini + rainw) < 0.15f)) {
            result.precip_type = PrecipType::rain;
        } else {
            result.precip_type = PrecipType::rain_ice_pellets;
        }
    } else if ((rainw == 0.0f) || (liquid < 0.15f)) {
        result.precip_type = PrecipType::ice_pellets;
    } else {
        result.precip_type = ((raini == 0.0f) || (liquid > 0.85f))
                                 ? PrecipType::freezing_rain
                                 : PrecipType::freezing_rain_ice_pellets;
        result.supercooled_liquid_height = 0.0f;
    }
    return result;
}

SpectralBinResult sbc_microphysics(const SpectralBinColumn& column,
                                   const SpectralBinDSD& dsd,
                                   const float ice_nucleation_temperature,
                                   float liquid_fraction_profile[]) {
    const std::ptrdiff_t nbins = dsd.nbins();
    SBCColumnState state;
    state.temp_nuc = ice_nucleation_temperature;
    state.rime_factor = dsd.rime_factor();
    state.snow_volume_min = 1.81e-5f * std::pow(state.rime_factor, 3.26f);
    std::fill_n(state.refrz_lvl_flag.begin(), nbins, false);
    sbc_crossings(column, ice_nucleation_temperature, state);

    // Levels K - 1 and K alternate between the buffers. The cloud top keeps
    // its own record.
    std::array<SBCLevelState, 2> buffers;
    SBCLevelState* prev = &state.top;
    SBCLevelState* next = &state.top;
    const float tw_top = column.wetbulb[column.top];
    std::ptrdiff_t K = 0;
    std::ptrdiff_t k_above = -1;
    for (std::ptrdiff_t k = column.top; k >= column.surface; --k) {
        if (!column.valid(k)) continue;
        const float tw = column.wetbulb[k];
        const SBCLevel level{
            k,
            k_above,
            K,
            tw,
            column.pressure[k] / SBC_DRY_AIR_GAS_CONSTANT / tw * 1.0e-3f,
            (K == 0) ? 0.0f : column.height[k_above] - column.height[k]};
        if (K == 0) {
            next->clear(nbins);
            if (tw_top <= state.temp_nuc) {
                sbc_frozen_cloud_top(state, *next, level, column, dsd);
            } else {
                sbc_liquid_cloud_top(state, *next, level, column, dsd);
            }
        } else {
            const bool above_crossing = (K < state.cross_index);
            const float tice = state.temp_nuc;
            if (above_crossing && (tw_top <= tice) && (tw < ZEROCNK)) {
                sbc_frozen_above_crossing(state, *prev, *next, level, column,
                                          dsd);
            } else {
                next->clear(nbins);
                if (above_crossing && (tice < tw_top) && (tw_top <= ZEROCNK) &&
                    (tw > tice) && (state.ice_nuc_top >= state.cross_index)) {
                    sbc_supercooled_above_crossing(state, *prev, *next, level,
                                                   column, dsd);
                } else if (above_crossing && (tw_top > ZEROCNK) &&
                           (tw > tice)) {
                    sbc_warm_above_crossing(state, *prev, *next, level,
                                            column, dsd);
                } else if (tw >= ZEROCNK) {
                    sbc_melting(state, *prev, *next, level, column, dsd);
                } else {
                    sbc_subfreezing(state, *prev, *next, level, column, dsd);
                }
            }
        }

        if (liquid_fraction_profile != nullptr) {
            std::copy_n(next->water_fraction.begin(), nbins,
                        liquid_fraction_profile + k * nbins);
        }
        if (K == 0) state.set_refrz_lvl(0, state.top, nbins);

        prev = next;
        next = (prev == &buffers[0]) ? &buffers[1] : &buffers[0];
        k_above = k;
        ++K;
    }
    return sbc_surface(state, *prev, column.wetbulb[column.surface], dsd);
}
}  // namespace

// ---------------------------------------------------------------------------
// Microphysics: refreezing
// ---------------------------------------------------------------------------

namespace {
inline void sbc_subfreezing(SBCColumnState& state, const SBCLevelState& prev,
                            SBCLevelState& next, const SBCLevel& level,
                            const SpectralBinColumn& column,
                            const SpectralBinDSD& dsd) {
    // tnuc_alt: Tice below a level where every bin is liquid (K)
    constexpr float TICE_ALT = 263.15f;
    // conduct_ice (J m^-1 s^-1 K^-1)
    constexpr float ICE_CONDUCTIVITY = 2.26f;
    // lh_melt (J kg^-1)
    constexpr float LATENT_HEAT_MELTING = 3.35e5f;
    // 2.85e6 / 461.5, the latent heat of sublimation over the gas constant
    // of water vapor (K)
    constexpr float SUBLIMATION_OVER_RV = 6175.51465f;
    // 2.85e6 * 2.0e-5, the latent heat of sublimation times the vapor
    // diffusivity
    constexpr float SUBLIMATION_DIFFUSION = 57.0f;
    // 0.308 Pr^(1/3), with the Prandtl number Pr = 1.8e-5 / 1.91e-5
    constexpr float PRANDTL_COEFF = 0.301969975f;
    // pi / 6
    constexpr float PI_OVER_SIX = 0.52359879f;
    // 1000 / 0.917, mm^3 of ice per g
    constexpr float ICE_VOLUME_PER_MASS = 1090.51257f;

    const std::ptrdiff_t nbins = dsd.nbins();
    state.refrz_detected_flag = true;
    bool melted = true;
    for (std::ptrdiff_t j = 0; j < nbins; ++j) {
        melted &= (prev.water_fraction[j] == 1.0f);
    }
    if (melted) state.temp_nuc = TICE_ALT;
    // Tw is below 0 C here, so a bin that holds ice refreezes, and at or
    // below Tice every bin does.
    const bool nucleates = (level.wetbulb <= state.temp_nuc);

    const float* D = dsd.diameter().data();
    bool supercooled = false;

    // classify:500-539: a liquid bin above Tice falls unchanged.
    if (!nucleates) {
        for (std::ptrdiff_t j = 0; j < nbins; ++j) {
            if (prev.water_fraction[j] != 1.0f) continue;
            next.water_fraction[j] = prev.water_fraction[j];
            next.velocity_melt_snow[j] = prev.velocity_melt_snow[j];
            next.mass_water[j] = prev.mass_water[j];
            next.mass_snow[j] = prev.mass_snow[j];
            next.volume_liq[j] = prev.volume_liq[j];
            next.volume_ice[j] = prev.volume_ice[j];
            next.volume_snow[j] = prev.volume_snow[j];
            next.diam_melt_snow[j] = prev.diam_melt_snow[j];
            const SBCClass ptype = prev.psd_ptype[j];
            const bool rain = (ptype == SBCClass::rain);
            const bool mixed = (ptype == SBCClass::rain_snow) ||
                               (ptype == SBCClass::rain_ice_pellets);
            next.psd_ptype[j] =
                (rain || mixed) ? sbc_freezing_class(D[j], mixed) : ptype;
            supercooled |= (rain || mixed);
        }
    }

    // classify:412-497: the other bins refreeze.
    if (nucleates || !melted) {
        const float tw = level.wetbulb;
        const float tw_c = sbc_celsius(tw);
        const float undercooling = -tw_c;
        // deriv_rho_ice (kg m^-3 K^-1)
        const float deriv_rho_ice = (3.8f + 0.25f * tw_c) * 1.0e-4f;
        // svp_wrt_ice (Pa)
        const float svp_ice =
            611.0f * std::exp(SUBLIMATION_OVER_RV * tw_c / (ZEROCNK * tw));
        // abs_humid_ice, with the dry-air density p / (R_d T) (kg m^-3)
        const float abs_humid_ice =
            svp_ice * 0.622f /
            (SBC_DRY_AIR_GAS_CONSTANT * column.temperature[level.k]);
        // unknwn_xsi and the vapor term of the numerator, over the
        // refreezing Prandtl number
        const float xsi_factor =
            SBC_AIR_CONDUCTIVITY + SUBLIMATION_DIFFUSION * deriv_rho_ice;
        const float vapor_factor = SUBLIMATION_DIFFUSION *
                                   (1.0f - column.relh[level.k]) *
                                   abs_humid_ice;
        const float latent_factor =
            LATENT_HEAT_MELTING * (1.0f + 0.012f * sbc_celsius(state.temp_nuc));
        const float density_correction = level.density_correction();
        const bool surface = (level.k == column.surface);
        next.sn_core_scale = 1.0f;

        const float* m0 = dsd.mass().data();
        const float* N = dsd.concentration().data();
        const float* v_pellet = dsd.pellet_fall_speed().data();
        for (std::ptrdiff_t j = 0; j < nbins; ++j) {
            const float fw_prev = prev.water_fraction[j];
            if (!nucleates && (fw_prev == 1.0f)) continue;
            // A bin that starts refreezing moves the refreeze level of every
            // bin to K - 1. prev holds level K - 1 for the bins before it too.
            if (!state.refrz_lvl_flag[j]) {
                state.refrz_lvl_flag[j] = true;
                if (state.refrz_lvl != level.K - 1) {
                    state.set_refrz_lvl(level.K - 1, prev, nbins);
                }
            }

            const float fw_refrz = state.refrz_water_fraction[j];
            const float v_ground = v_pellet[j] * density_correction;
            const float ratio =
                (fw_prev == fw_refrz) ? 1.0f : fw_prev / fw_refrz;
            const float v =
                v_ground +
                (state.refrz_velocity_melt_snow[j] - v_ground) * ratio;
            const float dm_prev = prev.diam_melt_snow[j];
            const float prandtl =
                0.78f + PRANDTL_COEFF *
                            std::sqrt(dm_prev * SBC_REYNOLDS_COEFF * v);
            const float xsi = prandtl * xsi_factor;
            const float numerator = std::max(
                2.0e-3f * dm_prev * ICE_CONDUCTIVITY *
                    (undercooling * xsi + prandtl * vapor_factor),
                0.0f);
            // The reference adds xsi D / D_w, with the diameter D of this
            // level read before it is set, so the term is 0.
            const float denominator =
                v * latent_factor * (ICE_CONDUCTIVITY - xsi);
            const float change =
                1000.0f * numerator / denominator * level.layer_depth;

            const float mw = std::max(prev.mass_water[j] - change, 0.0f);
            const float fw = mw / m0[j];
            // rfrz_volume_snow, rfrz_volume_water, and rfrz_volume_ice (mm^3)
            const float dc = state.refrz_in_diam_sn_core[j];
            const float vs = PI_OVER_SIX * dc * dc * dc;
            const float vw = std::max(1000.0f * (m0[j] * fw_refrz - change),
                                      0.0f);
            const float vi = std::max(ICE_VOLUME_PER_MASS * change, 0.0f);

            next.water_fraction[j] = fw;
            next.velocity_melt_snow[j] = v;
            next.mass_water[j] = mw;
            next.mass_snow[j] = m0[j] - mw;
            // Its cube root follows, in a loop of its own.
            next.diam_melt_snow[j] = SBC_SIX_OVER_PI * (vs + vw + vi);
            next.sn_core_cube[j] = SBC_SIX_OVER_PI * (vs + vw);
        }
        for (std::ptrdiff_t j = 0; j < nbins; ++j) {
            if (!nucleates && (prev.water_fraction[j] == 1.0f)) continue;
            // Nothing reads the diameter at the surface.
            if (!surface) {
                next.diam_melt_snow[j] = std::cbrt(next.diam_melt_snow[j]);
            }
            const float fw = next.water_fraction[j];
            const float v = next.velocity_melt_snow[j];
            const float rw = m0[j] * fw * v * N[j];
            const float ri =
                m0[j] * (1.0f - fw) * v * N[j] / SBC_CLASS_ICE_DENSITY;
            const float total = ri + rw;
            const SBCClass ptype =
                ((ri == 0.0f) || (ri / total < SBC_CLASS_FRACTION))
                    ? sbc_freezing_class(D[j], false)
                : ((rw == 0.0f) || (rw / total < SBC_CLASS_FRACTION))
                    ? SBCClass::ice_pellets
                    : sbc_freezing_class(D[j], true);
            next.psd_ptype[j] = ptype;
            supercooled |= (ptype != SBCClass::ice_pellets);
        }
    }

    if (supercooled) state.slw_hgt = column.height_agl(level.k);
}
}  // namespace

// ---------------------------------------------------------------------------
// Microphysics: liquid cloud tops
// ---------------------------------------------------------------------------

namespace {
// Foote and du Toit (1969) scale their fall speed by exp(z / this) (m).
constexpr float SBC_FOOTE_DU_TOIT_SCALE_HEIGHT = 20000.0f;

// Branches D and E: the drops of a liquid cloud top fall unchanged, at the
// Foote and du Toit fall speed. Each branch sets its own classes.
inline void sbc_liquid_top_drops(const SBCLevelState& top, SBCLevelState& next,
                                 const SBCLevel& level,
                                 const SpectralBinColumn& column,
                                 const SpectralBinDSD& dsd) {
    const float height_factor =
        std::exp(column.height[level.k] / SBC_FOOTE_DU_TOIT_SCALE_HEIGHT);
    const float* diameter = dsd.diameter().data();
    const float* v_drop = dsd.foote_du_toit_fall_speed().data();
    for (std::ptrdiff_t j = 0; j < dsd.nbins(); ++j) {
        next.water_fraction[j] = top.water_fraction[j];
        next.velocity_melt_snow[j] = v_drop[j] * height_factor;
        next.mass_water[j] = top.mass_water[j];
        next.volume_liq[j] = top.volume_liq[j];
        next.volume_ice[j] = top.volume_ice[j];
        next.volume_snow[j] = top.volume_snow[j];
        next.diam_melt_snow[j] = diameter[j];
    }
}

// Branch C (classify:186): a liquid cloud top
inline void sbc_liquid_cloud_top(SBCColumnState& state, SBCLevelState& next,
                                 const SBCLevel& level,
                                 const SpectralBinColumn& column,
                                 const SpectralBinDSD& dsd) {
    state.uknwn_aa = dsd.liquid_aa().data();
    state.uknwn_bb = dsd.liquid_bb().data();
    const float density_correction = level.density_correction();
    const float* diameter = dsd.diameter().data();
    const float* m0 = dsd.mass().data();
    const float* v_rain = dsd.rain_fall_speed().data();
    for (std::ptrdiff_t j = 0; j < dsd.nbins(); ++j) {
        next.water_fraction[j] = 1.0f;
        next.velocity_melt_snow[j] = v_rain[j] * density_correction;
        next.mass_water[j] = m0[j];
        // A liquid density of 1 g cm^-3
        next.volume_liq[j] = m0[j];
        next.diam_melt_snow[j] = diameter[j];
        next.psd_ptype[j] = sbc_freezing_class(diameter[j], false);
    }
    state.slw_hgt = column.height_agl(level.k);
}

// Branch D (classify:217): the drops of a cloud top at or below 0 C, above
// the first crossing
inline void sbc_supercooled_above_crossing(
    SBCColumnState& state, [[maybe_unused]] const SBCLevelState& prev,
    SBCLevelState& next, const SBCLevel& level,
    const SpectralBinColumn& column, const SpectralBinDSD& dsd) {
    sbc_liquid_top_drops(state.top, next, level, column, dsd);
    const float* diameter = dsd.diameter().data();
    for (std::ptrdiff_t j = 0; j < dsd.nbins(); ++j) {
        next.psd_ptype[j] = sbc_freezing_class(diameter[j], false);
    }
    state.slw_hgt = column.height_agl(level.k);
}

// Branch E (classify:243): the drops of a cloud top above 0 C, above the
// first crossing
inline void sbc_warm_above_crossing(
    SBCColumnState& state, [[maybe_unused]] const SBCLevelState& prev,
    SBCLevelState& next, const SBCLevel& level,
    const SpectralBinColumn& column, const SpectralBinDSD& dsd) {
    sbc_liquid_top_drops(state.top, next, level, column, dsd);
    std::fill_n(next.psd_ptype.begin(), dsd.nbins(), SBCClass::rain);
}
}  // namespace

// ---------------------------------------------------------------------------
// Precipitation type from a full sounding
// ---------------------------------------------------------------------------

SpectralBinResult spectral_bin_classifier(
    const float pressure[], const float height[], const float temperature[],
    const float dewpoint[], const float relh[], const float wetbulb[],
    const std::ptrdiff_t N, const SpectralBinDSD& dsd,
    const float ice_nucleation_temperature, float liquid_fraction_profile[]) {
    const float cloud_top =
        spectral_bin_cloud_top(pressure, height, temperature, dewpoint, relh, N);
    return spectral_bin_classifier(pressure, height, temperature, dewpoint,
                                   relh, wetbulb, N, cloud_top, dsd,
                                   ice_nucleation_temperature,
                                   liquid_fraction_profile);
}

}  // namespace sharp
