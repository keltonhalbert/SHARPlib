#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/params/winter.h>
#include <SHARPlib/thermo.h>

#include <cmath>
#include <limits>
#include <vector>

#include "doctest.h"

// ===========================================================================
// Precipitation type: the modified Bourgouin method (Birk et al. 2021)
// ===========================================================================

TEST_CASE("Testing precipitation-type result types default to MISSING") {
    const sharp::BourgouinEnergy energy;
    CHECK(energy.melting_energy_total == sharp::MISSING);
    CHECK(energy.melting_energy_aloft == sharp::MISSING);
    CHECK(energy.refreezing_energy == sharp::MISSING);

    const sharp::PrecipTypeProbabilities probs;
    CHECK(probs.rain == sharp::MISSING);
    CHECK(probs.snow == sharp::MISSING);
    CHECK(probs.freezing_rain == sharp::MISSING);
    CHECK(probs.ice_pellets == sharp::MISSING);
}

// ---------------------------------------------------------------------------
// Wet-bulb melting and refreezing energies from a sounding
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Precipitation generation layer from a sounding
// ---------------------------------------------------------------------------

namespace {
// Air temperature (K) of the warm test soundings, where relative humidity
// is over liquid water.
constexpr float WARM_AIR = 283.15f;

// Relative humidity as precipitation_generation_layer computes it: over ice
// below 0 C, over liquid otherwise.
float switched_relh(const float pressure, const float temperature,
                    const float dewpoint) {
    return (temperature < sharp::ZEROCNK)
               ? sharp::relative_humidity_ice(pressure, temperature, dewpoint)
               : sharp::relative_humidity(pressure, temperature, dewpoint);
}

struct RelhSounding {
    std::vector<float> pressure;
    std::vector<float> height;
    std::vector<float> temperature;
    std::vector<float> dewpoint;

    sharp::HeightLayer generation_layer(const float min_depth = 0.0f) const {
        return sharp::precipitation_generation_layer(
            pressure.data(), height.data(), temperature.data(), dewpoint.data(),
            static_cast<std::ptrdiff_t>(height.size()), min_depth);
    }
};

// A sounding whose relative humidity at each height is exactly relh, so
// levels meant to be at 75 % are exactly at the threshold. One float step in
// temperature or dewpoint moves the relative humidity by several float
// steps, so each level nudges its temperature up from near_temperature
// until some dewpoint gives relh exactly.
RelhSounding relh_sounding(const std::vector<float>& height,
                           const std::vector<float>& relh,
                           const std::vector<float>& near_temperature) {
    RelhSounding snd;
    for (std::size_t k = 0; k < height.size(); ++k) {
        const float pres = 100000.0f - 10.0f * height[k];
        float tmpk = near_temperature[k];
        float dwpk = sharp::MISSING;
        for (int tries = 0; (tries < 10000) && (dwpk == sharp::MISSING);
             ++tries) {
            // Bisect for the lowest dewpoint whose relative humidity is at
            // or above relh.
            float lo = tmpk - 80.0f;
            float hi = tmpk;
            for (float mid = 0.5f * (lo + hi); (mid != lo) && (mid != hi);
                 mid = 0.5f * (lo + hi)) {
                if (switched_relh(pres, tmpk, mid) < relh[k]) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            if (switched_relh(pres, tmpk, hi) == relh[k]) {
                dwpk = hi;
            } else {
                tmpk = std::nextafter(tmpk, 400.0f);
            }
        }
        REQUIRE(dwpk != sharp::MISSING);
        snd.pressure.push_back(pres);
        snd.height.push_back(height[k]);
        snd.temperature.push_back(tmpk);
        snd.dewpoint.push_back(dwpk);
    }
    return snd;
}

// Alternating moist and dry layers of the given depths (m), from the
// surface up, at WARM_AIR. Each layer has a level in its middle, at 0.9 if
// moist and 0.5 if dry. The boundaries between layers are levels at exactly
// 0.75, which continue the layer below them, so every layer ends exactly at
// its boundary level.
RelhSounding layered_sounding(bool moist, const std::vector<float>& depths) {
    std::vector<float> height = {0.0f};
    std::vector<float> relh = {moist ? 0.9f : 0.5f};
    float bottom = 0.0f;
    for (std::size_t i = 0; i < depths.size(); ++i) {
        const float top = bottom + depths[i];
        const float layer_relh = moist ? 0.9f : 0.5f;
        height.push_back(0.5f * (bottom + top));
        relh.push_back(layer_relh);
        height.push_back(top);
        relh.push_back((i + 1 < depths.size()) ? 0.75f : layer_relh);
        bottom = top;
        moist = !moist;
    }
    return relh_sounding(height, relh,
                         std::vector<float>(height.size(), WARM_AIR));
}

void check_layer(const sharp::HeightLayer& layer, const float bottom,
                 const float top) {
    CHECK(layer.bottom == bottom);
    CHECK(layer.top == top);
}
}  // namespace

TEST_CASE("Testing precipitation_generation_layer depth thresholds") {
    // A generation layer is a moist layer deeper than 1000 m: 1100 m is one,
    // and 900 m and exactly 1000 m are not.
    auto snd = layered_sounding(false, {500.0f, 1100.0f, 500.0f});
    check_layer(snd.generation_layer(), 500.0f, 1600.0f);
    snd = layered_sounding(false, {500.0f, 900.0f, 500.0f});
    check_layer(snd.generation_layer(), sharp::MISSING, sharp::MISSING);
    snd = layered_sounding(false, {500.0f, 1000.0f, 500.0f});
    check_layer(snd.generation_layer(), sharp::MISSING, sharp::MISSING);

    // A dry layer deeper than 1500 m eliminates the layers above it. 1400 m
    // and exactly 1500 m don't, and 1600 m does.
    snd = layered_sounding(true, {1200.0f, 1400.0f, 1200.0f});
    check_layer(snd.generation_layer(), 2600.0f, 3800.0f);
    snd = layered_sounding(true, {1200.0f, 1500.0f, 1200.0f});
    check_layer(snd.generation_layer(), 2700.0f, 3900.0f);
    snd = layered_sounding(true, {1200.0f, 1600.0f, 1200.0f});
    check_layer(snd.generation_layer(), 0.0f, 1200.0f);
}

TEST_CASE("Testing precipitation_generation_layer picks the highest layer") {
    // Generation layers at 0-1200 m and 1700-2800 m, then an 800 m moist
    // layer that is too shallow. The highest generation layer wins.
    auto snd =
        layered_sounding(true, {1200.0f, 500.0f, 1100.0f, 300.0f, 800.0f});
    check_layer(snd.generation_layer(), 1700.0f, 2800.0f);

    // A 1300 m moist layer on top, above a 1600 m dry layer, is eliminated.
    // Above a 1400 m dry layer it is the highest generation layer.
    snd = layered_sounding(
        true, {1200.0f, 500.0f, 1100.0f, 300.0f, 800.0f, 1600.0f, 1300.0f});
    check_layer(snd.generation_layer(), 1700.0f, 2800.0f);
    snd = layered_sounding(
        true, {1200.0f, 500.0f, 1100.0f, 300.0f, 800.0f, 1400.0f, 1300.0f});
    check_layer(snd.generation_layer(), 5300.0f, 6600.0f);

    // A deep dry layer at the surface eliminates everything above it.
    snd = layered_sounding(false, {1600.0f, 1200.0f});
    check_layer(snd.generation_layer(), sharp::MISSING, sharp::MISSING);
    snd = layered_sounding(false, {1400.0f, 1200.0f});
    check_layer(snd.generation_layer(), 1400.0f, 2600.0f);
}

TEST_CASE("Testing precipitation_generation_layer relative humidity phase") {
    const float pres[] = {100000.0f, 90000.0f, 80000.0f};
    const float hght[] = {0.0f, 1000.0f, 2000.0f};

    // At or above 0 C, relative humidity is over liquid. T 283.15 K with
    // Td 280 K is moist over liquid, though dry over ice.
    const float warm_tmpk[] = {WARM_AIR, WARM_AIR, WARM_AIR};
    const float warm_dwpk[] = {280.0f, 280.0f, 280.0f};
    CHECK(sharp::relative_humidity(pres[0], WARM_AIR, 280.0f) ==
          doctest::Approx(0.808).epsilon(1e-3));
    CHECK(sharp::relative_humidity_ice(pres[0], WARM_AIR, 280.0f) ==
          doctest::Approx(0.733).epsilon(1e-3));
    check_layer(sharp::precipitation_generation_layer(pres, hght, warm_tmpk,
                                                      warm_dwpk, 3),
                0.0f, 2000.0f);

    // Below 0 C, relative humidity is over ice. T 263.15 K with Td 259.15 K
    // is moist over ice, though dry over liquid.
    const float cold_tmpk[] = {263.15f, 263.15f, 263.15f};
    const float cold_dwpk[] = {259.15f, 259.15f, 259.15f};
    CHECK(sharp::relative_humidity(pres[0], 263.15f, 259.15f) ==
          doctest::Approx(0.725).epsilon(1e-3));
    CHECK(sharp::relative_humidity_ice(pres[0], 263.15f, 259.15f) ==
          doctest::Approx(0.801).epsilon(1e-3));
    check_layer(sharp::precipitation_generation_layer(pres, hght, cold_tmpk,
                                                      cold_dwpk, 3),
                0.0f, 2000.0f);
}

TEST_CASE("Testing precipitation_generation_layer levels at exactly 75 %") {
    // Levels at exactly 75 % continue the moist layer, so this is one
    // 2100 m generation layer. Read strictly, the paper would give two
    // 100 m moist layers and no generation layer.
    const auto snd = relh_sounding({0.0f, 100.0f, 2000.0f, 2100.0f},
                                   {0.8f, 0.75f, 0.75f, 0.8f},
                                   {WARM_AIR, WARM_AIR, WARM_AIR, WARM_AIR});
    check_layer(snd.generation_layer(), 0.0f, 2100.0f);
}

TEST_CASE("Testing precipitation_generation_layer min_depth") {
    // A 50 m dry sliver splits a moist layer into two 600 m layers, neither
    // deep enough. With min_depth = 100 m it is absorbed, leaving one
    // 1250 m generation layer.
    auto snd = layered_sounding(true, {600.0f, 50.0f, 600.0f});
    check_layer(snd.generation_layer(), sharp::MISSING, sharp::MISSING);
    check_layer(snd.generation_layer(100.0f), 0.0f, 1250.0f);

    // A 50 m moist sliver splits a dry layer into 1400 m and 200 m. With
    // min_depth = 100 m it is absorbed, and the merged 1650 m dry layer
    // eliminates the moist layer above it.
    snd = layered_sounding(true, {1200.0f, 1400.0f, 50.0f, 200.0f, 1200.0f});
    check_layer(snd.generation_layer(), 2850.0f, 4050.0f);
    check_layer(snd.generation_layer(100.0f), 0.0f, 1200.0f);

    // Shallow runs between a dry layer and the moist layer above it go to
    // the side with more of their depth, so the dry layer's depth is only
    // known once the moist layer is min_depth deep. 50 m moist and then
    // 80 m dry make the dry layer 1530 m, which eliminates the moist layer.
    snd = layered_sounding(true, {1200.0f, 1400.0f, 50.0f, 80.0f, 1200.0f});
    check_layer(snd.generation_layer(), 2730.0f, 3930.0f);
    check_layer(snd.generation_layer(100.0f), 0.0f, 1200.0f);

    // 80 m moist and then 50 m dry go to the moist layer instead, leaving a
    // 1400 m dry layer and a 1330 m generation layer.
    snd = layered_sounding(true, {1200.0f, 1400.0f, 80.0f, 50.0f, 1200.0f});
    check_layer(snd.generation_layer(100.0f), 2600.0f, 3930.0f);
}

TEST_CASE("Testing precipitation_generation_layer with no or one level") {
    // N = 0 reads no element.
    check_layer(sharp::precipitation_generation_layer(nullptr, nullptr, nullptr,
                                                      nullptr, 0),
                sharp::MISSING, sharp::MISSING);

    const float pres = 100000.0f;
    const float hght = 0.0f;
    const float tmpk = WARM_AIR;
    check_layer(
        sharp::precipitation_generation_layer(&pres, &hght, &tmpk, &tmpk, 1),
        sharp::MISSING, sharp::MISSING);
}

#ifndef NO_QC
TEST_CASE("Testing precipitation_generation_layer skips missing levels") {
    // Below 0 C: relative humidity over ice 0.6 and 0.7, then 0.9 from
    // 1000 m up. With the 1000 m temperature missing, the walk joins 500 m
    // to 1500 m and crosses 75 % at 750 m instead of 625 m.
    auto snd = relh_sounding({0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f},
                             {0.6f, 0.7f, 0.9f, 0.9f, 0.9f, 0.9f},
                             {258.0f, 258.0f, 258.0f, 268.0f, 268.0f, 268.0f});
    auto layer = snd.generation_layer();
    CHECK(layer.bottom == doctest::Approx(625.0f));
    CHECK(layer.top == 2500.0f);

    for (const float bad :
         {sharp::MISSING, std::numeric_limits<float>::quiet_NaN()}) {
        CAPTURE(bad);
        snd.temperature[2] = bad;
        layer = snd.generation_layer();
        CHECK(layer.bottom == doctest::Approx(750.0f));
        CHECK(layer.top == 2500.0f);
    }

    // A missing temperature inside the layer as well changes nothing.
    snd.temperature[4] = sharp::MISSING;
    layer = snd.generation_layer();
    CHECK(layer.bottom == doctest::Approx(750.0f));
    CHECK(layer.top == 2500.0f);
}
#endif

// ---------------------------------------------------------------------------
// Probability of ice, and precipitation-type probabilities from energies
// ---------------------------------------------------------------------------

TEST_CASE("Testing probability_of_ice") {
    // The constant pieces of Eq. 2. -15 C and -7 C are exact in float K.
    CHECK(sharp::probability_of_ice(sharp::ZEROCNK - 15.0f) == 1.0f);
    CHECK(sharp::probability_of_ice(sharp::ZEROCNK - 40.0f) == 1.0f);
    CHECK(sharp::probability_of_ice(sharp::ZEROCNK - 7.0f) == 0.0f);
    CHECK(sharp::probability_of_ice(sharp::ZEROCNK + 10.0f) == 0.0f);

    // The Eq. 2b polynomial: about 58 % at -11 C (Fig. 2), and 260.5 K
    // (-12.65 C).
    CHECK(sharp::probability_of_ice(sharp::ZEROCNK - 11.0f) ==
          doctest::Approx(0.583474));
    CHECK(sharp::probability_of_ice(260.5f) == doctest::Approx(0.7286607));

    // The documented jumps. Just inside the polynomial's interval it gives
    // 98.3 % at -15 C and 0.81 % at -7 C, not 100 % and 0 %.
    const float just_above_m15 =
        std::nextafter(sharp::ZEROCNK - 15.0f, sharp::ZEROCNK);
    const float just_below_m7 = std::nextafter(sharp::ZEROCNK - 7.0f, 0.0f);
    CHECK(sharp::probability_of_ice(just_above_m15) ==
          doctest::Approx(0.98325).epsilon(1e-4));
    CHECK(sharp::probability_of_ice(just_below_m7) ==
          doctest::Approx(0.008082).epsilon(1e-4));

    // Celsius input is misread, as documented.
    CHECK(sharp::probability_of_ice(-12.0f) == sharp::MISSING);
    CHECK(sharp::probability_of_ice(5.0f) == 1.0f);

    // MISSING and NaN give MISSING in every build.
    CHECK(sharp::probability_of_ice(sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::probability_of_ice(std::numeric_limits<float>::quiet_NaN()) ==
          sharp::MISSING);
}

namespace {
constexpr float COLD_SFC = sharp::ZEROCNK - 2.0f;
constexpr float WARM_SFC = sharp::ZEROCNK + 2.0f;

void check_all_missing(const sharp::PrecipTypeProbabilities& probs) {
    CHECK(probs.rain == sharp::MISSING);
    CHECK(probs.snow == sharp::MISSING);
    CHECK(probs.freezing_rain == sharp::MISSING);
    CHECK(probs.ice_pellets == sharp::MISSING);
}
}  // namespace

TEST_CASE("Testing modified_bourgouin clamps Eq. 7 before the taper") {
    // RE 0, ME 3, ProbIce 1. Eq. 7 gives 458.6 %, clamped to 100 %, then
    // tapered by 0.2 * 3 to 60 %. Tapering before the clamp would give 100 %.
    const auto probs =
        sharp::modified_bourgouin({3.0f, 0.0f, 0.0f}, 1.0f, WARM_SFC);
    CHECK(probs.rain == doctest::Approx(0.60));
    CHECK(probs.freezing_rain == 0.0f);
    CHECK(probs.snow == 1.0f);  // Eq. 9 gives 645 %
    CHECK(probs.ice_pellets == 0.0f);
}

TEST_CASE("Testing modified_bourgouin freezing rain uses the total ME") {
    // An absorbed warm sliver: ME_total 4, ME_aloft 2.01, RE 180, ProbIce 1,
    // cold surface. Eq. 7 with ME_total gives 80.8 %, tapered by 0.2 * 4 to
    // 64.64 %. Using ME_aloft instead would give about 32.32 %.
    const auto probs =
        sharp::modified_bourgouin({4.0f, 2.01f, 180.0f}, 1.0f, COLD_SFC);
    CHECK(probs.freezing_rain == doctest::Approx(0.6464));
    CHECK(probs.rain == 0.0f);
}

TEST_CASE("Testing the modified_bourgouin ice pellet gate") {
    const float tiny = std::numeric_limits<float>::min();

    // No melting aloft: no ice pellets, however much refreezing.
    auto probs =
        sharp::modified_bourgouin({10.0f, 0.0f, 10.0f}, 1.0f, COLD_SFC);
    CHECK(probs.ice_pellets == 0.0f);

    // A tiny positive ME_aloft jumps to the Eq. 8 value, 2.3 * 10 + 3 = 26 %.
    probs = sharp::modified_bourgouin({10.0f, tiny, 10.0f}, 1.0f, COLD_SFC);
    CHECK(probs.ice_pellets == doctest::Approx(0.26));

    // No refreezing: no ice pellets. Eq. 8 alone would give 3 %.
    probs = sharp::modified_bourgouin({10.0f, tiny, 0.0f}, 1.0f, COLD_SFC);
    CHECK(probs.ice_pellets == 0.0f);
}

TEST_CASE("Testing modified_bourgouin branches") {
    // Unclamped values. ME 15 and RE 200: snow 19.8765 %, Eq. 7 gives 41 %
    // with no taper, and Eq. 8 gives 346.6 %, clamped to 100 %.
    auto probs =
        sharp::modified_bourgouin({15.0f, 15.0f, 200.0f}, 1.0f, COLD_SFC);
    CHECK(probs.snow == doctest::Approx(0.198765));
    CHECK(probs.freezing_rain == doctest::Approx(0.41));
    CHECK(probs.ice_pellets == 1.0f);
    CHECK(probs.rain == 0.0f);

    // Snow uses the total melting energy. With ME_aloft (2 J/kg) it would
    // be clamped to 100 %.
    probs = sharp::modified_bourgouin({15.0f, 2.0f, 50.0f}, 1.0f, COLD_SFC);
    CHECK(probs.snow == doctest::Approx(0.198765));

    // ProbIce 0.5 scales snow and ice pellets, and Eq. 7 becomes
    // 50 + 0.5 * 41 = 70.5 %.
    probs = sharp::modified_bourgouin({15.0f, 15.0f, 200.0f}, 0.5f, COLD_SFC);
    CHECK(probs.snow == doctest::Approx(0.0993825));
    CHECK(probs.freezing_rain == doctest::Approx(0.705));
    CHECK(probs.ice_pellets == doctest::Approx(0.5));

    // ProbIce 0: no snow or ice pellets, and all liquid.
    probs = sharp::modified_bourgouin({15.0f, 15.0f, 200.0f}, 0.0f, COLD_SFC);
    CHECK(probs.snow == 0.0f);
    CHECK(probs.ice_pellets == 0.0f);
    CHECK(probs.freezing_rain == 1.0f);

    // Eq. 8 unclamped: 2.3 * 50 - 42 ln(16) + 3 = 1.55127 %. Eq. 7 gives
    // 356 %, clamped to 100 %.
    probs = sharp::modified_bourgouin({15.0f, 15.0f, 50.0f}, 1.0f, COLD_SFC);
    CHECK(probs.ice_pellets == doctest::Approx(0.0155127));
    CHECK(probs.freezing_rain == 1.0f);

    // Lower clamps. Eq. 8 gives 2.3 * 10 - 42 ln(101) + 3 = -167.8 %, and
    // Eq. 7 gives -2.1 * 300 + 0.2 * 100 + 458 = -152 %.
    probs = sharp::modified_bourgouin({100.0f, 100.0f, 10.0f}, 1.0f, COLD_SFC);
    CHECK(probs.ice_pellets == 0.0f);
    probs = sharp::modified_bourgouin({100.0f, 100.0f, 300.0f}, 1.0f, COLD_SFC);
    CHECK(probs.freezing_rain == 0.0f);

    // The taper applies below 5 J/kg only, and is continuous there.
    probs = sharp::modified_bourgouin({5.0f, 0.0f, 0.0f}, 1.0f, WARM_SFC);
    CHECK(probs.rain == 1.0f);
    probs = sharp::modified_bourgouin({4.9f, 0.0f, 0.0f}, 1.0f, WARM_SFC);
    CHECK(probs.rain == doctest::Approx(0.98));

    // Surface wet-bulb: above 0 C is rain; exactly 0 C and below is
    // freezing rain.
    const sharp::BourgouinEnergy energy = {15.0f, 15.0f, 200.0f};
    probs = sharp::modified_bourgouin(energy, 1.0f,
                                      std::nextafter(sharp::ZEROCNK, WARM_SFC));
    CHECK(probs.rain == doctest::Approx(0.41));
    CHECK(probs.freezing_rain == 0.0f);
    probs = sharp::modified_bourgouin(energy, 1.0f, sharp::ZEROCNK);
    CHECK(probs.rain == 0.0f);
    CHECK(probs.freezing_rain == doctest::Approx(0.41));

    // A warm surface does not suppress ice pellets:
    // 2.3 * 50 - 42 ln(11) + 3 = 17.2884 %.
    probs = sharp::modified_bourgouin({200.0f, 10.0f, 50.0f}, 1.0f, WARM_SFC);
    CHECK(probs.rain == 1.0f);
    CHECK(probs.ice_pellets == doctest::Approx(0.172884));
}

TEST_CASE("Testing modified_bourgouin MISSING and NaN inputs") {
    const sharp::BourgouinEnergy valid = {15.0f, 15.0f, 200.0f};
    for (const float bad :
         {sharp::MISSING, std::numeric_limits<float>::quiet_NaN()}) {
        CAPTURE(bad);
        check_all_missing(
            sharp::modified_bourgouin({bad, 15.0f, 200.0f}, 1.0f, COLD_SFC));
        check_all_missing(
            sharp::modified_bourgouin({15.0f, bad, 200.0f}, 1.0f, COLD_SFC));
        check_all_missing(
            sharp::modified_bourgouin({15.0f, 15.0f, bad}, 1.0f, COLD_SFC));
        check_all_missing(sharp::modified_bourgouin(valid, bad, COLD_SFC));
        check_all_missing(sharp::modified_bourgouin(valid, 1.0f, bad));
    }
}

TEST_CASE("Testing the modified_bourgouin reconstruction check") {
    // Along the paper's zone boundaries (Eq. 6 and the 20 % lines of section
    // 4a), the equations should give the probabilities the paper quotes. The
    // ranges are the reconstructed values in the plan, to the whole percent,
    // for ME from 5 J/kg (no taper) to 200 J/kg. Cold surface, ProbIce 1,
    // and ME_total = ME_aloft = ME.
    for (float me = 5.0f; me <= 200.0f; me += 5.0f) {
        CAPTURE(me);
        const float ln_me = std::log(me + 1.0f);

        // FZRA at RE = 170 + 0.08 ME. Paper 100 %; 101-107 %, clamped.
        auto probs = sharp::modified_bourgouin({me, me, 170.0f + 0.08f * me},
                                               1.0f, COLD_SFC);
        CHECK(probs.freezing_rain == 1.0f);

        // FZRA at RE = 208 + 0.08 ME. Paper about 20 %; 21-28 %.
        probs = sharp::modified_bourgouin({me, me, 208.0f + 0.08f * me}, 1.0f,
                                          COLD_SFC);
        CHECK(probs.freezing_rain >= 0.205f);
        CHECK(probs.freezing_rain < 0.285f);

        // PL at RE = 41 + 17.9 ln(ME + 1). Paper about 100 %; 93-97 %.
        probs = sharp::modified_bourgouin({me, me, 41.0f + 17.9f * ln_me}, 1.0f,
                                          COLD_SFC);
        CHECK(probs.ice_pellets >= 0.925f);
        CHECK(probs.ice_pellets < 0.975f);

        // PL at RE = 7 + 17.9 ln(ME + 1). Paper about 20 %; 15-18 %.
        probs = sharp::modified_bourgouin({me, me, 7.0f + 17.9f * ln_me}, 1.0f,
                                          COLD_SFC);
        CHECK(probs.ice_pellets >= 0.145f);
        CHECK(probs.ice_pellets < 0.185f);
    }

    // SN at ME = 9 and 15. Paper about 100 % and 20 %; 113 % (clamped) and
    // 19.9 %.
    auto probs = sharp::modified_bourgouin({9.0f, 9.0f, 0.0f}, 1.0f, COLD_SFC);
    CHECK(probs.snow == 1.0f);
    probs = sharp::modified_bourgouin({15.0f, 15.0f, 0.0f}, 1.0f, COLD_SFC);
    CHECK(probs.snow == doctest::Approx(0.199).epsilon(1e-3));
}

// ---------------------------------------------------------------------------
// Precipitation type from a full sounding
// ---------------------------------------------------------------------------
