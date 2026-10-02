#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/params/winter.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
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

namespace {
// Eq. 1: the energy (J/kg) of an area between Tw and 0 C (K m).
constexpr float energy_of_area(const float area) {
    return sharp::GRAVITY / sharp::ZEROCNK * area;
}

// A test sounding of up to 32 levels. Pressure falls 7 Pa per meter from
// 1000 hPa, so every level below 10.7 km is at or above 250 hPa.
struct TwSounding {
    static constexpr std::ptrdiff_t capacity = 32;
    float pres[capacity];
    float hght[capacity];
    float wetbulb[capacity];
    std::ptrdiff_t N = 0;

    // Adds a level at height z (m) with the wet-bulb dtw (K) from 0 C.
    void add(const float z, const float dtw) {
        REQUIRE(N < capacity);
        pres[N] = 100000.0f - 7.0f * z;
        hght[N] = z;
        wetbulb[N] = sharp::ZEROCNK + dtw;
        ++N;
    }

    sharp::BourgouinEnergy energy(const float min_energy = 0.0f) const {
        return sharp::bourgouin_energy(pres, hght, wetbulb, N, min_energy);
    }
};

// Builds a sounding from the energy of each layer between 0 C crossings,
// from the surface up (J/kg; negative for refreezing). Each layer is a
// triangle that peaks 4 K from 0 C and is as deep as its energy needs, so
// every crossing falls on a level. With surface_peak, the lowest layer
// peaks at the surface instead of starting at 0 C.
TwSounding tw_triangles(const std::initializer_list<float> energies,
                        const bool surface_peak = false) {
    constexpr float peak = 4.0f;
    TwSounding snd;
    float z = 0.0f;
    bool first = true;
    for (const float energy : energies) {
        const float dtw = std::copysign(peak, energy);
        // area = peak * depth / 2 for a full or a half triangle
        const float depth = 2.0f * std::fabs(energy) / energy_of_area(peak);
        if (first && surface_peak) {
            snd.add(z, dtw);
        } else {
            if (first) snd.add(z, 0.0f);
            snd.add(z + depth / 2.0f, dtw);
        }
        z += depth;
        snd.add(z, 0.0f);
        first = false;
    }
    return snd;
}

void check_energies(const sharp::BourgouinEnergy& energy,
                    const float melting_total, const float melting_aloft,
                    const float refreezing) {
    CHECK(energy.melting_energy_total == doctest::Approx(melting_total));
    CHECK(energy.melting_energy_aloft == doctest::Approx(melting_aloft));
    CHECK(energy.refreezing_energy == doctest::Approx(refreezing));
}

void check_energies_missing(const sharp::BourgouinEnergy& energy) {
    CHECK(energy.melting_energy_total == sharp::MISSING);
    CHECK(energy.melting_energy_aloft == sharp::MISSING);
    CHECK(energy.refreezing_energy == sharp::MISSING);
}
}  // namespace

TEST_CASE("Testing bourgouin_energy on the paper's Fig. 1 profiles") {
    // Linear Tw between levels, with every 0 C crossing on a level, so each
    // layer is a triangle or a trapezoid. Heights are those of Fig. 1.
    {
        INFO("(a) all below 0 C: no energy");
        TwSounding snd;
        snd.add(0.0f, -1.0f);
        snd.add(1000.0f, -3.0f);
        snd.add(3000.0f, -10.0f);
        check_energies(snd.energy(), 0.0f, 0.0f, 0.0f);
    }
    {
        INFO("(b) surface melting below 300 m: ME_total only");
        TwSounding snd;
        snd.add(0.0f, 2.0f);
        snd.add(300.0f, 0.0f);
        snd.add(3000.0f, -10.0f);
        check_energies(snd.energy(), energy_of_area(0.5f * 2.0f * 300.0f), 0.0f,
                       0.0f);
    }
    {
        INFO("(c) melting at 1000-2500 m over surface refreezing");
        TwSounding snd;
        snd.add(0.0f, -3.0f);
        snd.add(1000.0f, 0.0f);
        snd.add(1750.0f, 4.0f);
        snd.add(2500.0f, 0.0f);
        snd.add(4000.0f, -10.0f);
        const float melting = energy_of_area(0.5f * 4.0f * 1500.0f);
        check_energies(snd.energy(), melting, melting,
                       energy_of_area(0.5f * 3.0f * 1000.0f));
    }
    {
        INFO("(d) as (c), with surface melting below 180 m");
        TwSounding snd;
        snd.add(0.0f, 1.0f);
        snd.add(180.0f, 0.0f);
        snd.add(490.0f, -2.0f);
        snd.add(800.0f, 0.0f);
        snd.add(1300.0f, 3.0f);
        snd.add(1800.0f, 0.0f);
        snd.add(3000.0f, -8.0f);
        const float surface = 0.5f * 1.0f * 180.0f;
        const float aloft = 0.5f * 3.0f * 1000.0f;
        check_energies(snd.energy(), energy_of_area(surface + aloft),
                       energy_of_area(aloft),
                       energy_of_area(0.5f * 2.0f * 620.0f));
    }
}

TEST_CASE("Testing bourgouin_energy with a crossing between levels") {
    // A crossing between levels splits the segment at the interpolated
    // height: -2 K at 0 m to 6 K at 1000 m crosses at 250 m.
    TwSounding snd;
    snd.add(0.0f, -2.0f);
    snd.add(1000.0f, 6.0f);
    snd.add(2000.0f, 2.0f);
    const float melting = 0.5f * 6.0f * 750.0f + 0.5f * (6.0f + 2.0f) * 1000.0f;
    check_energies(snd.energy(), energy_of_area(melting),
                   energy_of_area(melting),
                   energy_of_area(0.5f * 2.0f * 250.0f));
}

TEST_CASE("Testing bourgouin_energy with several warm layers (A2)") {
    // RE comes from the near-surface cold layer only, and ME_aloft adds up
    // every warm layer above it. The cold 80 J/kg layer never counts.
    check_energies(tw_triangles({-100.0f, 30.0f, -80.0f, 20.0f}).energy(),
                   50.0f, 50.0f, 100.0f);

    // A warm surface layer below them adds to ME_total only (Fig. 1d).
    check_energies(
        tw_triangles({40.0f, -100.0f, 30.0f, -80.0f, 20.0f}, true).energy(),
        90.0f, 50.0f, 100.0f);

    // With no warm layer above it, a cold layer is not a near-surface cold
    // layer, so RE is 0 however cold the column gets aloft.
    check_energies(tw_triangles({40.0f, -100.0f}, true).energy(), 40.0f, 0.0f,
                   0.0f);
}

TEST_CASE("Testing bourgouin_energy under a warm surface layer (A3)") {
    // The documented example: warm 150 at the surface, then cold 51 and
    // warm 10 J/kg.
    const TwSounding snd = tw_triangles({150.0f, -51.0f, 10.0f}, true);
    const sharp::BourgouinEnergy energy = snd.energy();
    check_energies(energy, 160.0f, 10.0f, 51.0f);

    // As in the paper, the warm surface does not suppress ice pellets:
    // rain 100 %, and 2.3 * 51 - 42 ln(11) + 3 = 19.5884 % ice pellets.
    REQUIRE(snd.wetbulb[0] > sharp::ZEROCNK);
    const auto probs = sharp::modified_bourgouin(energy, 1.0f, snd.wetbulb[0]);
    CHECK(probs.rain == 1.0f);
    CHECK(probs.freezing_rain == 0.0f);
    CHECK(probs.ice_pellets == doctest::Approx(0.195884).epsilon(1e-4));
}

TEST_CASE("Testing bourgouin_energy with a 0 C isothermal layer") {
    {
        INFO("an isothermal 0 C column has zero energy, not MISSING");
        TwSounding snd;
        snd.add(0.0f, 0.0f);
        snd.add(1000.0f, 0.0f);
        snd.add(2000.0f, 0.0f);
        check_energies(snd.energy(), 0.0f, 0.0f, 0.0f);
    }
    {
        INFO("0 C from 500 to 1500 m between cold and warm");
        TwSounding snd;
        snd.add(0.0f, -2.0f);
        snd.add(500.0f, 0.0f);
        snd.add(1500.0f, 0.0f);
        snd.add(2000.0f, 2.0f);
        const float triangle = energy_of_area(0.5f * 2.0f * 500.0f);
        check_energies(snd.energy(), triangle, triangle, triangle);
    }
    {
        // A 0 C stretch inside a cold layer does not split it, so RE
        // counts the cold air on both sides of it. Splitting it would give
        // RE from the lowest triangle only.
        INFO("cold, 0 C from 500 to 1500 m, cold, then warm");
        TwSounding snd;
        snd.add(0.0f, -2.0f);
        snd.add(500.0f, 0.0f);
        snd.add(1500.0f, 0.0f);
        snd.add(2000.0f, -2.0f);
        snd.add(2500.0f, 0.0f);
        snd.add(3000.0f, 2.0f);
        const float triangle = energy_of_area(0.5f * 2.0f * 500.0f);
        check_energies(snd.energy(), triangle, triangle, 3.0f * triangle);
    }
}

TEST_CASE("Testing the bourgouin_energy pressure cap") {
    // Cold below 500 m and warm through the top level. In each case the
    // top kept segment is warm, so a lost top level shows up in ME_total.
    // Melting area through level 1: 250 K m; through level 2: 1250 K m;
    // through level 3: 3250 K m.
    constexpr std::ptrdiff_t N = 4;
    constexpr float hght[N] = {0.0f, 1000.0f, 2000.0f, 3000.0f};
    constexpr float K = sharp::ZEROCNK;
    constexpr float tw[N] = {K - 1.0f, K + 1.0f, K + 1.0f, K + 3.0f};
    const float refreezing = energy_of_area(250.0f);
    {
        INFO("every level at or above 250 hPa: the top level is kept");
        constexpr float pres[N] = {100000.0f, 90000.0f, 80000.0f, 30000.0f};
        const auto energy = sharp::bourgouin_energy(pres, hght, tw, N);
        check_energies(energy, energy_of_area(3250.0f), energy_of_area(3250.0f),
                       refreezing);
    }
    {
        INFO("a top level exactly at the cap is kept");
        constexpr float pres[N] = {100000.0f, 90000.0f, 80000.0f, 25000.0f};
        const auto energy = sharp::bourgouin_energy(pres, hght, tw, N);
        check_energies(energy, energy_of_area(3250.0f), energy_of_area(3250.0f),
                       refreezing);
    }
    constexpr float pres[N] = {100000.0f, 90000.0f, 80000.0f, 20000.0f};
    {
        INFO("the cap between levels drops the top level");
        const auto energy = sharp::bourgouin_energy(pres, hght, tw, N);
        check_energies(energy, energy_of_area(1250.0f), energy_of_area(1250.0f),
                       refreezing);
    }
    {
        INFO("pressure_min = 0 keeps every level");
        const auto energy =
            sharp::bourgouin_energy(pres, hght, tw, N, 0.0f, 0.0f);
        check_energies(energy, energy_of_area(3250.0f), energy_of_area(3250.0f),
                       refreezing);
    }
    {
        INFO("pressure_min is configurable, and a level at it is kept");
        auto energy =
            sharp::bourgouin_energy(pres, hght, tw, N, 0.0f, 80000.0f);
        check_energies(energy, energy_of_area(1250.0f), energy_of_area(1250.0f),
                       refreezing);
        energy = sharp::bourgouin_energy(pres, hght, tw, N, 0.0f, 85000.0f);
        check_energies(energy, energy_of_area(250.0f), energy_of_area(250.0f),
                       refreezing);
    }
}

TEST_CASE("Testing bourgouin_energy ignores a warm layer above the cap") {
    // Below the cap: cold to 500 m, then warm to the top kept level at
    // 2000 m (melting area 1250 K m). Above it: cold, then warm again.
    constexpr std::ptrdiff_t N = 5;
    constexpr float pres[N] = {100000.0f, 90000.0f, 80000.0f, 20000.0f,
                               15000.0f};
    constexpr float hght[N] = {0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f};
    constexpr float K = sharp::ZEROCNK;
    constexpr float tw[N] = {K - 1.0f, K + 1.0f, K + 1.0f, K - 4.0f, K + 6.0f};
    check_energies(sharp::bourgouin_energy(pres, hght, tw, N),
                   energy_of_area(1250.0f), energy_of_area(1250.0f),
                   energy_of_area(250.0f));

    // Without the cap, the segments above add warm 100 K m (2000-2200 m)
    // and 1800 K m (3400-4000 m), with cold air between them.
    const float melting = energy_of_area(1250.0f + 100.0f + 1800.0f);
    check_energies(sharp::bourgouin_energy(pres, hght, tw, N, 0.0f, 0.0f),
                   melting, melting, energy_of_area(250.0f));
}

TEST_CASE("Testing bourgouin_energy never reads wet-bulb above the cap") {
    // The documented caller tip: fill the wet-bulb above the cap with
    // MISSING. A warm value there would add melting energy if it were read,
    // in every build, and MISSING would wreck the energies in NO_QC builds.
    constexpr std::ptrdiff_t N = 4;
    constexpr float pres[N] = {100000.0f, 90000.0f, 80000.0f, 20000.0f};
    constexpr float hght[N] = {0.0f, 1000.0f, 2000.0f, 3000.0f};
    constexpr float K = sharp::ZEROCNK;
    for (const float above_cap : {sharp::MISSING, K + 50.0f}) {
        CAPTURE(above_cap);
        const float tw[N] = {K - 1.0f, K + 1.0f, K + 1.0f, above_cap};
        check_energies(sharp::bourgouin_energy(pres, hght, tw, N),
                       energy_of_area(1250.0f), energy_of_area(1250.0f),
                       energy_of_area(250.0f));
    }
}

TEST_CASE("Testing bourgouin_energy with fewer than 2 levels under the cap") {
    // N < 2 returns before reading any element: these null arrays would
    // crash otherwise.
    check_energies_missing(
        sharp::bourgouin_energy(nullptr, nullptr, nullptr, 0));
    check_energies_missing(
        sharp::bourgouin_energy(nullptr, nullptr, nullptr, 1));

    // Valid, warm levels above the cap don't count toward the 2.
    constexpr float hght[3] = {0.0f, 1000.0f, 2000.0f};
    constexpr float K = sharp::ZEROCNK;
    constexpr float tw[3] = {K + 5.0f, K + 5.0f, K + 5.0f};
    {
        INFO("one level below the cap");
        constexpr float pres[3] = {30000.0f, 20000.0f, 15000.0f};
        check_energies_missing(sharp::bourgouin_energy(pres, hght, tw, 3));
    }
    {
        INFO("no level below the cap");
        constexpr float pres[3] = {24000.0f, 20000.0f, 15000.0f};
        check_energies_missing(sharp::bourgouin_energy(pres, hght, tw, 3));
    }
    {
        INFO("two levels below the cap");
        constexpr float pres[3] = {30000.0f, 25000.0f, 15000.0f};
        const float melting = energy_of_area(5000.0f);
        check_energies(sharp::bourgouin_energy(pres, hght, tw, 3), melting,
                       0.0f, 0.0f);
    }
}

#ifndef NO_QC
TEST_CASE("Testing bourgouin_energy with MISSING and NaN levels") {
    constexpr float MISSING = sharp::MISSING;
    constexpr float nanval = std::numeric_limits<float>::quiet_NaN();
    constexpr float K = sharp::ZEROCNK;
    constexpr float pres[4] = {100000.0f, 90000.0f, 80000.0f, 70000.0f};
    {
        // Each profile reduces to -2 K at 0 m and 2 K at 1000 m.
        INFO("missing levels are bridged");
        constexpr float hght[4] = {0.0f, 250.0f, 750.0f, 1000.0f};
        const float triangle = energy_of_area(0.5f * 2.0f * 500.0f);
        constexpr float tw[4] = {K - 2.0f, MISSING, nanval, K + 2.0f};
        check_energies(sharp::bourgouin_energy(pres, hght, tw, 4), triangle,
                       triangle, triangle);
    }
    {
        INFO("one valid level gives MISSING, not zeros");
        constexpr float hght[4] = {0.0f, 250.0f, 750.0f, 1000.0f};
        constexpr float tw[4] = {MISSING, K + 2.0f, nanval, MISSING};
        check_energies_missing(sharp::bourgouin_energy(pres, hght, tw, 4));
    }
    {
        INFO("no valid level gives MISSING");
        constexpr float hght[4] = {0.0f, 250.0f, 750.0f, 1000.0f};
        constexpr float tw[4] = {MISSING, nanval, nanval, MISSING};
        check_energies_missing(sharp::bourgouin_energy(pres, hght, tw, 4));
    }
}
#endif

TEST_CASE("Testing bourgouin_energy with min_energy") {
    {
        INFO("a weak warm layer aloft gives ME_aloft 0");
        const TwSounding snd = tw_triangles({-100.0f, 1.99f, -80.0f});
        check_energies(snd.energy(2.0f), 1.99f, 0.0f, 0.0f);
        check_energies(snd.energy(), 1.99f, 1.99f, 100.0f);

        // The same at the top of the profile.
        const TwSounding top = tw_triangles({-100.0f, 1.99f});
        check_energies(top.energy(2.0f), 1.99f, 0.0f, 0.0f);
        check_energies(top.energy(), 1.99f, 1.99f, 100.0f);
    }
    {
        INFO("a weak warm layer merges two cold layers into one");
        const TwSounding snd = tw_triangles({-100.0f, 1.99f, -80.0f, 20.0f});
        check_energies(snd.energy(2.0f), 21.99f, 20.0f, 180.0f);
        check_energies(snd.energy(), 21.99f, 21.99f, 100.0f);
    }
    {
        // ME_total counts the merged warm layer, and ME_aloft does not.
        INFO("ME_total and ME_aloft differ over a cold surface");
        const TwSounding snd = tw_triangles({-100.0f, 1.99f, -80.0f, 2.01f});
        check_energies(snd.energy(2.0f), 4.0f, 2.01f, 180.0f);
        check_energies(snd.energy(), 4.0f, 4.0f, 100.0f);
    }
}

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

namespace {
// A sounding for the full-column modified_bourgouin. Pressure falls 10 Pa
// per meter from 1000 hPa, so every level below 7.5 km is at or above
// 250 hPa.
struct ColumnSounding {
    std::vector<float> pressure;
    std::vector<float> height;
    std::vector<float> temperature;
    std::vector<float> dewpoint;
    std::vector<float> wetbulb;

    // Adds a level at height z (m) with the temperature, dewpoint, and
    // wet-bulb temperature in C.
    void add(const float z, const float tmpc, const float dwpc,
             const float wbc) {
        pressure.push_back(100000.0f - 10.0f * z);
        height.push_back(z);
        temperature.push_back(sharp::ZEROCNK + tmpc);
        dewpoint.push_back(sharp::ZEROCNK + dwpc);
        wetbulb.push_back(sharp::ZEROCNK + wbc);
    }

    // Adds a saturated level, where all three temperatures are equal.
    void add_saturated(const float z, const float tmpc) {
        add(z, tmpc, tmpc, tmpc);
    }

    std::ptrdiff_t size() const {
        return static_cast<std::ptrdiff_t>(height.size());
    }

    sharp::PrecipTypeProbabilities probs() const {
        return sharp::modified_bourgouin(pressure.data(), height.data(),
                                         temperature.data(), dewpoint.data(),
                                         wetbulb.data(), size());
    }

    sharp::PrecipTypeProbabilities probs(const float min_depth,
                                         const float min_energy,
                                         const float pressure_min) const {
        return sharp::modified_bourgouin(
            pressure.data(), height.data(), temperature.data(), dewpoint.data(),
            wetbulb.data(), size(), min_depth, min_energy, pressure_min);
    }

    // The steps of the full-column function, called one at a time.
    sharp::PrecipTypeProbabilities composed(
        const float min_depth = 0.0f, const float min_energy = 0.0f,
        const float pressure_min = sharp::BOURGOUIN_PRESSURE_MIN) const {
        const sharp::HeightLayer layer = sharp::precipitation_generation_layer(
            pressure.data(), height.data(), temperature.data(), dewpoint.data(),
            size(), min_depth);
        if (layer.bottom == sharp::MISSING) {
            return sharp::PrecipTypeProbabilities{};
        }
        const float prob_ice = sharp::probability_of_ice(
            sharp::layer_min(layer, height.data(), temperature.data(), size()));
        const sharp::BourgouinEnergy energy = sharp::bourgouin_energy(
            pressure.data(), height.data(), wetbulb.data(), size(), min_energy,
            pressure_min);
        return sharp::modified_bourgouin(energy, prob_ice, wetbulb[0]);
    }
};

// A ColumnSounding from a RelhSounding and a wet-bulb profile (K).
ColumnSounding with_wetbulb(const RelhSounding& snd,
                            const std::vector<float>& wetbulb) {
    REQUIRE(wetbulb.size() == snd.height.size());
    return {snd.pressure, snd.height, snd.temperature, snd.dewpoint, wetbulb};
}

void check_probs(const sharp::PrecipTypeProbabilities& probs,
                 const sharp::PrecipTypeProbabilities& expected) {
    CHECK(probs.rain == doctest::Approx(expected.rain));
    CHECK(probs.snow == doctest::Approx(expected.snow));
    CHECK(probs.freezing_rain == doctest::Approx(expected.freezing_rain));
    CHECK(probs.ice_pellets == doctest::Approx(expected.ice_pellets));
}

void check_same_probs(const sharp::PrecipTypeProbabilities& probs,
                      const sharp::PrecipTypeProbabilities& expected) {
    CHECK(probs.rain == expected.rain);
    CHECK(probs.snow == expected.snow);
    CHECK(probs.freezing_rain == expected.freezing_rain);
    CHECK(probs.ice_pellets == expected.ice_pellets);
}

bool same_probs(const sharp::PrecipTypeProbabilities& a,
                const sharp::PrecipTypeProbabilities& b) {
    return (a.rain == b.rain) && (a.snow == b.snow) &&
           (a.freezing_rain == b.freezing_rain) &&
           (a.ice_pellets == b.ice_pellets);
}

// ProbIce at -10 C, where the Eq. 2b polynomial gives 51 %.
constexpr float PROB_ICE_M10 = 0.51f;

// Saturated soundings, where the whole column is one generation layer and
// ProbIce comes from its coldest level.
ColumnSounding all_snow_sounding() {
    // All below 0 C, with -16 C at the top: ProbIce 1 and no melting.
    ColumnSounding snd;
    snd.add_saturated(0.0f, -2.0f);
    snd.add_saturated(1500.0f, -8.0f);
    snd.add_saturated(3000.0f, -16.0f);
    return snd;
}

ColumnSounding surface_melting_sounding() {
    // Fig. 1b: melting below 300 m, then cold to -16 C.
    ColumnSounding snd;
    snd.add_saturated(0.0f, 2.0f);
    snd.add_saturated(300.0f, 0.0f);
    snd.add_saturated(1500.0f, -8.0f);
    snd.add_saturated(3000.0f, -16.0f);
    return snd;
}

ColumnSounding warm_nose_sounding() {
    // Fig. 1c: refreezing below 1400 m, a 2 C warm nose at 1400-1800 m,
    // then cold to -10 C at the top.
    ColumnSounding snd;
    snd.add_saturated(0.0f, -8.0f);
    snd.add_saturated(1400.0f, 0.0f);
    snd.add_saturated(1600.0f, 2.0f);
    snd.add_saturated(1800.0f, 0.0f);
    snd.add_saturated(3000.0f, -10.0f);
    return snd;
}

ColumnSounding warm_above_generation_layer_sounding() {
    // A saturated generation layer from the surface to about 2.2 km, down
    // to -10 C. Above it, a dry layer is warmer than 0 C from about 2.4
    // to 3.5 km and reaches -16 C at the top.
    ColumnSounding snd;
    snd.add_saturated(0.0f, -8.0f);
    snd.add_saturated(1000.0f, -9.0f);
    snd.add_saturated(2000.0f, -10.0f);
    snd.add(2500.0f, 4.0f, -20.0f, 2.0f);
    snd.add(3000.0f, 4.0f, -20.0f, 2.0f);
    snd.add(3500.0f, 1.0f, -25.0f, 0.0f);
    snd.add(5000.0f, -16.0f, -30.0f, -18.0f);
    return snd;
}
}  // namespace

TEST_CASE("Testing modified_bourgouin on idealized soundings") {
    // The expected probabilities come from the overload that takes
    // energies, given closed-form energies, the ProbIce of the coldest level
    // of the generation layer, and the surface wet-bulb temperature.
    {
        INFO("all snow");
        const auto snd = all_snow_sounding();
        const auto probs = snd.probs();
        check_probs(probs, sharp::modified_bourgouin({0.0f, 0.0f, 0.0f}, 1.0f,
                                                     snd.wetbulb[0]));
        CHECK(probs.snow == 1.0f);
        CHECK(probs.rain == 0.0f);
        CHECK(probs.freezing_rain == 0.0f);
        CHECK(probs.ice_pellets == 0.0f);
    }
    {
        INFO("surface melting (Fig. 1b): rain, with some snow");
        const auto snd = surface_melting_sounding();
        const float melting = energy_of_area(0.5f * 2.0f * 300.0f);
        const auto probs = snd.probs();
        check_probs(probs, sharp::modified_bourgouin({melting, 0.0f, 0.0f},
                                                     1.0f, snd.wetbulb[0]));
        // 1540 exp(-0.29 * 10.77) = 67.8 %
        CHECK(probs.snow == doctest::Approx(0.678).epsilon(1e-3));
        CHECK(probs.rain == 1.0f);
        CHECK(probs.freezing_rain == 0.0f);
        CHECK(probs.ice_pellets == 0.0f);
    }
    {
        INFO("a warm nose over surface refreezing (Fig. 1c)");
        const auto snd = warm_nose_sounding();
        const float melting = energy_of_area(0.5f * 2.0f * 400.0f);
        const float refreezing = energy_of_area(0.5f * 8.0f * 1400.0f);
        const auto probs = snd.probs();
        check_probs(probs,
                    sharp::modified_bourgouin({melting, melting, refreezing},
                                              PROB_ICE_M10, snd.wetbulb[0]));
        // Eq. 8 is clamped to 100 %, so ice pellets are ProbIce.
        CHECK(probs.ice_pellets == doctest::Approx(PROB_ICE_M10));
        CHECK(probs.rain == 0.0f);
    }
}

TEST_CASE("Testing modified_bourgouin with a warm layer above the cloud") {
    // ProbIce comes from the generation layer only: its coldest level is
    // -10 C, and the -16 C level aloft doesn't count. The melting energy
    // comes from the whole column, including the warm layer above the
    // generation layer.
    const auto snd = warm_above_generation_layer_sounding();
    const sharp::HeightLayer layer = sharp::precipitation_generation_layer(
        snd.pressure.data(), snd.height.data(), snd.temperature.data(),
        snd.dewpoint.data(), snd.size());
    REQUIRE(layer.bottom == 0.0f);
    REQUIRE(layer.top > 2000.0f);
    REQUIRE(layer.top < 2500.0f);

    // The wet-bulb crosses 0 C at 2416.67 m.
    const float cross = 2000.0f + 500.0f * 10.0f / 12.0f;
    const float refreezing = energy_of_area(0.5f * (8.0f + 9.0f) * 1000.0f +
                                            0.5f * (9.0f + 10.0f) * 1000.0f +
                                            0.5f * 10.0f * (cross - 2000.0f));
    const float melting = energy_of_area(0.5f * 2.0f * (2500.0f - cross) +
                                         2.0f * 500.0f + 0.5f * 2.0f * 500.0f);
    const auto probs = snd.probs();
    check_probs(probs, sharp::modified_bourgouin({melting, melting, refreezing},
                                                 PROB_ICE_M10, snd.wetbulb[0]));

    // Without the warm layer aloft, there would be no ice pellets, and snow
    // would be ProbIce. With ProbIce from the whole column, ice pellets
    // would be 1.
    CHECK(probs.ice_pellets == doctest::Approx(PROB_ICE_M10));
    CHECK(probs.freezing_rain == doctest::Approx(1.0f - PROB_ICE_M10));
    CHECK(probs.snow < 1e-5f);
}

TEST_CASE("Testing modified_bourgouin without a generation layer") {
    {
        // A warm cloud 800 m deep under dry air. The energies are valid,
        // but a cloud 1 km deep or less is not a generation layer.
        INFO("a drizzle cloud");
        ColumnSounding snd;
        snd.add_saturated(0.0f, 5.0f);
        snd.add_saturated(800.0f, 4.0f);
        snd.add(1000.0f, 3.0f, -15.0f, -2.0f);
        snd.add(3000.0f, -10.0f, -30.0f, -14.0f);
        check_all_missing(snd.probs());

        // The documented fallback: the energies with ProbIce 0.
        const auto energy =
            sharp::bourgouin_energy(snd.pressure.data(), snd.height.data(),
                                    snd.wetbulb.data(), snd.size());
        const auto probs =
            sharp::modified_bourgouin(energy, 0.0f, snd.wetbulb[0]);
        CHECK(probs.rain == 1.0f);
        CHECK(probs.snow == 0.0f);
        CHECK(probs.freezing_rain == 0.0f);
        CHECK(probs.ice_pellets == 0.0f);
    }
    {
        // A 1.4 km cloud above a 1.6 km dry layer at the surface, which
        // eliminates it.
        INFO("virga");
        ColumnSounding snd;
        snd.add(0.0f, -2.0f, -20.0f, -6.0f);
        snd.add(1600.0f, -8.0f, -25.0f, -11.0f);
        snd.add_saturated(1650.0f, -8.5f);
        snd.add_saturated(3000.0f, -16.0f);
        check_all_missing(snd.probs());
    }
}

TEST_CASE("Testing modified_bourgouin with no or one level") {
    // N < 2 reads no element: these null arrays would crash otherwise.
    check_all_missing(sharp::modified_bourgouin(nullptr, nullptr, nullptr,
                                                nullptr, nullptr, 0));
    check_all_missing(sharp::modified_bourgouin(nullptr, nullptr, nullptr,
                                                nullptr, nullptr, 1));
}

TEST_CASE("Testing modified_bourgouin matches its steps called by hand") {
    for (const auto& snd :
         {all_snow_sounding(), surface_melting_sounding(), warm_nose_sounding(),
          warm_above_generation_layer_sounding()}) {
        const auto probs = snd.probs();
        check_same_probs(probs, snd.composed());
        check_same_probs(probs, snd.probs(0.0f, 0.0f, 25000.0f));
    }
}

TEST_CASE("Testing modified_bourgouin passes each option to its own step") {
    {
        // A 50 m dry sliver splits a cloud at -10 C into two 600 m moist
        // layers, so there is no generation layer. min_depth = 100 m
        // absorbs it. min_energy and pressure_min don't reach the
        // generation layer.
        INFO("min_depth");
        const auto relh = relh_sounding(
            {0.0f, 300.0f, 600.0f, 625.0f, 650.0f, 950.0f, 1250.0f, 1500.0f,
             3500.0f},
            {0.9f, 0.9f, 0.75f, 0.5f, 0.75f, 0.9f, 0.75f, 0.5f, 0.5f},
            std::vector<float>(9, sharp::ZEROCNK - 10.0f));
        const auto snd =
            with_wetbulb(relh, std::vector<float>(9, sharp::ZEROCNK - 11.0f));
        check_all_missing(snd.probs());
        check_all_missing(snd.probs(0.0f, 100.0f, 25000.0f));
        check_all_missing(snd.probs(0.0f, 0.0f, 100.0f));

        const auto probs = snd.probs(100.0f, 0.0f, 25000.0f);
        check_same_probs(probs, snd.composed(100.0f));
        // No melting, so snow is ProbIce, about 51 % at -10 C.
        CHECK(probs.snow == doctest::Approx(PROB_ICE_M10).epsilon(1e-4));
    }
    {
        // Saturated, from the surface up: cold 100, warm 1.99, cold 80, and
        // warm 20 J/kg, then -16 C at the top. min_energy = 2 J/kg merges
        // the two cold layers. min_depth doesn't reach the energies.
        INFO("min_energy");
        const TwSounding tw = tw_triangles({-100.0f, 1.99f, -80.0f, 20.0f});
        ColumnSounding snd;
        for (std::ptrdiff_t k = 0; k < tw.N; ++k) {
            snd.add_saturated(tw.hght[k], tw.wetbulb[k] - sharp::ZEROCNK);
        }
        snd.add_saturated(tw.hght[tw.N - 1] + 2000.0f, -16.0f);

        const auto probs = snd.probs(0.0f, 2.0f, 25000.0f);
        check_same_probs(probs, snd.composed(0.0f, 2.0f));
        CHECK_FALSE(same_probs(probs, snd.probs()));
        check_same_probs(snd.probs(2.0f, 0.0f, 25000.0f), snd.probs());

        // RE 180 instead of 100: Eq. 7 gives 84.4 % instead of 100 %.
        CHECK(probs.freezing_rain == doctest::Approx(0.844).epsilon(1e-3));
        CHECK(snd.probs().freezing_rain == 1.0f);
    }
    {
        // Saturated and cold, with a warm layer above 250 hPa that only
        // pressure_min = 0 brings in.
        INFO("pressure_min");
        ColumnSounding snd;
        snd.add_saturated(0.0f, -2.0f);
        snd.add_saturated(3000.0f, -16.0f);
        snd.add_saturated(7000.0f, -30.0f);
        snd.add_saturated(8000.0f, 5.0f);
        snd.add_saturated(9000.0f, -40.0f);
        REQUIRE(snd.pressure[3] < sharp::BOURGOUIN_PRESSURE_MIN);

        const auto probs = snd.probs(0.0f, 0.0f, 0.0f);
        check_same_probs(probs, snd.composed(0.0f, 0.0f, 0.0f));
        CHECK_FALSE(same_probs(probs, snd.probs()));
        CHECK(snd.probs().snow == 1.0f);
        CHECK(probs.ice_pellets == 1.0f);
    }
}

#ifndef NO_QC
TEST_CASE("Testing modified_bourgouin with a leading missing wet-bulb") {
    // The surface wet-bulb is the lowest valid level's, 2 C. Liquid
    // relative humidity 0.93 everywhere makes the column a generation
    // layer, and its 280 K gives ProbIce 0, so everything is rain.
    for (const float bad :
         {sharp::MISSING, std::numeric_limits<float>::quiet_NaN()}) {
        CAPTURE(bad);
        ColumnSounding snd;
        for (const float z : {0.0f, 500.0f, 1500.0f, 2500.0f}) {
            snd.add(z, 280.0f - sharp::ZEROCNK, 279.0f - sharp::ZEROCNK, 2.0f);
        }
        snd.wetbulb[0] = bad;
        REQUIRE(snd.wetbulb[1] == 275.15f);

        const auto probs = snd.probs();
        CHECK(probs.rain == 1.0f);
        CHECK(probs.snow == 0.0f);
        CHECK(probs.freezing_rain == 0.0f);
        CHECK(probs.ice_pellets == 0.0f);

        // Taking wetbulb[0] as the surface would give MISSING.
        check_all_missing(snd.composed());

        // The surface is the lowest valid level, not one above it: 2 C at
        // 500 m gives rain, where -2 C at 1500 m would give freezing rain.
        snd.wetbulb[2] = sharp::ZEROCNK - 2.0f;
        snd.wetbulb[3] = sharp::ZEROCNK - 2.0f;
        CHECK(snd.probs().rain == 1.0f);
        CHECK(snd.probs().freezing_rain == 0.0f);

        // With no valid wet-bulb at all, everything is MISSING.
        snd.wetbulb.assign(snd.wetbulb.size(), bad);
        check_all_missing(snd.probs());
    }
}

TEST_CASE("Testing modified_bourgouin with NaN next to the generation layer") {
    // Relative humidity over ice 0.6 and 0.7, then 0.9 from 1000 m up, with
    // the 1000 m temperature NaN. The generation layer starts where the
    // walk crosses 75 % between 500 and 1500 m, at 750 m, and its bottom
    // temperature bridges the NaN: 258 K + 0.25 * 10 K = 260.5 K, or
    // -12.65 C.
    auto relh =
        relh_sounding({0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f},
                      {0.6f, 0.7f, 0.9f, 0.9f, 0.9f, 0.9f},
                      {258.0f, 258.0f, 258.0f, 268.0f, 268.0f, 268.0f});
    relh.temperature[2] = std::numeric_limits<float>::quiet_NaN();
    const auto snd = with_wetbulb(relh, std::vector<float>(6, 255.0f));

    const sharp::HeightLayer layer = relh.generation_layer();
    CHECK(layer.bottom == doctest::Approx(750.0f));
    CHECK(layer.top == 2500.0f);
    const float tmin = sharp::layer_min(layer, snd.height.data(),
                                        snd.temperature.data(), snd.size());
    CHECK(tmin == doctest::Approx(260.5f));
    CHECK(sharp::probability_of_ice(tmin) ==
          doctest::Approx(0.7287).epsilon(1e-4));

    // No melting, so snow is ProbIce and freezing rain the rest.
    const auto probs = snd.probs();
    CHECK(probs.snow == doctest::Approx(0.7287).epsilon(1e-4));
    CHECK(probs.freezing_rain == doctest::Approx(0.2713).epsilon(1e-3));
    CHECK(probs.rain == 0.0f);
    CHECK(probs.ice_pellets == 0.0f);
}
#endif

#ifndef NO_QC
namespace {
constexpr float M = sharp::MISSING;

void check_missing_wind(const sharp::WindComponents wind) {
    CHECK(wind.u == M);
    CHECK(wind.v == M);
}

void check_wind(const sharp::WindComponents wind, const float u,
                const float v) {
    CHECK(wind.u == doctest::Approx(u));
    CHECK(wind.v == doctest::Approx(v));
}

constexpr std::ptrdiff_t KN = 5;
constexpr float k_pres[KN] = {100000, 95000, 90000, 85000, 80000};
constexpr float k_hght[KN] = {0, 500, 1000, 1500, 2000};
constexpr float k_uwin[KN] = {0, 5, 10, 15, 20};
constexpr float k_vwin[KN] = {0, 2, 4, 6, 8};
constexpr float d_pres[KN] = {100000, 80000, 62000, 47000, 35000};
constexpr float d_hght[KN] = {0, 2000, 4000, 6000, 8000};
constexpr float d_uwin[KN] = {0, 10, 20, 30, 40};
constexpr float d_vwin[KN] = {0, 2, 4, 6, 8};
}  // namespace

TEST_CASE("Testing effective_bulk_wind_difference with a MISSING layer") {
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, k_hght, k_uwin, k_vwin, KN, {105000, 95000}, 90000));
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, k_hght, k_uwin, k_vwin, KN, {100000, 95000}, 70000));

    check_wind(sharp::effective_bulk_wind_difference(
                   k_pres, k_hght, k_uwin, k_vwin, KN, {100000, 95000}, 85000),
               7.5f, 3.0f);
}

TEST_CASE("Testing storm_motion_bunkers with a MISSING layer") {
    check_missing_wind(sharp::storm_motion_bunkers(
        k_pres, k_hght, k_uwin, k_vwin, KN, {0, 3000}, {0, 2000}));
    check_missing_wind(sharp::storm_motion_bunkers(
        k_pres, k_hght, k_uwin, k_vwin, KN, {0, 2000}, {0, 3000}));

    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, d_hght, d_uwin, d_vwin, KN, {0, 6000}, sharp::HeightLayer()));
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, d_hght, d_uwin, d_vwin, KN, sharp::HeightLayer(), {0, 6000}));

    check_wind(sharp::storm_motion_bunkers(k_pres, k_hght, k_uwin, k_vwin, KN,
                                           {0, 2000}, {0, 2000}),
               12.78543f, -2.96357536f);
}

TEST_CASE("Testing effective-inflow storm_motion_bunkers fallback") {
    constexpr float vwin[KN] = {0, 0, 0, 0, 0};
    sharp::Parcel mupcl;
    mupcl.eql_pressure = 80000;
    for (const bool left : {false, true}) {
        CAPTURE(left);
        const sharp::WindComponents fallback = sharp::storm_motion_bunkers(
            d_pres, d_hght, d_uwin, vwin, KN, {0, 6000}, {0, 6000}, left);
        const sharp::WindComponents motion = sharp::storm_motion_bunkers(
            d_pres, d_hght, d_uwin, vwin, KN, {105000, 95000}, mupcl, left);
        CHECK(motion.u == fallback.u);
        CHECK(motion.v == fallback.v);
        check_wind(motion, 14.0566034f, left ? 7.5f : -7.5f);
    }

    mupcl.eql_pressure = 40000;
    check_wind(sharp::storm_motion_bunkers(d_pres, d_hght, d_uwin, d_vwin, KN,
                                           {100000, 90000}, mupcl),
               11.077548f, -5.43301964f);
}

TEST_CASE("Testing mcs_motion_corfidi with a MISSING layer") {
    constexpr float hght[KN] = {0, 500, 1000, 2000, 3000};
    constexpr float pres[KN] = {100000, 95000, 90000, 80000, 70000};
    const auto check_missing_pair = [](const auto vectors) {
        check_missing_wind(vectors.first);
        check_missing_wind(vectors.second);
    };
    constexpr float hght_low[KN] = {0, 250, 500, 750, 1000};
    check_missing_pair(
        sharp::mcs_motion_corfidi(k_pres, hght_low, k_uwin, k_vwin, KN));

    const auto vectors =
        sharp::mcs_motion_corfidi(pres, hght, k_uwin, k_vwin, KN);
    check_wind(vectors.first, 9.16666603f, 3.66666651f);
    check_wind(vectors.second, 25.4044037f, 10.1617622f);
}

TEST_CASE("Testing large_hail_parameter with a MISSING layer") {
    constexpr std::ptrdiff_t N = 6;
    constexpr float hght[N] = {0, 1500, 3000, 4500, 5500, 7000};
    constexpr float pres[N] = {100000, 85000, 70000, 59000, 51000, 40000};
    constexpr float uwin[N] = {0, 6, 12, 18, 24, 30};
    constexpr float vwin[N] = {0, 2, 4, 6, 8, 10};
    sharp::Parcel mu_pcl;
    mu_pcl.cape = 3000;
    mu_pcl.eql_pressure = 55000;
    const sharp::WindComponents storm = {5, 5};
    const sharp::PressureLayer hgz = {65000, 52000};

    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, {M, M}, storm, pres, hght,
                                      uwin, vwin, N) == M);
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, {65000, 35000}, storm, pres,
                                      hght, uwin, vwin, N) == M);
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N - 1) == M);
    sharp::Parcel low_el = mu_pcl;
    low_el.eql_pressure = 90000;
    CHECK(sharp::large_hail_parameter(low_el, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N) == M);
    sharp::Parcel high_el = mu_pcl;
    high_el.eql_pressure = 30000;
    CHECK(sharp::large_hail_parameter(high_el, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N) == M);

    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin,
                                      N) == doctest::Approx(70.2362289f));
}

TEST_CASE("Testing wind parameters with MISSING wind data") {
    constexpr float none[KN] = {M, M, M, M, M};
    constexpr float u_low[KN] = {0, 10, 20, M, M};
    constexpr float v_low[KN] = {0, 2, 4, M, M};
    check_missing_wind(sharp::storm_motion_bunkers(d_pres, d_hght, none, none,
                                                   KN, {0, 6000}, {0, 6000}));
    check_missing_wind(sharp::storm_motion_bunkers(d_pres, d_hght, u_low, v_low,
                                                   KN, {0, 3000}, {0, 6000}));
    check_missing_wind(sharp::storm_motion_bunkers(d_pres, d_hght, u_low, v_low,
                                                   KN, {0, 6000}, {0, 6000}));

    const auto vectors =
        sharp::mcs_motion_corfidi(d_pres, d_hght, none, none, KN);
    check_missing_wind(vectors.first);
    check_missing_wind(vectors.second);

    constexpr std::ptrdiff_t N = 6;
    constexpr float hght[N] = {0, 1500, 3000, 4500, 5500, 7000};
    constexpr float pres[N] = {100000, 85000, 70000, 59000, 51000, 40000};
    constexpr float uwin[N] = {0, 6, 12, 18, 24, 30};
    constexpr float vwin[N] = {0, 2, 4, 6, 8, 10};
    constexpr float u_none[N] = {M, M, M, M, M, M};
    sharp::Parcel mu_pcl;
    mu_pcl.cape = 3000;
    mu_pcl.eql_pressure = 55000;
    const sharp::PressureLayer hgz = {65000, 52000};
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, {5, 5}, pres, hght,
                                      u_none, u_none, N) == M);
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, {M, M}, pres, hght,
                                      uwin, vwin, N) == M);
}
#endif

namespace {
constexpr std::ptrdiff_t EN = 6;
constexpr float e_pres[EN] = {100000, 90000, 80000, 70000, 60000, 50000};
constexpr float e_uwin[EN] = {0, 10, 12, 13, 13, 13};

sharp::WindComponents ebwd_shifted(const float shift, const float vwin[],
                                   const sharp::PressureLayer eil,
                                   const float eql_pres) {
    float hght[EN];
    for (std::ptrdiff_t k = 0; k < EN; ++k) hght[k] = 1000.0f * k + shift;
    return sharp::effective_bulk_wind_difference(e_pres, hght, e_uwin, vwin, EN,
                                                 eil, eql_pres);
}

constexpr float shifts[] = {0.0f, 1000.0f, 1234.5f, 762.3f};
}  // namespace

TEST_CASE("Testing effective_bulk_wind_difference against a known value") {
    constexpr float vwin[EN] = {0, -2, 1, 4, 6, 7};
    for (const float shift : shifts) {
        CAPTURE(shift);
        const sharp::WindComponents ebwd =
            ebwd_shifted(shift, vwin, {90000, 80000}, 60000);
        CHECK(ebwd.u == doctest::Approx(2.5f));
        CHECK(ebwd.v == doctest::Approx(4.5f));
    }
}

namespace {
constexpr std::ptrdiff_t BN = 33;

struct BunkersSounding {
    float pres[BN];
    float hght[BN];
    float uwin[BN];
    float vwin[BN];

    explicit BunkersSounding(const float elevation) {
        for (std::ptrdiff_t k = 0; k < BN; ++k) {
            const float z = 500.0f * k;
            hght[k] = elevation + z;
            pres[k] = 100000.0f * std::exp(-z / 8000.0f);
            uwin[k] = 30.0f * (1.0f - std::exp(-z / 4000.0f));
            vwin[k] = 10.0f * std::sin(z / 3000.0f);
        }
    }

    float pres_at(const float z_agl) const {
        return pres[static_cast<std::ptrdiff_t>(z_agl / 500.0f)];
    }

    sharp::WindComponents effective(const float base_agl, const float el_agl,
                                    const bool left) const {
        sharp::Parcel mupcl;
        mupcl.eql_pressure = pres_at(el_agl);
        const sharp::PressureLayer eil = {pres_at(base_agl),
                                          pres_at(base_agl + 1000.0f)};
        return sharp::storm_motion_bunkers(pres, hght, uwin, vwin, BN, eil,
                                           mupcl, left);
    }

    sharp::WindComponents classic(const sharp::HeightLayer mean_wind_agl,
                                  const bool left, const bool weighted) const {
        return sharp::storm_motion_bunkers(pres, hght, uwin, vwin, BN,
                                           mean_wind_agl, {0, 6000}, left,
                                           weighted);
    }
};

struct BunkersCase {
    float base_agl;
    float el_agl;
    float mw_top_agl;  // sharp::MISSING for the 0-6 km fallback
    float u;
    float v;
};

void check_bunkers(const BunkersCase c) {
    CAPTURE(c.base_agl);
    CAPTURE(c.el_agl);
    for (const float elevation : {0.0f, 1000.0f, 762.3f}) {
        CAPTURE(elevation);
        const BunkersSounding snd(elevation);
        for (const bool left : {false, true}) {
            CAPTURE(left);
            const sharp::WindComponents motion =
                snd.effective(c.base_agl, c.el_agl, left);
            const sharp::WindComponents expected =
                (c.mw_top_agl != sharp::MISSING)
                    ? snd.classic({c.base_agl, c.mw_top_agl}, left, true)
                    : snd.classic({0, 6000}, left, false);
            CHECK(motion.u == doctest::Approx(expected.u).epsilon(1e-6));
            CHECK(motion.v == doctest::Approx(expected.v).epsilon(1e-6));
            if (!left) {
                CHECK(motion.u == doctest::Approx(c.u));
                CHECK(motion.v == doctest::Approx(c.v));
            }
        }
    }
}
}  // namespace

TEST_CASE("Testing the effective-inflow storm_motion_bunkers mean wind layer") {
    for (const BunkersCase c : {
             BunkersCase{0, 10000, 6500, 14.8269749f, -1.06337833f},
             BunkersCase{1000, 12000, 7800, 18.9706116f, 0.523952484f},
         }) {
        check_bunkers(c);
    }
}

TEST_CASE("Testing the effective-inflow storm_motion_bunkers 3 km fallback") {
    for (const BunkersCase c : {
             BunkersCase{2000, 7000, sharp::MISSING, 15.8496647f,
                         -0.509417534f},
             BunkersCase{6000, 14000, 9100, 27.9163494f, -0.846437931f},
         }) {
        check_bunkers(c);
    }
}

// ===========================================================================
// Precipitation type: the spectral bin classifier (Reeves et al. 2016)
// ===========================================================================

// ---------------------------------------------------------------------------
// Result types and the drop-size distribution
// ---------------------------------------------------------------------------

static_assert(std::is_trivially_copyable_v<sharp::SpectralBinDSD>);
static_assert(std::is_trivially_copyable_v<sharp::SpectralBinResult>);

TEST_CASE("Testing the PrecipType encoding") {
    using sharp::PrecipType;
    CHECK(static_cast<int>(PrecipType::missing) == -9999);
    CHECK(static_cast<float>(PrecipType::missing) == sharp::MISSING);
    CHECK(static_cast<int>(PrecipType::rain) == 1);
    CHECK(static_cast<int>(PrecipType::snow) == 2);
    CHECK(static_cast<int>(PrecipType::rain_snow) == 3);
    CHECK(static_cast<int>(PrecipType::freezing_rain) == 4);
    CHECK(static_cast<int>(PrecipType::ice_pellets) == 5);
    CHECK(static_cast<int>(PrecipType::freezing_rain_ice_pellets) == 6);
    CHECK(static_cast<int>(PrecipType::rain_ice_pellets) == 7);
}

TEST_CASE("Testing SpectralBinResult defaults to missing") {
    const sharp::SpectralBinResult result;
    CHECK(result.precip_type == sharp::PrecipType::missing);
    CHECK(result.liquid_fraction == sharp::MISSING);
    CHECK(result.supercooled_liquid_height == sharp::MISSING);
}

TEST_CASE("Testing the spectral bin classifier constants") {
    CHECK(sharp::SBC_MAX_BINS == 64);
    CHECK(sharp::SBC_ICE_NUCLEATION_TEMPERATURE == 267.15f);
}

namespace {
struct DiameterConstants {
    std::vector<double> mass;
    std::vector<double> rain_fall_speed;
    std::vector<double> pellet_fall_speed;
    std::vector<double> foote_du_toit_fall_speed;
    std::vector<double> liquid_capacitance_factor;
    std::vector<double> liquid_length_factor;
};

struct RimeConstants {
    std::vector<double> snow_diameter;
    std::vector<double> snow_density;
    std::vector<double> snow_aa;
    std::vector<double> snow_bb;
    std::vector<double> liquid_aa;
    std::vector<double> liquid_bb;
};

// A drop-size distribution and its per-bin constants from the reference
// formulas (sbc_alg_2023Aug31.py), evaluated in float64 numpy on the float32
// diameters, at rime factors 1 and 5.
struct ExpectedDSD {
    std::vector<float> diameter;
    std::vector<float> concentration;
    DiameterConstants constants;
    RimeConstants rime_1;
    RimeConstants rime_5;
};

// The Python default (run_sbc.py, deld = 0.7)
const ExpectedDSD PYTHON_DSD{
    {0.05f, 0.75f, 1.45f, 2.15f},
    {55.1843f, 146.647f, 11.6891f, 3.60886f},
    {
        // mass
        {6.54498499e-08, 0.000220893233, 0.00159625647, 0.00520372167},
        // rain_fall_speed
        {0.142122154, 3.09237202, 5.27257807, 6.82459228},
        // pellet_fall_speed
        {0.30790029, 1.43337744, 2.51680444, 3.55818123},
        // foote_du_toit_fall_speed
        {0.0527470786, 3.04237812, 5.2708923, 6.8547722},
        // liquid_capacitance_factor
        {0.39782573, 0.397916111, 0.40139199, 0.407350743},
        // liquid_length_factor
        {1.00423449, 1.00377655, 0.986744087, 0.960031483},
    },
    {
        // snow_diameter
        {0.0629960534, 1.51198988, 3.91461928, 6.91105275},
        // snow_density
        {0.5, 0.121583861, 0.0505778374, 0.0299474772},
        // snow_aa
        {0.957678002, 2.54338381, 3.40709147, 4.05742466},
        // snow_bb
        {-0.0211609991, 0.771691903, 1.20354574, 1.52871233},
        // liquid_aa
        {0.892032466, 2.05037221, 2.5108705, 2.83400081},
        // liquid_bb
        {-0.053983767, 0.525186103, 0.755435252, 0.917000406},
    },
    {
        // snow_diameter
        {0.0629960534, 0.944940787, 1.82688558, 3.19182157},
        // snow_density
        {0.5, 0.5, 0.5, 0.305257549},
        // snow_aa
        {0.560053493, 1.28730529, 1.57642446, 1.8713215},
        // snow_bb
        {-0.219973254, 0.143652647, 0.288212228, 0.435660751},
        // liquid_aa
        {0.521663751, 1.19906494, 1.46836598, 1.65733373},
        // liquid_bb
        {-0.239168125, 0.0995324698, 0.234182989, 0.328666864},
    },
};

// The C++ MRMS code, version 2.0.3
const ExpectedDSD CXX_2_0_3_DSD{
    {0.05f, 0.65f, 1.25f, 1.85f},
    {55.1843f, 206.606f, 25.4924f, 3.60886f},
    {
        // mass
        {6.54498499e-08, 0.000143793298, 0.00102265386, 0.00331523123},
        // rain_fall_speed
        {0.142122154, 2.72153178, 4.71971152, 6.22782749},
        // pellet_fall_speed
        {0.30790029, 1.27516945, 2.21154466, 3.11702582},
        // foote_du_toit_fall_speed
        {0.0527470786, 2.66460368, 4.70504687, 6.24743003},
        // liquid_capacitance_factor
        {0.39782573, 0.397670598, 0.400111625, 0.404540105},
        // liquid_length_factor
        {1.00423449, 1.00502234, 0.992888514, 0.972254874},
    },
    {
        // snow_diameter
        {0.0629960534, 1.22989896, 3.15992395, 5.56370756},
        // snow_density
        {0.5, 0.147082291, 0.0616194969, 0.0365758243},
        // snow_aa
        {0.957678002, 2.386989, 3.19005249, 3.79582199},
        // snow_bb
        {-0.0211609991, 0.693494498, 1.09502624, 1.397911},
        // liquid_aa
        {0.892032466, 1.96215169, 2.39891147, 2.70608431},
        // liquid_bb
        {-0.053983767, 0.481075843, 0.699455737, 0.853042157},
    },
    {
        // snow_diameter
        {0.0629960534, 0.818948652, 1.57490131, 2.56955958},
        // snow_density
        {0.5, 0.5, 0.5, 0.372820936},
        // snow_aa
        {0.560053493, 1.23191694, 1.50613212, 1.75066795},
        // snow_bb
        {-0.219973254, 0.11595847, 0.253066061, 0.375333975},
        // liquid_aa
        {0.521663751, 1.14747327, 1.40289194, 1.58252771},
        // liquid_bb
        {-0.239168125, 0.0737366333, 0.20144597, 0.291263854},
    },
};

// The Python construction with deld = 0.1
const ExpectedDSD DELD_0_1_DSD{
    {0.05f, 0.15f, 0.25f, 0.35f, 0.45f, 0.55f, 0.65f, 0.75f, 0.85f, 0.95f,
     1.05f, 1.15f, 1.25f, 1.35f, 1.45f, 1.55f, 1.65f, 1.75f, 1.85f},
    {55.1843f, 66.0695f, 130.272f, 154.556f, 203.649f, 171.814f, 206.606f,
     146.647f, 94.9404f, 79.4013f, 61.0083f, 35.6567f, 25.4924f, 16.2522f,
     11.6891f, 7.49152f, 3.60886f, 3.60886f, 3.60886f},
    {
        // mass
        {6.54498499e-08, 1.76714608e-06, 8.18123087e-06, 2.24492964e-05,
         4.77129346e-05, 8.7113752e-05, 0.000143793298, 0.000220893233,
         0.000321555125, 0.000448920483, 0.00060613095, 0.000796328238,
         0.00102265386, 0.00128824941, 0.00159625647, 0.00194981621,
         0.00235207105, 0.00280616219, 0.00331523123},
        // rain_fall_speed
        {0.142122154, 0.616476787, 1.0724364, 1.51046562, 1.93102338,
         2.33456302, 2.72153178, 3.09237202, 3.44751975, 3.78740533, 4.11245389,
         4.42308483, 4.71971152, 5.00274186, 5.27257807, 5.52961642, 5.77424839,
         6.0068589, 6.22782749},
        // pellet_fall_speed
        {0.30790029, 0.471257251, 0.633756027, 0.795396634, 0.956179074,
         1.11610339, 1.27516945, 1.43337744, 1.59072725, 1.74721881, 1.9028522,
         2.05762751, 2.21154466, 2.36460363, 2.51680444, 2.6681469, 2.81863138,
         2.96825769, 3.11702582},
        // foote_du_toit_fall_speed
        {0.0527470786, 0.530851053, 0.991384375, 1.4346867, 1.86109763,
         2.27095687, 2.66460368, 3.04237812, 3.40461956, 3.75166738, 4.08386142,
         4.40154145, 4.70504687, 4.99471729, 5.2708923, 5.5339112, 5.78411422,
         6.02184063, 6.24743003},
        // liquid_capacitance_factor
        {0.39782573, 0.397586652, 0.397438505, 0.397376825, 0.397397473,
         0.397496596, 0.397670598, 0.397916111, 0.398229973, 0.398609202,
         0.399050981, 0.399552636, 0.400111625, 0.400725518, 0.40139199,
         0.402108804, 0.402873807, 0.403684914, 0.404540105},
        // liquid_length_factor
        {1.00423449, 1.00544962, 1.00620533, 1.00652059, 1.00641501, 1.00590875,
         1.00502234, 1.00377655, 1.0021923, 1.00029056, 0.998092224,
         0.995618031, 0.992888514, 0.989923906, 0.986744087, 0.983368528,
         0.979816228, 0.976105691, 0.972254874},
    },
    {
        // snow_diameter
        {0.0629960534, 0.188988165, 0.314980262, 0.503414429, 0.723470993,
         0.966448813, 1.22989896, 1.51198988, 1.81128591, 2.12662364,
         2.45703652, 2.80170531, 3.15992395, 3.53107646, 3.91461928, 4.31006782,
         4.71698836, 5.13498745, 5.56370756},
        // snow_density
        {0.5, 0.5, 0.5, 0.335154106, 0.239901922, 0.183689823, 0.147082291,
         0.121583861, 0.10293332, 0.0887747044, 0.0777070632, 0.0688488199,
         0.0616194969, 0.0556223888, 0.0505778374, 0.046283443, 0.042589354,
         0.0393824462, 0.0365758243},
        // snow_aa
        {0.957678002, 1.34231604, 1.57049403, 1.81393769, 2.02780397,
         2.21653969, 2.386989, 2.54338381, 2.68855269, 2.82449493, 2.95268452,
         3.07424431, 3.19005249, 3.30081106, 3.40709147, 3.50936596, 3.60803045,
         3.70342013, 3.79582199},
        // snow_bb
        {-0.0211609991, 0.171158021, 0.285247017, 0.406968843, 0.513901986,
         0.608269847, 0.693494498, 0.771691903, 0.844276345, 0.912247467,
         0.976342258, 1.03712216, 1.09502624, 1.15040553, 1.20354574,
         1.25468298, 1.30401523, 1.35171006, 1.397911},
        // liquid_aa
        {0.892032466, 1.25030489, 1.46284207, 1.62221142, 1.7524724, 1.86395467,
         1.96215169, 2.05037221, 2.1307801, 2.20487648, 2.27375002, 2.33821806,
         2.39891147, 2.45632856, 2.5108705, 2.5628655, 2.61258608, 2.66026094,
         2.70608431},
        // liquid_bb
        {-0.053983767, 0.125152445, 0.231421033, 0.31110571, 0.3762362,
         0.431977337, 0.481075843, 0.525186103, 0.565390048, 0.602438238,
         0.636875011, 0.669109032, 0.699455737, 0.728164281, 0.755435252,
         0.781432752, 0.806293038, 0.83013047, 0.853042157},
    },
    {
        // snow_diameter
        {0.0629960534, 0.188988165, 0.314980262, 0.44097236, 0.566964457,
         0.692956592, 0.818948652, 0.944940787, 1.07093292, 1.19692498,
         1.32291704, 1.44890918, 1.57490131, 1.70089345, 1.82688558, 1.9905748,
         2.17850822, 2.37155819, 2.56955958},
        // snow_density
        {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
         0.5, 0.471771638, 0.434117428, 0.401429105, 0.372820936},
        // snow_aa
        {0.560053493, 0.784991184, 0.918430483, 1.01848891, 1.10027194,
         1.17026495, 1.23191694, 1.28730529, 1.33778857, 1.38430923, 1.42755078,
         1.46802638, 1.50613212, 1.54218086, 1.57642446, 1.6185518, 1.66405677,
         1.70805136, 1.75066795},
        // snow_bb
        {-0.219973254, -0.107504408, -0.0407847586, 0.00924445397, 0.0501359682,
         0.0851324729, 0.11595847, 0.143652647, 0.168894284, 0.192154613,
         0.21377539, 0.234013191, 0.253066061, 0.271090428, 0.288212228,
         0.309275899, 0.332028384, 0.354025681, 0.375333975},
        // liquid_aa
        {0.521663751, 0.731182736, 0.855475229, 0.948674993, 1.02485208,
         1.09004731, 1.14747327, 1.19906494, 1.24608776, 1.28941959, 1.32969708,
         1.36739822, 1.40289194, 1.43646966, 1.46836598, 1.49877284, 1.52784961,
         1.55573004, 1.58252771},
        // liquid_bb
        {-0.239168125, -0.134408632, -0.0722623853, -0.0256625033, 0.0124260381,
         0.0450236527, 0.0737366333, 0.0995324698, 0.12304388, 0.144709793,
         0.16484854, 0.183699109, 0.20144597, 0.218234829, 0.234182989,
         0.249386419, 0.263924803, 0.277865018, 0.291263854},
    },
};

// Drops small enough for the 0.01 m/s floor and a negative Foote-du Toit
// fall speed, and large enough for the non-decreasing fix
const ExpectedDSD EDGE_DSD{
    {0.01f, 0.02f, 0.05f, 5.95f, 7.0f, 8.0f, 12.0f},
    {10.0f, 20.0f, 30.0f, 40.0f, 30.0f, 20.0f, 10.0f},
    {
        // mass
        {5.2359874e-10, 4.18878992e-09, 6.54498499e-08, 0.110293388, 0.17959438,
         0.268082573, 0.904778684},
        // rain_fall_speed
        {0.01, 0.01, 0.142122154, 9.17834173, 9.17834173, 9.17834173, 9.634028},
        // pellet_fall_speed
        {0.242317221, 0.25872586, 0.30790029, 8.47763546, 9.61844772,
         10.6169732, 13.7529075},
        // foote_du_toit_fall_speed
        {-0.143490345, -0.0941611494, 0.0527470786, 9.23763988, 9.6448, 10.6102,
         26.9558},
        // liquid_capacitance_factor
        {0.39794789, 0.397915896, 0.39782573, 0.456563843, 0.469568003,
         0.482577119, 1.11613785},
        // liquid_length_factor
        {1.00361571, 1.00377763, 1.00423449, 0.827499954, 0.8107243,
         0.799020908, 1.40744565},
    },
    {
        // snow_diameter
        {0.0125992102, 0.0251984204, 0.0629960534, 30.0236469, 37.9587545,
         46.0250568, 82.6216798},
        // snow_density
        {0.5, 0.5, 0.5, 0.00773034198, 0.00622722246, 0.00521361823,
         0.00303990013},
        // snow_aa
        {0.583986389, 0.722635474, 0.957678002, 6.37241879, 6.84866229,
         7.2664809, 8.69796068},
        // snow_bb
        {-0.208006805, -0.138682263, -0.0211609991, 2.68620939, 2.92433115,
         3.13324045, 3.84898034},
        // liquid_aa
        {0.543956129, 0.673101296, 0.892032466, 3.87494391, 4.07340265,
         4.24404715, 4.80727442},
        // liquid_bb
        {-0.228021935, -0.163449352, -0.053983767, 1.43747196, 1.53670133,
         1.62202358, 1.90363721},
    },
    {
        // snow_diameter
        {0.0125992102, 0.0251984204, 0.0629960534, 13.8662122, 17.5309864,
         21.2563519, 38.1582473},
        // snow_density
        {0.5, 0.5, 0.5, 0.0787961279, 0.0634746844, 0.0531429179, 0.0309859977},
        // snow_aa
        {0.341517312, 0.422599789, 0.560053493, 2.93901805, 3.15866592,
         3.3513677, 4.01157932},
        // snow_bb
        {-0.329241344, -0.288700105, -0.219973254, 0.969509024, 1.07933296,
         1.17568385, 1.50578966},
        // liquid_aa
        {0.318107474, 0.393632026, 0.521663751, 2.26608095, 2.38214032,
         2.48193383, 2.81131113},
        // liquid_bb
        {-0.340946263, -0.303183987, -0.239168125, 0.633040474, 0.691070161,
         0.740966915, 0.905655567},
    },
};

// float32 stays within 2e-6 of the float64 reference.
constexpr double SBC_RTOL = 1e-5;

void check_bins(const char* name,
                const std::array<float, sharp::SBC_MAX_BINS>& actual,
                const std::vector<double>& expected) {
    for (std::size_t j = 0; j < expected.size(); ++j) {
        CAPTURE(name);
        CAPTURE(j);
        CHECK(actual[j] ==
              doctest::Approx(expected[j]).epsilon(SBC_RTOL).scale(0.0));
    }
}

void check_dsd(const sharp::SpectralBinDSD& dsd, const ExpectedDSD& expected,
               const float rime_factor) {
    REQUIRE(dsd.nbins() ==
            static_cast<std::ptrdiff_t>(expected.diameter.size()));
    CHECK(dsd.rime_factor() == rime_factor);
    for (std::size_t j = 0; j < expected.diameter.size(); ++j) {
        CHECK(dsd.diameter()[j] == expected.diameter[j]);
        CHECK(dsd.concentration()[j] == expected.concentration[j]);
    }

    const DiameterConstants& c = expected.constants;
    check_bins("mass", dsd.mass(), c.mass);
    check_bins("rain_fall_speed", dsd.rain_fall_speed(), c.rain_fall_speed);
    check_bins("pellet_fall_speed", dsd.pellet_fall_speed(),
               c.pellet_fall_speed);
    check_bins("foote_du_toit_fall_speed", dsd.foote_du_toit_fall_speed(),
               c.foote_du_toit_fall_speed);
    check_bins("liquid_capacitance_factor", dsd.liquid_capacitance_factor(),
               c.liquid_capacitance_factor);
    check_bins("liquid_length_factor", dsd.liquid_length_factor(),
               c.liquid_length_factor);

    const RimeConstants& r =
        (rime_factor == 5.0f) ? expected.rime_5 : expected.rime_1;
    check_bins("snow_diameter", dsd.snow_diameter(), r.snow_diameter);
    check_bins("snow_density", dsd.snow_density(), r.snow_density);
    check_bins("snow_aa", dsd.snow_aa(), r.snow_aa);
    check_bins("snow_bb", dsd.snow_bb(), r.snow_bb);
    check_bins("liquid_aa", dsd.liquid_aa(), r.liquid_aa);
    check_bins("liquid_bb", dsd.liquid_bb(), r.liquid_bb);
}

sharp::SpectralBinDSD build_dsd(const std::vector<float>& diameter,
                                const std::vector<float>& concentration,
                                const float rime_factor = 1.0f) {
    return sharp::spectral_bin_dsd(
        diameter.data(), concentration.data(),
        static_cast<std::ptrdiff_t>(diameter.size()), rime_factor);
}

void check_invalid(const sharp::SpectralBinDSD& dsd) {
    CHECK(dsd.nbins() == 0);
    CHECK(dsd.rime_factor() == sharp::MISSING);
    CHECK(dsd.diameter()[0] == 0.0f);
    CHECK(dsd.mass()[0] == 0.0f);
}
}  // namespace

TEST_CASE("Testing spectral_bin_dsd_default") {
    const sharp::SpectralBinDSD dsd = sharp::spectral_bin_dsd_default();
    REQUIRE(dsd.nbins() == 4);
    CHECK(dsd.diameter()[0] == 0.05f);
    CHECK(dsd.diameter()[1] == 0.75f);
    CHECK(dsd.diameter()[2] == 1.45f);
    CHECK(dsd.diameter()[3] == 2.15f);
    CHECK(dsd.concentration()[0] == 55.1843f);
    CHECK(dsd.concentration()[1] == 146.647f);
    CHECK(dsd.concentration()[2] == 11.6891f);
    CHECK(dsd.concentration()[3] == 3.60886f);
    check_dsd(dsd, PYTHON_DSD, 1.0f);
}

TEST_CASE("Testing spectral_bin_dsd per-bin constants") {
    for (const ExpectedDSD* expected :
         {&PYTHON_DSD, &CXX_2_0_3_DSD, &DELD_0_1_DSD, &EDGE_DSD}) {
        for (const float rime_factor : {1.0f, 5.0f}) {
            CAPTURE(expected->diameter.size());
            CAPTURE(rime_factor);
            check_dsd(build_dsd(expected->diameter, expected->concentration,
                                rime_factor),
                      *expected, rime_factor);
        }
    }
}

TEST_CASE("Testing spectral_bin_dsd at its limits") {
    const std::vector<float>& D = CXX_2_0_3_DSD.diameter;
    const std::vector<float>& N = CXX_2_0_3_DSD.concentration;
    CHECK(build_dsd(D, N, 1.0f).nbins() == 4);
    CHECK(build_dsd(D, N, 5.0f).nbins() == 4);
    CHECK(build_dsd({D[0]}, {N[0]}).nbins() == 1);
    CHECK(build_dsd(D, {0.0f, N[1], 0.0f, 0.0f}).nbins() == 4);
    CHECK(build_dsd({0.05f, 12.15f}, {1.0f, 1.0f}).nbins() == 2);

    std::vector<float> D64(sharp::SBC_MAX_BINS);
    for (std::size_t j = 0; j < D64.size(); ++j) {
        D64[j] = 0.05f + 0.1f * static_cast<float>(j);
    }
    const std::vector<float> N64(sharp::SBC_MAX_BINS, 1.0f);
    const sharp::SpectralBinDSD dsd = build_dsd(D64, N64);
    REQUIRE(dsd.nbins() == sharp::SBC_MAX_BINS);
    CHECK(dsd.diameter()[sharp::SBC_MAX_BINS - 1] == D64.back());
}

TEST_CASE("Testing spectral_bin_dsd rejects invalid distributions") {
    const std::vector<float>& D = CXX_2_0_3_DSD.diameter;
    const std::vector<float>& N = CXX_2_0_3_DSD.concentration;
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();

    check_invalid(sharp::SpectralBinDSD{});
    check_invalid(sharp::spectral_bin_dsd(D.data(), N.data(), 0));
    check_invalid(sharp::spectral_bin_dsd(D.data(), N.data(), -1));
    check_invalid(sharp::spectral_bin_dsd(nullptr, nullptr, 0));

    std::vector<float> D65(sharp::SBC_MAX_BINS + 1);
    for (std::size_t j = 0; j < D65.size(); ++j) {
        D65[j] = 0.05f + 0.1f * static_cast<float>(j);
    }
    check_invalid(build_dsd(D65, std::vector<float>(D65.size(), 1.0f)));

    for (const std::vector<float>& bad : {
             std::vector<float>{0.05f, 0.65f, 0.65f, 1.85f},
             std::vector<float>{0.05f, 0.65f, 0.6f, 1.85f},
             std::vector<float>{0.0f, 0.65f, 1.25f, 1.85f},
             std::vector<float>{-0.05f, 0.65f, 1.25f, 1.85f},
             std::vector<float>{0.05f, NaN, 1.25f, 1.85f},
             std::vector<float>{0.05f, 0.65f, 1.25f, inf},
             std::vector<float>{0.05f, 0.65f, 1.25f, 12.16f},
             std::vector<float>{0.05f, 0.65f, 1.25f, sharp::MISSING},
         }) {
        CAPTURE(bad[0]);
        CAPTURE(bad[1]);
        CAPTURE(bad[2]);
        CAPTURE(bad[3]);
        check_invalid(build_dsd(bad, N));
    }

    for (const std::vector<float>& bad : {
             std::vector<float>{55.1843f, -206.606f, 25.4924f, 3.60886f},
             std::vector<float>{0.0f, 0.0f, 0.0f, 0.0f},
             std::vector<float>{55.1843f, NaN, 25.4924f, 3.60886f},
             std::vector<float>{55.1843f, 206.606f, inf, 3.60886f},
             std::vector<float>{55.1843f, 206.606f, 25.4924f, sharp::MISSING},
         }) {
        CAPTURE(bad[0]);
        CAPTURE(bad[1]);
        CAPTURE(bad[2]);
        CAPTURE(bad[3]);
        check_invalid(build_dsd(D, bad));
    }

    for (const float rime_factor : {0.99f, 5.01f, 0.0f, NaN, sharp::MISSING}) {
        CAPTURE(rime_factor);
        check_invalid(build_dsd(D, N, rime_factor));
    }
}

// ---------------------------------------------------------------------------
// Cloud top from a sounding
// ---------------------------------------------------------------------------

namespace {
struct CloudTopSounding {
    std::vector<float> pressure;
    std::vector<float> height;
    std::vector<float> temperature;
    std::vector<float> dewpoint;
    std::vector<float> relh;
};

float cloud_top(const CloudTopSounding& s) {
    return sharp::spectral_bin_cloud_top(
        s.pressure.data(), s.height.data(), s.temperature.data(),
        s.dewpoint.data(), s.relh.data(),
        static_cast<std::ptrdiff_t>(s.height.size()));
}

// Levels every 1000 m at 270 K, from the dewpoint depression (K) and the
// relative humidity (fraction) of each level
CloudTopSounding depression_sounding(const std::vector<float>& depression,
                                     const std::vector<float>& relh) {
    CloudTopSounding s;
    for (std::size_t k = 0; k < depression.size(); ++k) {
        const float z = 1000.0f * static_cast<float>(k);
        s.pressure.push_back(100000.0f * std::exp(-z / 8000.0f));
        s.height.push_back(z);
        s.temperature.push_back(270.0f);
        s.dewpoint.push_back(270.0f - depression[k]);
    }
    s.relh = relh;
    return s;
}

CloudTopSounding without_level(CloudTopSounding s, const std::size_t k) {
    for (std::vector<float>* v :
         {&s.pressure, &s.height, &s.temperature, &s.dewpoint, &s.relh}) {
        v->erase(v->begin() + static_cast<std::ptrdiff_t>(k));
    }
    return s;
}

// Named cases of data/sbc_reference (levels.parquet)
const CloudTopSounding TOP_AT_HIGHEST_LEVEL{
    {100000.0f, 88200.0f, 77900.0f, 68700.0f, 60700.0f},
    {0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f},
    {277.0f, 268.0f, 269.0f, 280.0f, 265.0f},
    {277.0f, 268.0f, 269.0f, 280.0f, 265.0f},
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
};

const CloudTopSounding NO_CLOUD{
    {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f},
    {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f},
    {278.15f, 276.15f, 274.15f, 272.15f, 270.15f},
    {266.15f, 264.15f, 262.15f, 260.15f, 258.15f},
    {0.3f, 0.3f, 0.3f, 0.3f, 0.3f},
};

const CloudTopSounding DRY_LAYER_RESTART{
    {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
     68700.0f, 64600.0f, 60700.0f, 53500.0f},
    {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
     3500.0f, 4000.0f, 5000.0f},
    {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
     262.15f, 260.15f, 256.15f},
    {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 249.15f,
     247.15f, 260.15f, 256.15f},
    {0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.25f, 0.25f, 0.95f,
     0.95f},
};

const CloudTopSounding RELH_FALLBACK{
    {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f, 68700.0f},
    {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f},
    {271.15f, 269.15f, 266.65f, 265.15f, 262.65f, 260.15f, 257.15f},
    {263.15f, 261.15f, 258.65f, 257.15f, 254.65f, 252.15f, 249.15f},
    {0.5f, 0.5f, 0.5f, 0.85f, 0.85f, 0.5f, 0.5f},
};
}  // namespace

TEST_CASE("Testing spectral_bin_cloud_top on the golden named cases") {
    CHECK(cloud_top(TOP_AT_HIGHEST_LEVEL) == 4000.0f);
    CHECK(cloud_top(NO_CLOUD) == sharp::MISSING);
    CHECK(cloud_top(DRY_LAYER_RESTART) == 2500.0f);
    CHECK(cloud_top(RELH_FALLBACK) == 2000.0f);

    const CloudTopSounding& s = DRY_LAYER_RESTART;
    CHECK(sharp::spectral_bin_cloud_top(
              nullptr, s.height.data(), s.temperature.data(),
              s.dewpoint.data(), s.relh.data(),
              static_cast<std::ptrdiff_t>(s.height.size())) == 2500.0f);
}

TEST_CASE("Testing spectral_bin_cloud_top thresholds") {
    // The highest level is tested, over a level that is neither cloud, dry,
    // nor moist enough for the fallback.
    const auto top = [](const float depression, const float relh) {
        return cloud_top(depression_sounding({8.0f, depression}, {0.5f, relh}));
    };
    CHECK(top(6.0f, 0.9f) == 1000.0f);
    CHECK(top(6.5f, 0.9f) == sharp::MISSING);
    CHECK(top(2.0f, 0.60f) == sharp::MISSING);
    CHECK(top(2.0f, 0.61f) == 1000.0f);

    // A level between two cloud levels, dry or not
    const auto dry = [](const float depression, const float relh) {
        return cloud_top(depression_sounding({1.0f, depression, 0.0f},
                                             {0.9f, relh, 0.9f}));
    };
    CHECK(dry(10.0f, 0.5f) == 2000.0f);
    CHECK(dry(10.5f, 0.5f) == 0.0f);
    CHECK(dry(8.0f, 0.40f) == 2000.0f);
    CHECK(dry(8.0f, 0.39f) == 0.0f);

    // No cloud level, so the fallback decides.
    const auto fallback = [](const float relh_0, const float relh_1,
                             const float relh_2) {
        return cloud_top(depression_sounding({8.0f, 8.0f, 8.0f},
                                             {relh_0, relh_1, relh_2}));
    };
    CHECK(fallback(0.80f, 0.5f, 0.5f) == 0.0f);
    CHECK(fallback(0.79f, 0.5f, 0.5f) == sharp::MISSING);
    CHECK(fallback(0.85f, 0.85f, 0.5f) == 1000.0f);
    CHECK(fallback(0.5f, 0.5f, 0.9f) == sharp::MISSING);
    CHECK(fallback(0.85f, 0.5f, 0.9f) == 0.0f);
}

TEST_CASE("Testing spectral_bin_cloud_top driest level") {
    // Tied driest levels: the highest one wins.
    CHECK(cloud_top(depression_sounding({3.0f, 1.0f, 12.0f, 1.0f, 12.0f, 2.0f},
                                        {0.5f, 0.9f, 0.3f, 0.9f, 0.3f,
                                         0.9f})) == 3000.0f);

    // Negative depressions are not raised to 0, which would make the
    // first cloud top the driest level.
    CHECK(cloud_top(depression_sounding({-3.0f, -0.5f, -1.0f, -2.0f},
                                        {0.3f, 0.9f, 0.9f, 0.9f})) ==
          1000.0f);

    // Dry by relative humidity, with the driest level above it
    CHECK(cloud_top(depression_sounding({1.0f, 5.0f, 2.0f, 9.0f, 2.0f},
                                        {0.9f, 0.3f, 0.9f, 0.5f, 0.9f})) ==
          2000.0f);

    // No cloud level at or below the driest level keeps the first cloud top,
    // even with a cloud level between them.
    CHECK(cloud_top(depression_sounding({8.0f, 15.0f, 2.0f},
                                        {0.5f, 0.2f, 0.9f})) == 2000.0f);
    CHECK(cloud_top(depression_sounding({8.0f, 15.0f, 1.0f, 2.0f},
                                        {0.5f, 0.2f, 0.9f, 0.9f})) == 3000.0f);
}

TEST_CASE("Testing spectral_bin_cloud_top with few levels") {
    CHECK(sharp::spectral_bin_cloud_top(nullptr, nullptr, nullptr, nullptr,
                                        nullptr, 0) == sharp::MISSING);
    CHECK(sharp::spectral_bin_cloud_top(nullptr, nullptr, nullptr, nullptr,
                                        nullptr, -1) == sharp::MISSING);
    CHECK(cloud_top(depression_sounding({2.0f}, {0.9f})) == 0.0f);
    CHECK(cloud_top(depression_sounding({8.0f}, {0.9f})) == sharp::MISSING);
}

#ifndef NO_QC
TEST_CASE("Testing spectral_bin_cloud_top skips missing levels") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    // A missing value gives the cloud top of the sounding without its level.
    for (const CloudTopSounding& base : {
             TOP_AT_HIGHEST_LEVEL,
             NO_CLOUD,
             DRY_LAYER_RESTART,
             RELH_FALLBACK,
             depression_sounding({3.0f, 1.0f, 12.0f, 1.0f, 12.0f, 2.0f},
                                 {0.5f, 0.9f, 0.3f, 0.9f, 0.3f, 0.9f}),
             depression_sounding({-3.0f, -0.5f, -1.0f, -2.0f},
                                 {0.3f, 0.9f, 0.9f, 0.9f}),
             depression_sounding({8.0f, 8.0f, 8.0f}, {0.85f, 0.85f, 0.5f}),
         }) {
        for (std::vector<float> CloudTopSounding::*field :
             {&CloudTopSounding::temperature, &CloudTopSounding::dewpoint,
              &CloudTopSounding::relh}) {
            for (const float bad : {sharp::MISSING, NaN}) {
                for (std::size_t k = 0; k < base.height.size(); ++k) {
                    CAPTURE(base.height.size());
                    CAPTURE(k);
                    CAPTURE(bad);
                    CloudTopSounding s = base;
                    (s.*field)[k] = bad;
                    CHECK(cloud_top(s) == cloud_top(without_level(base, k)));
                }
            }
        }
    }

    CHECK(cloud_top(depression_sounding({8.0f, 8.0f, 8.0f},
                                        {0.85f, 0.85f, NaN})) == 0.0f);

    CloudTopSounding all_missing = DRY_LAYER_RESTART;
    for (float& t : all_missing.temperature) t = sharp::MISSING;
    CHECK(cloud_top(all_missing) == sharp::MISSING);
    all_missing = DRY_LAYER_RESTART;
    for (float& rh : all_missing.relh) rh = NaN;
    CHECK(cloud_top(all_missing) == sharp::MISSING);
}
#endif

// ---------------------------------------------------------------------------
// Precipitation type from a given cloud top: pre-classifier
// ---------------------------------------------------------------------------

namespace {
// From the surface up.
struct SBCProfile {
    std::vector<float> pressure;
    std::vector<float> height;
    std::vector<float> temperature;
    std::vector<float> dewpoint;
    std::vector<float> relh;
    std::vector<float> wetbulb;
};

SBCProfile saturated_profile(const std::vector<float>& height,
                             const std::vector<float>& wetbulb) {
    SBCProfile snd;
    for (const float z : height) {
        snd.pressure.push_back(100000.0f * std::exp(-z / 8000.0f));
    }
    snd.height = height;
    snd.temperature = wetbulb;
    snd.dewpoint = wetbulb;
    snd.relh.assign(wetbulb.size(), 1.0f);
    snd.wetbulb = wetbulb;
    return snd;
}

struct SBCRun {
    sharp::SpectralBinResult result;
    std::vector<float> profile;
};

// Asking for the profile must not change the result.
SBCRun run_sbc(
    const SBCProfile& snd, const float cloud_top,
    const sharp::SpectralBinDSD& dsd = sharp::spectral_bin_dsd_default(),
    const float tice = sharp::SBC_ICE_NUCLEATION_TEMPERATURE) {
    const auto N = static_cast<std::ptrdiff_t>(snd.height.size());
    SBCRun run;
    run.profile.assign(
        snd.height.size() * static_cast<std::size_t>(dsd.nbins()), 0.0f);
    run.result = sharp::spectral_bin_classifier(
        snd.pressure.data(), snd.height.data(), snd.temperature.data(),
        snd.dewpoint.data(), snd.relh.data(), snd.wetbulb.data(), N,
        cloud_top, dsd, tice, run.profile.data());
    const sharp::SpectralBinResult without = sharp::spectral_bin_classifier(
        snd.pressure.data(), snd.height.data(), snd.temperature.data(),
        snd.dewpoint.data(), snd.relh.data(), snd.wetbulb.data(), N,
        cloud_top, dsd, tice);
    CHECK(without.precip_type == run.result.precip_type);
    CHECK(without.liquid_fraction == run.result.liquid_fraction);
    CHECK(without.supercooled_liquid_height ==
          run.result.supercooled_liquid_height);
    return run;
}

constexpr sharp::SpectralBinResult SBC_SN{sharp::PrecipType::snow, 0.0f,
                                          sharp::MISSING};
constexpr sharp::SpectralBinResult SBC_FZRA{sharp::PrecipType::freezing_rain,
                                            1.0f, 0.0f};
constexpr sharp::SpectralBinResult SBC_RA{sharp::PrecipType::rain, 1.0f,
                                          sharp::MISSING};
constexpr sharp::SpectralBinResult SBC_MISSING{};

void check_sbc(const SBCRun& run, const sharp::SpectralBinResult& expected) {
    CHECK(run.result.precip_type == expected.precip_type);
    CHECK(run.result.liquid_fraction == expected.liquid_fraction);
    CHECK(run.result.supercooled_liquid_height ==
          expected.supercooled_liquid_height);
    std::size_t not_missing = 0;
    for (const float lf : run.profile) not_missing += (lf != sharp::MISSING);
    CHECK(not_missing == 0);
}
}  // namespace

// The named cases of data/sbc_reference, one per pre-classifier path, with
// their reference cloud tops.

TEST_CASE("Testing spectral_bin_classifier pre-classifies snow") {
    // Case 0
    const SBCProfile snd{
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f, 68700.0f,
         64600.0f, 60700.0f, 57000.0f, 53500.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f, 3500.0f,
         4000.0f, 4500.0f, 5000.0f},
        {270.15f, 269.15f, 268.15f, 266.15f, 265.15f, 264.15f, 262.15f,
         260.15f, 258.15f, 256.15f, 253.15f},
        {270.15f, 269.15f, 268.15f, 266.15f, 265.15f, 264.15f, 262.15f,
         260.15f, 258.15f, 256.15f, 253.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {270.15f, 269.15f, 268.15f, 266.15f, 265.15f, 264.15f, 262.15f,
         260.15f, 258.15f, 256.15f, 253.15f},
    };
    check_sbc(run_sbc(snd, 5000.0f), SBC_SN);
}

TEST_CASE("Testing spectral_bin_classifier subfreezing column, warm top") {
    // Case 1
    const SBCProfile snd{
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f},
        {272.15f, 271.15f, 270.15f, 269.15f, 268.65f},
        {272.15f, 271.15f, 270.15f, 269.15f, 268.65f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {272.15f, 271.15f, 270.15f, 269.15f, 268.65f},
    };
    check_sbc(run_sbc(snd, 2000.0f), SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier subfreezing column, warm 3 km") {
    // Case 2
    const SBCProfile snd{
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f, 68700.0f,
         64600.0f, 60700.0f, 57000.0f, 53500.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f, 3500.0f,
         4000.0f, 4500.0f, 5000.0f},
        {272.15f, 271.65f, 271.15f, 270.65f, 270.15f, 269.65f, 269.15f,
         267.65f, 265.15f, 262.15f, 259.15f},
        {272.15f, 271.65f, 271.15f, 270.65f, 270.15f, 269.65f, 269.15f,
         267.65f, 265.15f, 262.15f, 259.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {272.15f, 271.65f, 271.15f, 270.65f, 270.15f, 269.65f, 269.15f,
         267.65f, 265.15f, 262.15f, 259.15f},
    };
    check_sbc(run_sbc(snd, 5000.0f), SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier warm top over a cold surface") {
    // Case 3
    const SBCProfile snd{
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f},
        {271.15f, 272.15f, 275.15f, 276.15f, 274.15f, 270.15f},
        {271.15f, 272.15f, 275.15f, 276.15f, 274.15f, 270.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {271.15f, 272.15f, 275.15f, 276.15f, 274.15f, 270.15f},
    };
    check_sbc(run_sbc(snd, 2500.0f), SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier pre-classifies rain") {
    // Case 4
    const SBCProfile snd{
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f},
        {281.15f, 279.15f, 277.15f, 276.15f, 275.15f},
        {281.15f, 279.15f, 277.15f, 276.15f, 275.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {281.15f, 279.15f, 277.15f, 276.15f, 275.15f},
    };
    check_sbc(run_sbc(snd, 2000.0f), SBC_RA);
}

TEST_CASE("Testing spectral_bin_classifier with a maximum Tw of 0 C") {
    const SBCProfile snd =
        saturated_profile({0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f},
                          {272.15f, sharp::ZEROCNK, 270.15f, 265.15f, 262.15f});
    check_sbc(run_sbc(snd, 4000.0f), SBC_SN);
}

TEST_CASE("Testing spectral_bin_classifier rule 2 at a surface below 0 C") {
    const SBCProfile snd =
        saturated_profile({0.0f, 1000.0f, 2000.0f, 3000.0f},
                          {273.1f, 275.15f, 276.15f, 274.15f});
    check_sbc(run_sbc(snd, 3000.0f), SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier with Tw at Tice") {
    const float tice = sharp::SBC_ICE_NUCLEATION_TEMPERATURE;
    const std::vector<float> height = {0.0f, 1000.0f, 2000.0f, 3000.0f,
                                       4000.0f};
    check_sbc(run_sbc(saturated_profile(
                          height, {270.15f, 268.15f, 262.15f, 260.15f, tice}),
                      4000.0f),
              SBC_FZRA);
    check_sbc(run_sbc(saturated_profile(height, {270.15f, 268.15f, 262.15f,
                                                 260.15f, 267.1f}),
                      4000.0f),
              SBC_SN);

    check_sbc(run_sbc(saturated_profile(
                          height, {270.15f, 268.15f, tice, 269.15f, 262.15f}),
                      4000.0f),
              SBC_FZRA);
    check_sbc(run_sbc(saturated_profile(height, {270.15f, 268.15f, 267.1f,
                                                 269.15f, 262.15f}),
                      4000.0f),
              SBC_SN);

    check_sbc(run_sbc(saturated_profile({0.0f, 1000.0f, 2000.0f, 3000.0f},
                                        {270.15f, 275.15f, 276.15f, 267.2f}),
                      3000.0f),
              SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier ties in the 3 km window") {
    const std::vector<float> tw = {272.15f, 271.15f, 270.15f,
                                   268.15f, 265.15f, 262.15f};
    for (const float surface : {0.0f, 1000.0f}) {
        CAPTURE(surface);
        const auto run = [&](const float below, const float above) {
            return run_sbc(
                saturated_profile({surface, surface + 1000.0f,
                                   surface + 2000.0f, surface + below,
                                   surface + above, surface + 4000.0f},
                                  tw),
                surface + 4000.0f);
        };
        check_sbc(run(2900.0f, 3100.0f), SBC_SN);
        check_sbc(run(2901.0f, 3100.0f), SBC_FZRA);
        check_sbc(run(2900.0f, 3099.0f), SBC_SN);
    }

    check_sbc(run_sbc(saturated_profile({0.0f, 6000.0f}, {270.15f, 262.15f}),
                      6000.0f),
              SBC_SN);
    check_sbc(run_sbc(saturated_profile({0.0f, 6001.0f}, {270.15f, 262.15f}),
                      6001.0f),
              SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier with a cloud top below 3 km") {
    const SBCProfile snd = saturated_profile(
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f, 3500.0f},
        {270.15f, 269.15f, 268.15f, 268.15f, 266.15f, 275.15f, 262.15f,
         260.15f});
    check_sbc(run_sbc(snd, 2000.0f), SBC_SN);
}

TEST_CASE("Testing spectral_bin_classifier with Tice -10 C") {
    const sharp::SpectralBinDSD dsd = sharp::spectral_bin_dsd_default();
    const std::vector<float> height = {0.0f, 1000.0f, 2000.0f, 3000.0f,
                                       4000.0f};
    const SBCProfile window = saturated_profile(
        height, {270.15f, 268.15f, 266.15f, 264.15f, 262.15f});
    check_sbc(run_sbc(window, 4000.0f, dsd, 267.15f), SBC_SN);
    check_sbc(run_sbc(window, 4000.0f, dsd, 263.15f), SBC_FZRA);

    const SBCProfile top = saturated_profile(
        height, {270.15f, 268.15f, 266.15f, 262.15f, 265.15f});
    check_sbc(run_sbc(top, 4000.0f, dsd, 267.15f), SBC_SN);
    check_sbc(run_sbc(top, 4000.0f, dsd, 263.15f), SBC_FZRA);
}

TEST_CASE("Testing spectral_bin_classifier cloud-top handling") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    const std::vector<float> tw = {272.15f, 270.15f, 266.15f, 268.15f,
                                   262.15f};
    for (const float surface : {0.0f, 1500.0f}) {
        CAPTURE(surface);
        const SBCProfile snd = saturated_profile(
            {surface, surface + 1000.0f, surface + 2000.0f, surface + 3000.0f,
             surface + 4000.0f},
            tw);
        check_sbc(run_sbc(snd, surface + 1000.0f), SBC_FZRA);
        check_sbc(run_sbc(snd, surface + 1500.0f), SBC_FZRA);
        check_sbc(run_sbc(snd, surface + 2000.0f), SBC_SN);
        check_sbc(run_sbc(snd, surface + 2999.0f), SBC_SN);
        check_sbc(run_sbc(snd, surface + 3000.0f), SBC_FZRA);
        check_sbc(run_sbc(snd, surface + 3500.0f), SBC_FZRA);
        check_sbc(run_sbc(snd, surface + 4000.0f), SBC_SN);
        check_sbc(run_sbc(snd, surface + 20000.0f), SBC_SN);

        check_sbc(run_sbc(snd, surface + 999.0f), SBC_MISSING);
        check_sbc(run_sbc(snd, surface), SBC_MISSING);
        check_sbc(run_sbc(snd, surface - 1.0f), SBC_MISSING);

        check_sbc(run_sbc(snd, sharp::MISSING), SBC_MISSING);
        check_sbc(run_sbc(snd, NaN), SBC_MISSING);
    }
}

TEST_CASE("Testing spectral_bin_classifier ignores data above the cloud top") {
    const SBCProfile snd =
        saturated_profile({0.0f, 1000.0f, 2000.0f, 3000.0f},
                          {280.15f, 278.15f, 276.15f, sharp::MISSING});
    check_sbc(run_sbc(snd, 2000.0f), SBC_RA);
    check_sbc(run_sbc(snd, 2999.0f), SBC_RA);
}

TEST_CASE("Testing spectral_bin_classifier with invalid inputs") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    const SBCProfile snd = saturated_profile(
        {0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f},
        {272.15f, 270.15f, 266.15f, 268.15f, 262.15f});
    const sharp::SpectralBinDSD dsd = sharp::spectral_bin_dsd_default();
    check_sbc(run_sbc(snd, 4000.0f, dsd), SBC_SN);

    const SBCRun invalid_dsd = run_sbc(snd, 4000.0f, sharp::SpectralBinDSD{});
    check_sbc(invalid_dsd, SBC_MISSING);
    CHECK(invalid_dsd.profile.empty());

    for (const float tice : {sharp::MISSING, NaN, 0.0f, -1.0f}) {
        CAPTURE(tice);
        check_sbc(run_sbc(snd, 4000.0f, dsd, tice), SBC_MISSING);
    }

    // N < 2 returns before reading an input array.
    for (const std::ptrdiff_t N : {-1, 0, 1}) {
        CAPTURE(N);
        std::array<float, 4> profile = {0.0f, 0.0f, 0.0f, 0.0f};
        const sharp::SpectralBinResult result = sharp::spectral_bin_classifier(
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, N, 4000.0f,
            dsd, sharp::SBC_ICE_NUCLEATION_TEMPERATURE, profile.data());
        CHECK(result.precip_type == sharp::PrecipType::missing);
        CHECK(result.liquid_fraction == sharp::MISSING);
        CHECK(result.supercooled_liquid_height == sharp::MISSING);
        for (const float lf : profile) {
            CHECK(lf == ((N == 1) ? sharp::MISSING : 0.0f));
        }
    }
}

#ifndef NO_QC
TEST_CASE("Testing spectral_bin_classifier skips missing levels") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    using Field = std::vector<float> SBCProfile::*;
    const Field fields[] = {&SBCProfile::temperature, &SBCProfile::dewpoint,
                            &SBCProfile::relh, &SBCProfile::wetbulb};

    const SBCProfile window = saturated_profile(
        {0.0f, 1000.0f, 2500.0f, 3000.0f, 4000.0f, 5000.0f},
        {272.15f, 270.15f, 268.15f, 264.15f, 262.15f, 260.15f});
    check_sbc(run_sbc(window, 5000.0f), SBC_SN);

    const SBCProfile surface = saturated_profile(
        {0.0f, 500.0f, 1000.0f, 3200.0f, 3600.0f, 5000.0f},
        {272.15f, 271.15f, 270.15f, 268.15f, 264.15f, 260.15f});
    check_sbc(run_sbc(surface, 5000.0f), SBC_FZRA);

    const SBCProfile top = saturated_profile(
        {0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f},
        {272.15f, 270.15f, 266.15f, 268.15f, 262.15f});
    check_sbc(run_sbc(top, 4000.0f), SBC_SN);
    check_sbc(run_sbc(top, 20000.0f), SBC_SN);

    for (const Field field : fields) {
        for (const float bad : {sharp::MISSING, NaN}) {
            CAPTURE(bad);
            SBCProfile snd = window;
            (snd.*field)[3] = bad;
            check_sbc(run_sbc(snd, 5000.0f), SBC_FZRA);

            snd = surface;
            (snd.*field)[0] = bad;
            check_sbc(run_sbc(snd, 5000.0f), SBC_SN);

            snd = top;
            (snd.*field)[4] = bad;
            check_sbc(run_sbc(snd, 4000.0f), SBC_FZRA);
            check_sbc(run_sbc(snd, 20000.0f), SBC_FZRA);
        }
    }
}

TEST_CASE("Testing spectral_bin_classifier with fewer than 2 valid levels") {
    const SBCProfile snd = saturated_profile(
        {0.0f, 1000.0f, 2000.0f}, {272.15f, 270.15f, 266.15f});
    check_sbc(run_sbc(snd, 2000.0f), SBC_SN);
    check_sbc(run_sbc(snd, 1000.0f), SBC_FZRA);

    SBCProfile one = snd;
    one.wetbulb[1] = sharp::MISSING;
    one.temperature[2] = sharp::MISSING;
    check_sbc(run_sbc(one, 2000.0f), SBC_MISSING);
    check_sbc(run_sbc(one, 1000.0f), SBC_MISSING);

    SBCProfile no_surface = snd;
    no_surface.relh[0] = sharp::MISSING;
    check_sbc(run_sbc(no_surface, 2000.0f), SBC_SN);
    check_sbc(run_sbc(no_surface, 1000.0f), SBC_MISSING);

    SBCProfile none = snd;
    none.dewpoint.assign(3, sharp::MISSING);
    check_sbc(run_sbc(none, 2000.0f), SBC_MISSING);
}
#endif

// ---------------------------------------------------------------------------
// Microphysics: frozen cloud tops and melting
// ---------------------------------------------------------------------------

namespace {
// The reference's results for a column, from data/sbc_reference or from the
// reference run on a hand-made column. The profile is levels x bins.
struct SBCGolden {
    SBCProfile snd;
    float cloud_top;
    const ExpectedDSD* dsd;
    float rime_factor;
    float tice;
    int crossings;
    sharp::PrecipType precip_type;
    double liquid_fraction;
    float supercooled_liquid_height;
    std::vector<float> profile;
};

SBCProfile with_relh(SBCProfile snd, const std::vector<float>& relh) {
    snd.relh = relh;
    return snd;
}

// The reference's decision tree on a liquid fraction
sharp::PrecipType sbc_decision(const float liquid, const int crossings,
                               const bool warm) {
    using sharp::PrecipType;
    const bool all_ice = (liquid == 0.0f);
    const bool all_liquid = (liquid == 1.0f);
    if (warm && (crossings == 1)) {
        if (all_liquid || (liquid > 0.85f)) return PrecipType::rain;
        if (all_ice || (liquid < 0.60f)) return PrecipType::snow;
        return PrecipType::rain_snow;
    }
    if (all_ice || (liquid < 0.15f)) return PrecipType::ice_pellets;
    if (warm) {
        return (all_liquid || (1.0f - liquid < 0.15f))
                   ? PrecipType::rain
                   : PrecipType::rain_ice_pellets;
    }
    return (all_liquid || (liquid > 0.85f))
               ? PrecipType::freezing_rain
               : PrecipType::freezing_rain_ice_pellets;
}

// The golden-data comparison rules for a case without flags, whose surface
// is level 0
void check_golden(const SBCGolden& golden) {
    const sharp::SpectralBinDSD dsd = build_dsd(
        golden.dsd->diameter, golden.dsd->concentration, golden.rime_factor);
    const SBCRun run = run_sbc(golden.snd, golden.cloud_top, dsd, golden.tice);
    const sharp::SpectralBinResult& result = run.result;
    CHECK(result.precip_type ==
          sbc_decision(result.liquid_fraction, golden.crossings,
                       golden.snd.wetbulb[0] > sharp::ZEROCNK));
    CHECK(result.precip_type == golden.precip_type);
    CHECK(std::abs(result.liquid_fraction - golden.liquid_fraction) <= 1e-4);
    CHECK(result.supercooled_liquid_height ==
          golden.supercooled_liquid_height);
    REQUIRE(run.profile.size() == golden.profile.size());
    for (std::size_t i = 0; i < golden.profile.size(); ++i) {
        CAPTURE(i);
        if (golden.profile[i] == sharp::MISSING) {
            CHECK(run.profile[i] == sharp::MISSING);
        } else {
            CHECK(std::abs(run.profile[i] - golden.profile[i]) <= 1e-3f);
        }
    }
}

// The reference's results for the column with a level inserted at index k,
// which the classifier skips
SBCGolden with_skipped_level(SBCGolden golden, const std::size_t k,
                             const float height, const float wetbulb) {
    SBCProfile& snd = golden.snd;
    const auto at = static_cast<std::ptrdiff_t>(k);
    snd.pressure.insert(snd.pressure.begin() + at,
                        100000.0f * std::exp(-height / 8000.0f));
    snd.height.insert(snd.height.begin() + at, height);
    snd.temperature.insert(snd.temperature.begin() + at, wetbulb);
    snd.dewpoint.insert(snd.dewpoint.begin() + at, wetbulb);
    snd.relh.insert(snd.relh.begin() + at, 1.0f);
    snd.wetbulb.insert(snd.wetbulb.begin() + at, wetbulb);
    const std::size_t nbins = golden.dsd->diameter.size();
    golden.profile.insert(
        golden.profile.begin() + static_cast<std::ptrdiff_t>(k * nbins), nbins,
        sharp::MISSING);
    return golden;
}

// data/sbc_reference case 5
const SBCGolden SBC_CORE_RA{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.95f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.95f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {273.95f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::rain,
    0.9121556184471886,
    sharp::MISSING,
    {1.0f, 1.0f, 1.0f, 0.7434978f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// data/sbc_reference case 6
const SBCGolden SBC_CORE_SN{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.35f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.35f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {273.35f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::snow,
    0.36193448304389836,
    sharp::MISSING,
    {1.0f, 0.61596805f, 0.2794438f, 0.18456055f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// data/sbc_reference case 7
const SBCGolden SBC_CORE_RASN{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::rain_snow,
    0.6903628867814539,
    sharp::MISSING,
    {1.0f, 1.0f, 0.63058996f, 0.41648102f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// data/sbc_reference case 16
const SBCGolden SBC_CORE_DRY_LAYER{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 249.15f,
         247.15f, 260.15f, 256.15f},
        {0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.95f, 0.25f, 0.25f, 0.95f,
         0.95f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f,
         258.71643f, 256.98294f, 260.15f, 256.15f},
    },
    2500.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::snow,
    0.24678250549865566,
    sharp::MISSING,
    {1.0f, 0.42426872f, 0.1924805f, 0.12712616f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING},
};

// data/sbc_reference case 18
const SBCGolden SBC_CORE_CXX_DSD{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
    },
    5000.0f,
    &CXX_2_0_3_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::rain_snow,
    0.7686941434936196,
    sharp::MISSING,
    {1.0f, 1.0f, 0.7457145f, 0.48528457f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// data/sbc_reference case 20
const SBCGolden SBC_CORE_RIME_5{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    5.0f,
    267.15f,
    1,
    sharp::PrecipType::snow,
    0.23181082327704985,
    sharp::MISSING,
    {1.0f, 0.4763523f, 0.13665824f, 0.08876687f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// data/sbc_reference case 23
const SBCGolden SBC_CORE_TICE_ALT{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {273.6f, 272.15f, 271.15f, 270.15f, 269.15f, 268.15f, 266.15f, 264.15f,
         262.15f, 260.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    1.0f,
    263.15f,
    1,
    sharp::PrecipType::rain_snow,
    0.6903628867814539,
    sharp::MISSING,
    {1.0f, 1.0f, 0.63058996f, 0.41648102f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// Hand-made columns, with the results of the reference on the same values
// in float64. 273.15f stands for 273.15 K exactly, which the golden data
// leave out.

// Melting just below a frozen cloud top: RA
const SBCGolden SBC_MELTING_BELOW_TOP{
    saturated_profile({0.0f, 500.0f, 1000.0f}, {273.5f, 273.25f, 266.0f}),
    1000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::rain,
    0.920670924307368,
    sharp::MISSING,
    {1.0f, 1.0f, 1.0f, 0.7377472f, 1.0f, 0.6033493f, 0.2741206f, 0.18115126f,
     0.0f, 0.0f, 0.0f, 0.0f},
};

// Two crossings and a surface at 0 C: PL
const SBCGolden SBC_ZERO_C_PL{
    saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.15f, 273.1875f, 271.5f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::ice_pellets,
    0.0657364244651753,
    sharp::MISSING,
    {1.0f, 0.11415317f, 0.051824972f, 0.03423811f, 1.0f, 0.11415317f,
     0.051824972f, 0.03423811f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// Two crossings and a surface at 0 C: FZRAPL
const SBCGolden SBC_ZERO_C_FZRAPL{
    saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.15f, 273.25f, 271.5f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::freezing_rain_ice_pellets,
    0.17979090352387,
    0.0f,
    {1.0f, 0.30462942f, 0.13830099f, 0.091368586f, 1.0f, 0.30462942f,
     0.13830099f, 0.091368586f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// Two crossings and a surface at 0 C: FZRA
const SBCGolden SBC_ZERO_C_FZRA{
    saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.15f, 274.0f, 271.5f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::freezing_rain,
    0.8854490179155798,
    0.0f,
    {1.0f, 1.0f, 1.0f, 0.7835615f, 1.0f, 1.0f, 1.0f, 0.7835615f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// Three crossings and a warm surface: RAPL
const SBCGolden SBC_ZERO_C_RAPL{
    saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.5f, 273.15f, 273.25f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    3,
    sharp::PrecipType::rain_ice_pellets,
    0.7436840013308715,
    sharp::MISSING,
    {1.0f, 1.0f, 0.6673701f, 0.45883283f, 1.0f, 0.6033493f, 0.2741206f,
     0.18115126f, 1.0f, 0.6033493f, 0.2741206f, 0.18115126f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f},
};

// Three crossings and a warm surface: RA
const SBCGolden SBC_ZERO_C_RA{
    saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.5f, 273.15f, 273.5f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    3,
    sharp::PrecipType::rain,
    0.9279022075581568,
    sharp::MISSING,
    {1.0f, 1.0f, 1.0f, 0.7704413f, 1.0f, 1.0f, 0.96224207f, 0.6358984f, 1.0f,
     1.0f, 0.96224207f, 0.6358984f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// Three crossings and a warm surface, with no melting in dry air: PL
const SBCGolden SBC_ZERO_C_DRY_PL{
    with_relh(saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.5f, 273.15f, 273.5f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
        {0.5f, 0.5f, 0.5f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    3,
    sharp::PrecipType::ice_pellets,
    0.0,
    sharp::MISSING,
    {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f},
};

// 0 C just above the first crossing melts: RASN
const SBCGolden SBC_ZERO_C_ABOVE_CROSSING{
    saturated_profile(
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f},
        {273.5f, 273.25f, 273.15f, 270.0f, 269.0f, 268.0f, 266.0f, 264.0f,
         263.0f, 262.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::rain_snow,
    0.6830516109364516,
    sharp::MISSING,
    {1.0f, 1.0f, 0.58917785f, 0.39497298f, 1.0f, 0.32282987f, 0.14180872f,
     0.0926089f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};
}  // namespace

TEST_CASE("Testing spectral_bin_classifier core: RA, SN, and RASN") {
    SUBCASE("RA") { check_golden(SBC_CORE_RA); }
    SUBCASE("SN") { check_golden(SBC_CORE_SN); }
    SUBCASE("RASN") { check_golden(SBC_CORE_RASN); }
}

TEST_CASE("Testing spectral_bin_classifier core: DSD, rime factor, and Tice") {
    SUBCASE("C++ 2.0.3 DSD") { check_golden(SBC_CORE_CXX_DSD); }
    SUBCASE("rime factor 5") { check_golden(SBC_CORE_RIME_5); }
    SUBCASE("Tice -10 C") { check_golden(SBC_CORE_TICE_ALT); }
}

TEST_CASE("Testing spectral_bin_classifier core below a dry layer") {
    check_golden(SBC_CORE_DRY_LAYER);
}

TEST_CASE("Testing spectral_bin_classifier core melting below the top") {
    check_golden(SBC_MELTING_BELOW_TOP);
}

TEST_CASE("Testing spectral_bin_classifier core surface decision") {
    SUBCASE("PL") { check_golden(SBC_ZERO_C_PL); }
    SUBCASE("FZRAPL") { check_golden(SBC_ZERO_C_FZRAPL); }
    SUBCASE("FZRA") { check_golden(SBC_ZERO_C_FZRA); }
    SUBCASE("RAPL") { check_golden(SBC_ZERO_C_RAPL); }
    SUBCASE("RA") { check_golden(SBC_ZERO_C_RA); }
    SUBCASE("PL, warm surface") { check_golden(SBC_ZERO_C_DRY_PL); }
}

TEST_CASE("Testing spectral_bin_classifier core with 0 C above a crossing") {
    check_golden(SBC_ZERO_C_ABOVE_CROSSING);
}

TEST_CASE("Testing spectral_bin_classifier core ignores data above the top") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    SBCGolden golden = SBC_CORE_RA;
    for (const float z : {6000.0f, 7000.0f}) {
        golden.snd.pressure.push_back(100000.0f * std::exp(-z / 8000.0f));
        golden.snd.height.push_back(z);
        golden.snd.temperature.push_back(NaN);
        golden.snd.dewpoint.push_back(sharp::MISSING);
        golden.snd.relh.push_back(NaN);
        golden.snd.wetbulb.push_back(300.0f);
        golden.profile.insert(golden.profile.end(), 4, sharp::MISSING);
    }
    check_golden(golden);
}

#ifndef NO_QC
TEST_CASE("Testing spectral_bin_classifier core skips missing levels") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    using Field = std::vector<float> SBCProfile::*;
    const Field fields[] = {&SBCProfile::temperature, &SBCProfile::dewpoint,
                            &SBCProfile::relh, &SBCProfile::wetbulb};
    for (const Field field : fields) {
        for (const float bad : {sharp::MISSING, NaN}) {
            CAPTURE(bad);
            SBCGolden golden = with_skipped_level(SBC_ZERO_C_ABOVE_CROSSING, 1,
                                                  125.0f, 300.0f);
            golden = with_skipped_level(golden, 4, 750.0f, 300.0f);
            (golden.snd.*field)[1] = bad;
            (golden.snd.*field)[4] = bad;
            check_golden(golden);
        }
    }
}
#endif

// ---------------------------------------------------------------------------
// Microphysics: refreezing
// ---------------------------------------------------------------------------

namespace {
// data/sbc_reference case 1000
const SBCGolden SBC_SAMPLE{
    {
        {98210.0f, 97500.0f, 95000.0f, 92500.0f, 90000.0f, 87500.0f, 85000.0f,
         82500.0f, 80000.0f, 77500.0f, 75000.0f, 72500.0f, 70000.0f, 67500.0f,
         65000.0f, 62500.0f, 60000.0f, 57500.0f, 55000.0f, 52500.0f, 50000.0f,
         47500.0f, 45000.0f, 42500.0f, 40000.0f, 37500.0f, 35000.0f, 32500.0f,
         30000.0f, 27500.0f, 25000.0f, 22500.0f, 20000.0f, 17500.0f, 15000.0f,
         12500.0f, 10000.0f, 7500.0f, 5000.0f},
        {0.0f, 56.442383f, 260.70694f, 469.10275f, 681.9126f, 901.26404f,
         1130.5469f, 1369.6487f, 1616.9025f, 1871.5979f, 2133.8792f, 2403.5466f,
         2680.673f, 2966.7834f, 3261.1917f, 3565.882f, 3881.2122f, 4207.1167f,
         4545.4614f, 4895.5923f, 5261.2954f, 5641.041f, 6039.1777f, 6457.2085f,
         6895.511f, 7356.904f, 7843.8374f, 8359.449f, 8908.406f, 9495.977f,
         10127.756f, 10809.0205f, 11555.859f, 12411.309f, 13401.112f,
         14553.105f, 15933.92f, 17685.596f, 20181.354f},
        {269.5766f, 269.07043f, 267.5222f, 265.7832f, 265.22726f, 267.7831f,
         271.326f, 273.35565f, 273.753f, 272.99402f, 271.7636f, 270.32312f,
         268.6287f, 267.1676f, 265.66608f, 264.2079f, 262.53174f, 260.54184f,
         258.4745f, 256.37653f, 254.39236f, 252.43417f, 250.50314f, 248.35762f,
         245.80988f, 242.88681f, 239.50732f, 236.00835f, 232.58592f, 228.92673f,
         224.08499f, 218.03888f, 216.89758f, 220.19148f, 217.91052f, 214.21368f,
         208.7692f, 209.2101f, 211.9147f},
        {267.84027f, 267.2251f, 267.17688f, 265.34628f, 264.62677f, 267.3271f,
         270.72638f, 272.98447f, 273.4998f, 272.6873f, 271.4373f, 270.1248f,
         268.5623f, 267.1248f, 265.3748f, 263.4373f, 261.3748f, 259.1873f,
         256.9998f, 254.56229f, 252.43729f, 250.18729f, 248.12479f, 245.81229f,
         243.06229f, 239.81229f, 236.18729f, 232.49979f, 228.81229f, 224.62479f,
         217.56229f, 203.68729f, 199.31229f, 192.49979f, 192.12479f, 192.12479f,
         192.12479f, 192.12479f, 192.12479f},
        {0.874f, 0.8661301f, 0.9670676f, 0.9688999f, 0.9583809f, 0.9959251f,
         0.99636406f, 0.9963686f, 0.9968148f, 0.9961528f, 0.9861289f,
         0.9834164f, 0.9911664f, 0.98461264f, 0.9653304f, 0.93880683f,
         0.9168299f, 0.8952456f, 0.8777837f, 0.8567047f, 0.8419585f, 0.8228494f,
         0.8083167f, 0.7945054f, 0.775455f, 0.75370675f, 0.7282403f, 0.7048327f,
         0.6836666f, 0.6608954f, 0.6338591f, 0.6003386f, 0.4446618f,
         0.036687344f, 0.023816185f, 0.022249507f, 0.029207146f, 0.021175709f,
         0.010973916f},
        {268.94424f, 268.40955f, 267.4046f, 265.6432f, 265.03622f, 267.6174f,
         271.07037f, 273.18256f, 273.63135f, 272.84793f, 271.612f, 270.23395f,
         268.6001f, 267.1498f, 265.54956f, 263.91086f, 262.10638f, 260.0743f,
         257.99957f, 255.83368f, 253.84691f, 251.84987f, 249.92656f, 247.79037f,
         245.26207f, 242.35379f, 239.02518f, 235.58836f, 232.21382f, 228.5855f,
         223.70995f, 217.5039f, 216.24294f, 218.64377f, 216.56673f, 213.25244f,
         208.24313f, 208.49329f, 210.38112f},
    },
    9495.977f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::ice_pellets,
    0.0,
    1130.5469f,
    {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
     1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.62460726f, 0.41479364f, 1.0f, 1.0f,
     0.6179453f, 0.40897992f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING, sharp::MISSING, sharp::MISSING, sharp::MISSING,
     sharp::MISSING},
};

// data/sbc_reference case 10
const SBCGolden SBC_TICE_SWITCH{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {270.15f, 265.15f, 265.15f, 271.15f, 275.15f, 274.15f, 272.15f, 270.15f,
         265.15f, 261.15f, 256.15f},
        {270.15f, 265.15f, 265.15f, 271.15f, 275.15f, 274.15f, 272.15f, 270.15f,
         265.15f, 261.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {270.15f, 265.15f, 265.15f, 271.15f, 275.15f, 274.15f, 272.15f, 270.15f,
         265.15f, 261.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::freezing_rain,
    1.0,
    0.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// data/sbc_reference case 11
const SBCGolden SBC_REMELT{
    {
        {100000.0f, 96900.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f,
         68700.0f, 64600.0f, 60700.0f, 53500.0f},
        {0.0f, 250.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f,
         3500.0f, 4000.0f, 5000.0f},
        {274.15f, 272.15f, 269.15f, 265.15f, 271.15f, 273.55f, 272.15f, 270.15f,
         265.15f, 261.15f, 256.15f},
        {274.15f, 272.15f, 269.15f, 265.15f, 271.15f, 273.55f, 272.15f, 270.15f,
         265.15f, 261.15f, 256.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {274.15f, 272.15f, 269.15f, 265.15f, 271.15f, 273.55f, 272.15f, 270.15f,
         265.15f, 261.15f, 256.15f},
    },
    5000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    3,
    sharp::PrecipType::rain_ice_pellets,
    0.46042058216143616,
    1500.0f,
    {1.0f, 1.0f, 0.25688764f, 0.280271f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
     1.0f, 0.6902354f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
     0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// Hand-made columns, with the results of the reference on the same values
// in float64.

// Rain becomes FZ at 2000 m, and the FZ carried to 1000 m keeps the height: RA
const SBCGolden SBC_CARRIED_FZ{
    saturated_profile(
        {0.0f, 1000.0f, 2000.0f, 3000.0f, 4000.0f},
        {277.0f, 268.0f, 269.0f, 280.0f, 265.0f}),
    4000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    3,
    sharp::PrecipType::rain,
    1.0,
    2000.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},
};

// A cloud top at Tice is frozen, so rule 2 does not fire: FZRA
const SBCGolden SBC_RULE_2_AT_TICE{
    saturated_profile(
        {0.0f, 1000.0f, 2000.0f, 3000.0f},
        {270.15f, 275.15f, 276.15f,
         sharp::SBC_ICE_NUCLEATION_TEMPERATURE}),
    3000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::freezing_rain,
    1.0,
    0.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     0.0f, 0.0f, 0.0f, 0.0f},
};
}  // namespace

TEST_CASE("Testing spectral_bin_classifier core: the sample sounding") {
    check_golden(SBC_SAMPLE);
    const sharp::SpectralBinResult result =
        run_sbc(SBC_SAMPLE.snd, SBC_SAMPLE.cloud_top).result;
    CHECK(result.precip_type == sharp::PrecipType::ice_pellets);
    CHECK(result.liquid_fraction == 0.0f);
    CHECK(result.supercooled_liquid_height == 1130.5469f);
}

TEST_CASE("Testing spectral_bin_classifier core: the Tice switch") {
    check_golden(SBC_TICE_SWITCH);
}

TEST_CASE("Testing spectral_bin_classifier core: refrozen pellets melt again") {
    check_golden(SBC_REMELT);
}

TEST_CASE("Testing spectral_bin_classifier core: a carried FZ class") {
    check_golden(SBC_CARRIED_FZ);
}

TEST_CASE("Testing spectral_bin_classifier core: rule 2 at Tice") {
    check_golden(SBC_RULE_2_AT_TICE);
}

// ---------------------------------------------------------------------------
// Microphysics: liquid cloud tops
// ---------------------------------------------------------------------------

namespace {
// data/sbc_reference case 12
const SBCGolden SBC_SUPERCOOLED_TOP{
    {
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f, 68700.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f},
        {278.15f, 277.15f, 275.15f, 274.15f, 272.15f, 271.15f, 270.15f},
        {278.15f, 277.15f, 275.15f, 274.15f, 272.15f, 271.15f, 270.15f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {278.15f, 277.15f, 275.15f, 274.15f, 272.15f, 271.15f, 270.15f},
    },
    3000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::rain,
    1.0,
    2000.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f},
};

// A warm cloud top reaches the surface without refreezing only through
// levels at exactly 0 C. The reference ran on 273.15 K exactly.

// A warm cloud top over a surface at 0 C: FZRA
const SBCGolden SBC_WARM_TOP_0C_SURFACE{
    saturated_profile({0.0f, 1000.0f, 2000.0f, 3000.0f},
                      {273.15f, 275.15f, 276.15f, 274.15f}),
    3000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    1,
    sharp::PrecipType::freezing_rain,
    1.0,
    0.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f},
};

// A warm cloud top above a level at 0 C: RA
const SBCGolden SBC_WARM_TOP_0C_LEVEL{
    saturated_profile({0.0f, 1000.0f, 2000.0f, 3000.0f},
                      {278.15f, 273.15f, 275.15f, 274.15f}),
    3000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::rain,
    1.0,
    3000.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f},
};
}  // namespace

TEST_CASE("Testing spectral_bin_classifier core: supercooled cloud top") {
    SUBCASE("AGL") { check_golden(SBC_SUPERCOOLED_TOP); }
    SUBCASE("MSL") {
        SBCGolden golden = SBC_SUPERCOOLED_TOP;
        for (float& z : golden.snd.height) z += 1500.0f;
        golden.cloud_top += 1500.0f;
        check_golden(golden);
    }
}

TEST_CASE("Testing spectral_bin_classifier core: warm cloud top") {
    SUBCASE("surface at 0 C") { check_golden(SBC_WARM_TOP_0C_SURFACE); }
    SUBCASE("level at 0 C") { check_golden(SBC_WARM_TOP_0C_LEVEL); }
}

// ---------------------------------------------------------------------------
// Precipitation type from a full sounding
// ---------------------------------------------------------------------------

namespace {
// data/sbc_reference case 13
const SBCGolden SBC_WARM_TOP_REFREEZE{
    {
        {100000.0f, 93900.0f, 88200.0f, 82900.0f, 77900.0f, 73200.0f, 68700.0f},
        {0.0f, 500.0f, 1000.0f, 1500.0f, 2000.0f, 2500.0f, 3000.0f},
        {278.15f, 276.15f, 271.15f, 270.15f, 274.15f, 275.15f, 274.65f},
        {278.15f, 276.15f, 271.15f, 270.15f, 274.15f, 275.15f, 274.65f},
        {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f},
        {278.15f, 276.15f, 271.15f, 270.15f, 274.15f, 275.15f, 274.65f},
    },
    3000.0f,
    &PYTHON_DSD,
    1.0f,
    267.15f,
    2,
    sharp::PrecipType::rain,
    1.0,
    1500.0f,
    {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f,
     1.0f, 1.0f, 1.0f, 1.0f},
};

float sbc_cloud_top(const SBCProfile& snd) {
    return sharp::spectral_bin_cloud_top(
        snd.pressure.data(), snd.height.data(), snd.temperature.data(),
        snd.dewpoint.data(), snd.relh.data(),
        static_cast<std::ptrdiff_t>(snd.height.size()));
}

void check_same(const sharp::SpectralBinResult& result,
                const sharp::SpectralBinResult& expected) {
    CHECK(result.precip_type == expected.precip_type);
    CHECK(result.liquid_fraction == expected.liquid_fraction);
    CHECK(result.supercooled_liquid_height ==
          expected.supercooled_liquid_height);
}

// Asking for the profile must not change the result, and composing
// spectral_bin_cloud_top with the composed overload by hand must give the
// same result and profile.
SBCRun run_full_column(
    const SBCProfile& snd,
    const sharp::SpectralBinDSD& dsd = sharp::spectral_bin_dsd_default(),
    const float tice = sharp::SBC_ICE_NUCLEATION_TEMPERATURE) {
    const auto N = static_cast<std::ptrdiff_t>(snd.height.size());
    SBCRun run;
    run.profile.assign(
        snd.height.size() * static_cast<std::size_t>(dsd.nbins()), 0.0f);
    run.result = sharp::spectral_bin_classifier(
        snd.pressure.data(), snd.height.data(), snd.temperature.data(),
        snd.dewpoint.data(), snd.relh.data(), snd.wetbulb.data(), N, dsd, tice,
        run.profile.data());
    check_same(sharp::spectral_bin_classifier(
                   snd.pressure.data(), snd.height.data(),
                   snd.temperature.data(), snd.dewpoint.data(),
                   snd.relh.data(), snd.wetbulb.data(), N, dsd, tice),
               run.result);
    const SBCRun composed = run_sbc(snd, sbc_cloud_top(snd), dsd, tice);
    check_same(composed.result, run.result);
    CHECK(composed.profile == run.profile);
    return run;
}

// The column's own cloud top is the reference's, so the full column passes
// the golden-data comparison rules that check_golden applies to the
// composed overload.
void check_full_column(const SBCGolden& golden) {
    CHECK(sbc_cloud_top(golden.snd) == golden.cloud_top);
    run_full_column(golden.snd,
                    build_dsd(golden.dsd->diameter, golden.dsd->concentration,
                              golden.rime_factor),
                    golden.tice);
    check_golden(golden);
}
}  // namespace

TEST_CASE("Testing spectral_bin_classifier from a sounding: golden cases") {
    const std::pair<const char*, const SBCGolden*> goldens[] = {
        {"RA", &SBC_CORE_RA},
        {"SN", &SBC_CORE_SN},
        {"RASN", &SBC_CORE_RASN},
        {"dry layer", &SBC_CORE_DRY_LAYER},
        {"C++ 2.0.3 DSD", &SBC_CORE_CXX_DSD},
        {"rime factor 5", &SBC_CORE_RIME_5},
        {"Tice -10 C", &SBC_CORE_TICE_ALT},
        {"melting below the top", &SBC_MELTING_BELOW_TOP},
        {"PL, surface at 0 C", &SBC_ZERO_C_PL},
        {"FZRAPL, surface at 0 C", &SBC_ZERO_C_FZRAPL},
        {"FZRA, surface at 0 C", &SBC_ZERO_C_FZRA},
        {"RAPL", &SBC_ZERO_C_RAPL},
        {"RA, three crossings", &SBC_ZERO_C_RA},
        {"PL in dry air", &SBC_ZERO_C_DRY_PL},
        {"0 C above a crossing", &SBC_ZERO_C_ABOVE_CROSSING},
        {"sample", &SBC_SAMPLE},
        {"Tice switch", &SBC_TICE_SWITCH},
        {"remelt", &SBC_REMELT},
        {"carried FZ", &SBC_CARRIED_FZ},
        {"rule 2 at Tice", &SBC_RULE_2_AT_TICE},
        {"supercooled cloud top", &SBC_SUPERCOOLED_TOP},
        {"warm top, surface at 0 C", &SBC_WARM_TOP_0C_SURFACE},
        {"warm top, level at 0 C", &SBC_WARM_TOP_0C_LEVEL},
        {"warm top refreezes", &SBC_WARM_TOP_REFREEZE},
    };
    for (const auto& [name, golden] : goldens) {
        CAPTURE(name);
        check_full_column(*golden);
    }
}

TEST_CASE("Testing spectral_bin_classifier from a sounding: the sample") {
    const SBCRun run = run_full_column(SBC_SAMPLE.snd);
    CHECK(run.result.precip_type == sharp::PrecipType::ice_pellets);
    CHECK(run.result.liquid_fraction == 0.0f);
    CHECK(run.result.supercooled_liquid_height == 1130.5469f);
}

TEST_CASE("Testing spectral_bin_classifier from a sounding: defaults") {
    const SBCProfile& snd = SBC_SAMPLE.snd;
    const auto N = static_cast<std::ptrdiff_t>(snd.height.size());
    const sharp::SpectralBinDSD dsd = sharp::spectral_bin_dsd_default();
    check_same(sharp::spectral_bin_classifier(
                   snd.pressure.data(), snd.height.data(),
                   snd.temperature.data(), snd.dewpoint.data(),
                   snd.relh.data(), snd.wetbulb.data(), N, dsd),
               sharp::spectral_bin_classifier(
                   snd.pressure.data(), snd.height.data(),
                   snd.temperature.data(), snd.dewpoint.data(),
                   snd.relh.data(), snd.wetbulb.data(), N, dsd,
                   sharp::SBC_ICE_NUCLEATION_TEMPERATURE, nullptr));
}

TEST_CASE("Testing spectral_bin_classifier from a sounding: no cloud") {
    SBCProfile snd = saturated_profile({0.0f, 1000.0f, 2000.0f},
                                       {278.15f, 275.15f, 272.15f});
    snd.dewpoint = {258.15f, 255.15f, 252.15f};
    snd.relh = {0.2f, 0.2f, 0.2f};
    REQUIRE(sbc_cloud_top(snd) == sharp::MISSING);
    check_sbc(run_full_column(snd), SBC_MISSING);
}

TEST_CASE("Testing spectral_bin_classifier from a sounding with few levels") {
    const SBCProfile& snd = SBC_SAMPLE.snd;
    for (const std::size_t N : {0, 1}) {
        CAPTURE(N);
        const SBCProfile few{
            {snd.pressure.begin(), snd.pressure.begin() + N},
            {snd.height.begin(), snd.height.begin() + N},
            {snd.temperature.begin(), snd.temperature.begin() + N},
            {snd.dewpoint.begin(), snd.dewpoint.begin() + N},
            {snd.relh.begin(), snd.relh.begin() + N},
            {snd.wetbulb.begin(), snd.wetbulb.begin() + N},
        };
        check_sbc(run_full_column(few), SBC_MISSING);
    }
}

TEST_CASE("Testing spectral_bin_classifier from a sounding ignores Tw above "
          "the top") {
    SBCGolden golden = SBC_CORE_RA;
    for (const float z : {6000.0f, 7000.0f}) {
        golden.snd.pressure.push_back(100000.0f * std::exp(-z / 8000.0f));
        golden.snd.height.push_back(z);
        golden.snd.temperature.push_back(250.0f);
        golden.snd.dewpoint.push_back(230.0f);
        golden.snd.relh.push_back(0.2f);
        golden.snd.wetbulb.push_back(300.0f);
        golden.profile.insert(golden.profile.end(), 4, sharp::MISSING);
    }
    check_full_column(golden);
}

#ifndef NO_QC
TEST_CASE("Testing spectral_bin_classifier from a sounding: a cloud top "
          "without a wet-bulb temperature") {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    for (const float bad : {sharp::MISSING, NaN}) {
        CAPTURE(bad);
        const SBCProfile snd{
            {100000.0f, 88250.0f, 77880.0f},
            {0.0f, 1000.0f, 2000.0f},
            {280.0f, 280.0f, 280.0f},
            {280.0f, 270.0f, 280.0f},
            {0.9f, 0.5f, 0.9f},
            {279.0f, 279.0f, bad},
        };
        CHECK(sbc_cloud_top(snd) == 2000.0f);
        check_sbc(run_full_column(snd), SBC_RA);
    }
}
#endif
