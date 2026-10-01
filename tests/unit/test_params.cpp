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
        constexpr float hght_2[4] = {0.0f, MISSING, nanval, 1000.0f};
        constexpr float tw_2[4] = {K - 2.0f, K + 9.0f, K + 9.0f, K + 2.0f};
        check_energies(sharp::bourgouin_energy(pres, hght_2, tw_2, 4), triangle,
                       triangle, triangle);
    }
    {
        INFO("one valid level gives MISSING, not zeros");
        constexpr float hght[4] = {0.0f, 250.0f, 750.0f, 1000.0f};
        constexpr float tw[4] = {MISSING, K + 2.0f, nanval, MISSING};
        check_energies_missing(sharp::bourgouin_energy(pres, hght, tw, 4));
        constexpr float hght_2[4] = {nanval, 250.0f, MISSING, MISSING};
        constexpr float tw_2[4] = {K + 2.0f, K + 2.0f, K + 2.0f, K + 2.0f};
        check_energies_missing(sharp::bourgouin_energy(pres, hght_2, tw_2, 4));
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
