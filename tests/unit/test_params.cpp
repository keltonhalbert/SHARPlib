#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>

#include <cmath>
#include <limits>

#include "doctest.h"

// ===========================================================================
// Convective wind parameters whose layers convert to MISSING (SHARPlib-mut)
// ===========================================================================
//
// A layer conversion returns {MISSING, MISSING} when the layer extends past
// the profile, when an endpoint's data on its open side are MISSING, or when
// an AGL conversion has no valid height[0]. These routines used that layer
// in arithmetic: they threw std::range_error or returned garbage. QC builds
// now return MISSING (the effective-inflow Bunkers overload falls back to the
// non-parcel method). NO_QC builds don't change.

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

// 0 to 2 km
constexpr std::ptrdiff_t KN = 5;
constexpr float k_pres[KN] = {100000, 95000, 90000, 85000, 80000};
constexpr float k_hght[KN] = {0, 500, 1000, 1500, 2000};
constexpr float k_uwin[KN] = {0, 5, 10, 15, 20};
constexpr float k_vwin[KN] = {0, 2, 4, 6, 8};
// 0 to 8 km
constexpr float d_pres[KN] = {100000, 80000, 62000, 47000, 35000};
constexpr float d_hght[KN] = {0, 2000, 4000, 6000, 8000};
constexpr float d_uwin[KN] = {0, 10, 20, 30, 40};
constexpr float d_vwin[KN] = {0, 2, 4, 6, 8};
constexpr float d_hght_bot[KN] = {M, 2000, 4000, 6000, 8000};
}  // namespace

TEST_CASE("Testing effective_bulk_wind_difference with a MISSING layer") {
    constexpr float hght_top[KN] = {0, 500, 1000, 1500, M};
    // The inflow layer top has no valid height above it
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, hght_top, k_uwin, k_vwin, KN, {100000, 80000}, 95000));
    // The EL has no valid height above it
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, hght_top, k_uwin, k_vwin, KN, {100000, 95000}, 82000));

    // complete data, the inflow layer below the profile
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, k_hght, k_uwin, k_vwin, KN, {105000, 95000}, 90000));
    // complete data, the EL above the profile
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, k_hght, k_uwin, k_vwin, KN, {100000, 95000}, 70000));

    // height[0] MISSING or NaN leaves no ground to measure the AGL layer
    // from (SHARPlib-27b).
    for (const float bad : {M, std::numeric_limits<float>::quiet_NaN()}) {
        CAPTURE(bad);
        const float hght_bad[KN] = {bad, 500, 1000, 1500, 2000};
        check_missing_wind(sharp::effective_bulk_wind_difference(
            k_pres, hght_bad, k_uwin, k_vwin, KN, {95000, 90000}, 85000));
    }

    // inside the profile
    check_wind(sharp::effective_bulk_wind_difference(
                   k_pres, k_hght, k_uwin, k_vwin, KN, {100000, 95000}, 85000),
               7.5f, 3.0f);
}

TEST_CASE("Testing storm_motion_bunkers with a MISSING layer") {
    // the mean wind layer's top has no valid pressure above it
    constexpr float pres_top[KN] = {100000, 80000, 62000, 47000, M};
    check_missing_wind(sharp::storm_motion_bunkers(
        pres_top, d_hght, d_uwin, d_vwin, KN, {0, 7000}, {0, 6000}));
    // the upper end of the shear layer
    check_missing_wind(sharp::storm_motion_bunkers(
        pres_top, d_hght, d_uwin, d_vwin, KN, {0, 6000}, {0, 7000}));
    // height[0] MISSING
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, d_hght_bot, d_uwin, d_vwin, KN, {0, 6000}, {0, 6000}));

    // complete data, layers past a profile that ends at 2 km
    // mean wind layer
    check_missing_wind(sharp::storm_motion_bunkers(
        k_pres, k_hght, k_uwin, k_vwin, KN, {0, 3000}, {0, 2000}));
    // shear layer
    check_missing_wind(sharp::storm_motion_bunkers(
        k_pres, k_hght, k_uwin, k_vwin, KN, {0, 2000}, {0, 3000}));

    // a MISSING layer argument (SHARPlib-yni)
    // shear layer
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, d_hght, d_uwin, d_vwin, KN, {0, 6000}, sharp::HeightLayer()));
    // mean wind layer
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, d_hght, d_uwin, d_vwin, KN, sharp::HeightLayer(), {0, 6000}));

    // inside the profile
    check_wind(sharp::storm_motion_bunkers(k_pres, k_hght, k_uwin, k_vwin, KN,
                                           {0, 2000}, {0, 2000}),
               12.78543f, -2.96357536f);
}

TEST_CASE("Testing effective-inflow storm_motion_bunkers fallback") {
    // An inflow layer below the profile converts to MISSING, so the routine
    // takes its existing fallback, the non-parcel method with 0-6 km layers.
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

    // height[0] MISSING: the fallback is MISSING too.
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, d_hght_bot, d_uwin, d_vwin, KN, {95000, 85000}, mupcl));
    // An EL with no valid height above it
    constexpr float hght_top[KN] = {0, 2000, 4000, 6000, M};
    mupcl.eql_pressure = 40000;
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, hght_top, d_uwin, d_vwin, KN, {100000, 90000}, mupcl));

    // inside the profile
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
    // 1.5 km has no valid pressure above it
    constexpr float pres_top[KN] = {100000, 95000, 90000, M, M};
    check_missing_pair(
        sharp::mcs_motion_corfidi(pres_top, hght, k_uwin, k_vwin, KN));
    // height[0] MISSING
    constexpr float hght_bot[KN] = {M, 500, 1000, 2000, 3000};
    check_missing_pair(
        sharp::mcs_motion_corfidi(pres, hght_bot, k_uwin, k_vwin, KN));
    // complete data, a profile that ends at 1 km
    constexpr float hght_low[KN] = {0, 250, 500, 750, 1000};
    check_missing_pair(
        sharp::mcs_motion_corfidi(k_pres, hght_low, k_uwin, k_vwin, KN));
    // pressure[0] MISSING or NaN (SHARPlib-yni)
    for (const float bad : {M, std::numeric_limits<float>::quiet_NaN()}) {
        CAPTURE(bad);
        const float pres_bad[KN] = {bad, 95000, 90000, 80000, 70000};
        check_missing_pair(
            sharp::mcs_motion_corfidi(pres_bad, hght, k_uwin, k_vwin, KN));
    }

    // inside the profile
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

    // 6 km has no valid pressure above it
    constexpr float pres_top[N] = {100000, 85000, 70000, 59000, 51000, M};
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres_top, hght,
                                      uwin, vwin, N) == M);
    // complete data
    // a MISSING hail growth zone
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, {M, M}, storm, pres, hght,
                                      uwin, vwin, N) == M);
    // a hail growth zone above the profile
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, {65000, 35000}, storm, pres,
                                      hght, uwin, vwin, N) == M);
    // a profile that ends below 6 km
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N - 1) == M);
    // an EL less than 1500 m above the ground
    sharp::Parcel low_el = mu_pcl;
    low_el.eql_pressure = 90000;
    CHECK(sharp::large_hail_parameter(low_el, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N) == M);
    // an EL above the profile top (SHARPlib-yni)
    sharp::Parcel high_el = mu_pcl;
    high_el.eql_pressure = 30000;
    CHECK(sharp::large_hail_parameter(high_el, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N) == M);

    // inside the profile
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin,
                                      N) == doctest::Approx(70.2362289f));
}
#endif

// ===========================================================================
// The effective bulk wind difference and station height (SHARPlib-27b)
// ===========================================================================
//
// effective_bulk_wind_difference built its layer in meters MSL and passed it
// to wind_shear, which takes meters AGL and adds height[0] again. The shear
// came from a layer height[0] meters too high.

namespace {
// 1000 to 500 hPa, 1 km apart, as heights AGL
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

TEST_CASE("Testing effective_bulk_wind_difference ignores station height") {
    // The inflow layer base is the ground and the EL is 4000 m AGL, so the
    // layer is 0 to 2000 m AGL, where u goes from 0 to 12.
    constexpr float vwin[EN] = {0, 0, 0, 0, 0, 0};
    for (const float shift : shifts) {
        CAPTURE(shift);
        const sharp::WindComponents ebwd =
            ebwd_shifted(shift, vwin, {100000, 90000}, 60000);
        CHECK(ebwd.u == doctest::Approx(12.0f));
        CHECK(ebwd.v == doctest::Approx(0.0f));
    }
}

TEST_CASE("Testing effective_bulk_wind_difference against a known value") {
    // From the definition, with heights AGL:
    //   inflow base, 900 hPa:  1000 m
    //   EL, 600 hPa:           4000 m
    //   half the depth:        0.5 * (4000 - 1000) = 1500 m
    //   layer:                 1000 to 2500 m
    //   winds at 1000 m:       u = 10, v = -2 (a level)
    //   winds at 2500 m:       halfway from 2000 to 3000 m, so
    //                          u = (12 + 13) / 2 = 12.5, v = (1 + 4) / 2 = 2.5
    //   EBWD:                  (12.5 - 10, 2.5 - (-2)) = (2.5, 4.5)
    // The station height moves every level by the same amount, so it can't
    // change the answer.
    constexpr float vwin[EN] = {0, -2, 1, 4, 6, 7};
    for (const float shift : shifts) {
        CAPTURE(shift);
        const sharp::WindComponents ebwd =
            ebwd_shifted(shift, vwin, {90000, 80000}, 60000);
        CHECK(ebwd.u == doctest::Approx(2.5f));
        CHECK(ebwd.v == doctest::Approx(4.5f));
    }
}

// ===========================================================================
// The effective-inflow Bunkers mean wind layer (SHARPlib-efz)
// ===========================================================================
//
// Bunkers et al. (2014) take the pressure-weighted mean wind from the
// effective inflow base to 65% of the most-unstable EL height, both in
// meters AGL, with at least 3 km between them. The routine ended the layer
// at 0.65 * (EL - base) instead, and fell back to the 0-6 km method when
// that top was under 3 km or below the base. The two agree when the inflow
// base is the ground.

namespace {
// 0 to 16 km AGL every 500 m, so each inflow base and EL below is a level.
// Pressure falls off with an 8 km scale height, and the hodograph curves.
constexpr std::ptrdiff_t BN = 33;

struct BunkersSounding {
    float pres[BN];
    float hght[BN];
    float uwin[BN];
    float vwin[BN];

    // heights in meters MSL for a station at this elevation
    explicit BunkersSounding(const float elevation) {
        for (std::ptrdiff_t k = 0; k < BN; ++k) {
            const float z = 500.0f * k;  // m AGL
            hght[k] = elevation + z;
            pres[k] = 100000.0f * std::exp(-z / 8000.0f);
            uwin[k] = 30.0f * (1.0f - std::exp(-z / 4000.0f));
            vwin[k] = 10.0f * std::sin(z / 3000.0f);
        }
    }

    // the pressure of the level z meters AGL
    float pres_at(const float z) const {
        return pres[static_cast<std::ptrdiff_t>(z / 500.0f)];
    }

    // the effective-inflow method for an inflow base and an MU EL in m AGL
    sharp::WindComponents effective(const float base, const float el,
                                    const bool left) const {
        sharp::Parcel mupcl;
        mupcl.eql_pressure = pres_at(el);
        const sharp::PressureLayer eil = {pres_at(base),
                                          pres_at(base + 1000.0f)};
        return sharp::storm_motion_bunkers(pres, hght, uwin, vwin, BN, eil,
                                           mupcl, left);
    }

    // the classic method with the 0-6 km AGL shear
    sharp::WindComponents classic(const sharp::HeightLayer mean_wind_agl,
                                  const bool left, const bool weighted) const {
        return sharp::storm_motion_bunkers(pres, hght, uwin, vwin, BN,
                                           mean_wind_agl, {0, 6000}, left,
                                           weighted);
    }
};

struct BunkersCase {
    float base;    // effective inflow base, m AGL
    float el;      // MU EL, m AGL
    float mw_top;  // m AGL; 0 for the 0-6 km fallback
    float u;       // the right mover
    float v;
};

// The routine equals the classic method over {base, mw_top}, pressure
// weighted, or the 0-6 km fallback, and gives the same motion at every
// station elevation.
void check_bunkers(const BunkersCase c) {
    CAPTURE(c.base);
    CAPTURE(c.el);
    for (const float elevation : {0.0f, 1000.0f, 762.3f}) {
        CAPTURE(elevation);
        const BunkersSounding snd(elevation);
        for (const bool left : {false, true}) {
            CAPTURE(left);
            const sharp::WindComponents motion =
                snd.effective(c.base, c.el, left);
            const sharp::WindComponents expected =
                (c.mw_top > 0.0f) ? snd.classic({c.base, c.mw_top}, left, true)
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
             // the inflow base is the ground
             BunkersCase{0, 10000, 6500, 14.8269749f, -1.06337833f},
             BunkersCase{1000, 12000, 7800, 18.9706116f, 0.523952484f},
             BunkersCase{2000, 10000, 6500, 20.7405224f, 1.77062941f},
         }) {
        check_bunkers(c);
    }
}

TEST_CASE("Testing the effective-inflow storm_motion_bunkers 3 km fallback") {
    for (const BunkersCase c : {
             // 2550 m between the base and 0.65 * EL: falls back
             BunkersCase{2000, 7000, 0, 15.8496647f, -0.509417534f},
             // 3100 m: the layer {6000, 9100}
             BunkersCase{6000, 14000, 9100, 27.9163494f, -0.846437931f},
             // 2500 m: falls back
             BunkersCase{4000, 10000, 0, 15.8496647f, -0.509417534f},
             // the inflow base is the ground and 0.65 * EL is 2600 m: falls
             // back
             BunkersCase{0, 4000, 0, 15.8496647f, -0.509417534f},
         }) {
        check_bunkers(c);
    }
}
