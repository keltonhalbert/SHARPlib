#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>

#include <cmath>
#include <stdexcept>

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
// non-parcel method). Each "was" comment is the output measured before this
// change. NO_QC builds don't change.

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
}  // namespace

TEST_CASE("Testing effective_bulk_wind_difference with a MISSING layer") {
    constexpr float hght_top[KN] = {0, 500, 1000, 1500, M};
    // The inflow layer top has no valid height above it: was std::range_error
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, hght_top, k_uwin, k_vwin, KN, {100000, 80000}, 95000));
    // The EL has no valid height above it: was std::range_error
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, hght_top, k_uwin, k_vwin, KN, {100000, 95000}, 82000));

    // complete data, the inflow layer below the profile: was std::range_error
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, k_hght, k_uwin, k_vwin, KN, {105000, 95000}, 90000));
    // complete data, the EL above the profile: was std::range_error
    check_missing_wind(sharp::effective_bulk_wind_difference(
        k_pres, k_hght, k_uwin, k_vwin, KN, {100000, 95000}, 70000));

    // height[0] MISSING or NaN leaves no ground to measure the AGL layer
    // from (SHARPlib-27b). MISSING was (0.238117903, 0.0952471644). NaN was
    // already MISSING.
    for (const float bad : {M, std::numeric_limits<float>::quiet_NaN()}) {
        CAPTURE(bad);
        const float hght_bad[KN] = {bad, 500, 1000, 1500, 2000};
        check_missing_wind(sharp::effective_bulk_wind_difference(
            k_pres, hght_bad, k_uwin, k_vwin, KN, {95000, 90000}, 85000));
    }

    // inside the profile, unchanged
    check_wind(sharp::effective_bulk_wind_difference(
                   k_pres, k_hght, k_uwin, k_vwin, KN, {100000, 95000}, 85000),
               7.5f, 3.0f);
}

TEST_CASE("Testing storm_motion_bunkers with a MISSING layer") {
    // the mean wind layer's top has no valid pressure above it: was
    // std::range_error
    constexpr float pres_top[KN] = {100000, 80000, 62000, 47000, M};
    check_missing_wind(sharp::storm_motion_bunkers(
        pres_top, d_hght, d_uwin, d_vwin, KN, {0, 7000}, {0, 6000}));
    // the upper end of the shear layer: was (8.75354576, 8.11486626)
    check_missing_wind(sharp::storm_motion_bunkers(
        pres_top, d_hght, d_uwin, d_vwin, KN, {0, 6000}, {0, 7000}));
    // height[0] MISSING: was std::range_error
    constexpr float hght_bot[KN] = {M, 2000, 4000, 6000, 8000};
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, hght_bot, d_uwin, d_vwin, KN, {0, 6000}, {0, 6000}));

    // complete data, layers past a profile that ends at 2 km
    // mean wind layer: was (-9996.21484, -10005.9639)
    check_missing_wind(sharp::storm_motion_bunkers(
        k_pres, k_hght, k_uwin, k_vwin, KN, {0, 3000}, {0, 2000}));
    // shear layer: was (4.69709682, 9.30369854)
    check_missing_wind(sharp::storm_motion_bunkers(
        k_pres, k_hght, k_uwin, k_vwin, KN, {0, 2000}, {0, 3000}));

    // inside the profile, unchanged
    check_wind(sharp::storm_motion_bunkers(k_pres, k_hght, k_uwin, k_vwin, KN,
                                           {0, 2000}, {0, 2000}),
               12.78543f, -2.96357536f);
}

TEST_CASE("Testing effective-inflow storm_motion_bunkers fallback") {
    // An inflow layer below the profile converts to MISSING, so the routine
    // takes its existing fallback, the non-parcel method with 0-6 km layers.
    // Both were std::range_error.
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

    // height[0] MISSING: the fallback is MISSING too. Was std::range_error.
    constexpr float hght_bot[KN] = {M, 2000, 4000, 6000, 8000};
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, hght_bot, d_uwin, d_vwin, KN, {95000, 85000}, mupcl));
    // An EL with no valid height above it: was (NaN, NaN)
    constexpr float hght_top[KN] = {0, 2000, 4000, 6000, M};
    mupcl.eql_pressure = 40000;
    check_missing_wind(sharp::storm_motion_bunkers(
        d_pres, hght_top, d_uwin, d_vwin, KN, {100000, 90000}, mupcl));

    // inside the profile, unchanged
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
    // 1.5 km has no valid pressure above it: was std::range_error
    constexpr float pres_top[KN] = {100000, 95000, 90000, M, M};
    check_missing_pair(
        sharp::mcs_motion_corfidi(pres_top, hght, k_uwin, k_vwin, KN));
    // height[0] MISSING: was std::range_error
    constexpr float hght_bot[KN] = {M, 500, 1000, 2000, 3000};
    check_missing_pair(
        sharp::mcs_motion_corfidi(pres, hght_bot, k_uwin, k_vwin, KN));
    // complete data, a profile that ends at 1 km: was (10016.5, 10006) and
    // (10034, 10013)
    constexpr float hght_low[KN] = {0, 250, 500, 750, 1000};
    check_missing_pair(
        sharp::mcs_motion_corfidi(k_pres, hght_low, k_uwin, k_vwin, KN));

    // inside the profile, unchanged
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

    // 6 km has no valid pressure above it: was std::range_error
    constexpr float pres_top[N] = {100000, 85000, 70000, 59000, 51000, M};
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres_top, hght,
                                      uwin, vwin, N) == M);
    // complete data
    // a MISSING hail growth zone: was 117.714905
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, {M, M}, storm, pres, hght,
                                      uwin, vwin, N) == M);
    // a hail growth zone above the profile: was 117.714905
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, {65000, 35000}, storm, pres,
                                      hght, uwin, vwin, N) == M);
    // a profile that ends below 6 km: was 0
    CHECK(sharp::large_hail_parameter(mu_pcl, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N - 1) == M);
    // an EL less than 1500 m above the ground: was 0
    sharp::Parcel low_el = mu_pcl;
    low_el.eql_pressure = 90000;
    CHECK(sharp::large_hail_parameter(low_el, 8.0f, hgz, storm, pres, hght,
                                      uwin, vwin, N) == M);

    // inside the profile, unchanged
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
// came from a layer height[0] meters too high. Each was_ value is the output
// measured before the fix, in QC and NO_QC builds alike.

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
}  // namespace

TEST_CASE("Testing effective_bulk_wind_difference ignores station height") {
    // The inflow layer base is the ground and the EL is 4000 m AGL, so the
    // layer is 0 to 2000 m AGL, where u goes from 0 to 12.
    constexpr float vwin[EN] = {0, 0, 0, 0, 0, 0};
    struct Case {
        float shift;
        float was_u;
    };
    for (const Case c :
         {Case{0.0f, 12.0f}, Case{1000.0f, 3.0f}, Case{1234.5f, 2.53100014f},
          Case{762.3f, 5.13929987f}}) {
        CAPTURE(c.shift);
        CAPTURE(c.was_u);
        const sharp::WindComponents ebwd =
            ebwd_shifted(c.shift, vwin, {100000, 90000}, 60000);
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
    struct Case {
        float shift;
        float was_u;
        float was_v;
    };
    for (const Case c : {Case{0.0f, 2.5f, 4.5f}, Case{1000.0f, 1.0f, 4.0f},
                         Case{1234.5f, 0.765500069f, 3.76549983f},
                         Case{762.3f, 1.47539997f, 4.23769951f}}) {
        CAPTURE(c.shift);
        CAPTURE(c.was_u);
        CAPTURE(c.was_v);
        const sharp::WindComponents ebwd =
            ebwd_shifted(c.shift, vwin, {90000, 80000}, 60000);
        CHECK(ebwd.u == doctest::Approx(2.5f));
        CHECK(ebwd.v == doctest::Approx(4.5f));
    }
}
