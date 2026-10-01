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
