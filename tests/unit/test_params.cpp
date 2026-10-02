#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params/convective.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/winds.h>

#include <cmath>
#include <optional>

#include "doctest.h"

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
    std::optional<float> mw_top_agl;
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
                c.mw_top_agl
                    ? snd.classic({c.base_agl, *c.mw_top_agl}, left, true)
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
             BunkersCase{2000, 7000, std::nullopt, 15.8496647f, -0.509417534f},
             BunkersCase{6000, 14000, 9100, 27.9163494f, -0.846437931f},
         }) {
        check_bunkers(c);
    }
}
