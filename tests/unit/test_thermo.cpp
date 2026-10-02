#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>

#include <limits>

#include "doctest.h"

constexpr float nanval = std::numeric_limits<float>::quiet_NaN();

TEST_CASE("Testing theta") {
    constexpr float tmpk = 298.0f;
    constexpr float pres = 10000.0f;
    constexpr float exptexted_theta = 575.348;
    CHECK(sharp::theta(pres, tmpk) == doctest::Approx(exptexted_theta));

#ifndef NO_QC
    CHECK(sharp::theta(sharp::MISSING, tmpk, sharp::THETA_REF_PRESSURE) ==
          sharp::MISSING);
    CHECK(sharp::theta(pres, sharp::MISSING, sharp::THETA_REF_PRESSURE) ==
          sharp::MISSING);
    CHECK(sharp::theta(pres, tmpk, sharp::MISSING) == sharp::MISSING);

    CHECK(sharp::theta(sharp::MISSING, sharp::MISSING,
                       sharp::THETA_REF_PRESSURE) == sharp::MISSING);
    CHECK(sharp::theta(sharp::MISSING, tmpk, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta(pres, sharp::MISSING, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta(sharp::MISSING, sharp::MISSING, sharp::MISSING) ==
          sharp::MISSING);
#endif
}

TEST_CASE("Testing theta_level") {
#ifndef NO_QC
    constexpr float tmpk = 10.0 + sharp::ZEROCNK;
    constexpr float theta = 30.0 + sharp::ZEROCNK;

    CHECK(sharp::theta_level(sharp::MISSING, tmpk) == sharp::MISSING);
    CHECK(sharp::theta_level(theta, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta_level(sharp::MISSING, sharp::MISSING) == sharp::MISSING);
#endif
}

TEST_CASE("Testing temperature_at_mixratio") {
#ifndef NO_QC
    CHECK(sharp::temperature_at_mixratio(sharp::MISSING, 1000.0) ==
          sharp::MISSING);
    CHECK(sharp::temperature_at_mixratio(10.0, sharp::MISSING) ==
          sharp::MISSING);
    CHECK(sharp::temperature_at_mixratio(sharp::MISSING, sharp::MISSING) ==
          sharp::MISSING);
#endif
}

TEST_CASE("Testing lcl temperature and pressure") {
#ifndef NO_QC
    CHECK(sharp::lcl_temperature(sharp::MISSING, 10.0) == sharp::MISSING);
    CHECK(sharp::lcl_temperature(10.0, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::lcl_temperature(sharp::MISSING, sharp::MISSING) ==
          sharp::MISSING);
#endif

    static constexpr float pres = 101716.0f;
    static constexpr float tmpk = 273.09723f;
    static constexpr float dwpk = 264.5351f;
    static constexpr float expected_lcl_tmpk = 262.818f;
    static constexpr float expected_lcl_pres = 88934.5f;

    float lcl_temperature, lcl_pressure;
    sharp::drylift(pres, tmpk, dwpk, lcl_pressure, lcl_temperature);

    CHECK(lcl_temperature == doctest::Approx(expected_lcl_tmpk));
    CHECK(lcl_pressure == doctest::Approx(expected_lcl_pres));
}

TEST_CASE("Testing vapor_pressure") {
#ifndef NO_QC
    CHECK(sharp::vapor_pressure(100000.0f, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::vapor_pressure(sharp::MISSING, sharp::ZEROCNK) ==
          sharp::MISSING);
#endif

    static constexpr float pres = 100000.0f;
    static constexpr float dwpk = 25.0f + sharp::ZEROCNK;
    static constexpr float expected_es = 3167.0f;
    static constexpr float percent_tol = 0.0005f;
    CHECK(sharp::vapor_pressure(pres, dwpk) ==
          doctest::Approx(expected_es).epsilon(percent_tol));
}

TEST_CASE("Testing relative_humidity_ice") {
#ifndef NO_QC
    // MISSING in any argument gives MISSING, as in relative_humidity
    CHECK(sharp::relative_humidity(sharp::MISSING, 250.0f, 245.0f) ==
          sharp::MISSING);
    CHECK(sharp::relative_humidity_ice(sharp::MISSING, 250.0f, 245.0f) ==
          sharp::MISSING);
    CHECK(sharp::relative_humidity(50000.0f, sharp::MISSING, 245.0f) ==
          sharp::MISSING);
    CHECK(sharp::relative_humidity_ice(50000.0f, sharp::MISSING, 245.0f) ==
          sharp::MISSING);
    CHECK(sharp::relative_humidity(50000.0f, 250.0f, sharp::MISSING) ==
          sharp::MISSING);
    CHECK(sharp::relative_humidity_ice(50000.0f, 250.0f, sharp::MISSING) ==
          sharp::MISSING);
#endif

    // Expected values computed in float64 with numpy from Bolton (1980)
    // eq. 10 over liquid, 611.2 exp(17.67 t / (t + 243.5)), and
    // 611.2 exp(21.8745584 t / (t + 265.49)) over ice (t in C, Pa),
    // each floored at half the air pressure.
    // Saturated with respect to liquid at -10 C: supersaturated over ice.
    CHECK(sharp::relative_humidity_ice(80000.0f, 263.15f, 263.15f) ==
          doctest::Approx(1.1045472f));
    CHECK(sharp::relative_humidity_ice(50000.0f, 250.0f, 245.0f) ==
          doctest::Approx(0.8023845f));
    CHECK(sharp::relative_humidity_ice(70000.0f, 268.15f, 260.15f) ==
          doctest::Approx(0.5617494f));
    // At 300 Pa the ice saturation vapor pressure (195.5 Pa) is floored to
    // 150 Pa, which gives 0.6366 instead of 0.4885.
    CHECK(sharp::relative_humidity_ice(300.0f, 260.0f, 250.0f) ==
          doctest::Approx(0.6365938f));

    // Equal to relative humidity over liquid at exactly 0 C
    for (const float dwpk : {sharp::ZEROCNK, 268.15f, 255.0f}) {
        CHECK(sharp::relative_humidity_ice(100000.0f, sharp::ZEROCNK, dwpk) ==
              sharp::relative_humidity(100000.0f, sharp::ZEROCNK, dwpk));
    }
}

TEST_CASE("Testing wobf") {
#ifndef NO_QC
    CHECK(sharp::wobf(sharp::MISSING) == sharp::MISSING);
#endif
}

TEST_CASE("Testing Equivalent Potential Temperature") {
    static constexpr float pres = 100000.0f;
    static constexpr float tmpk = 293.0f;
    static constexpr float dwpk = 280.0f;
    static constexpr float expected_thetae = 312.385;

    const float theta_e = sharp::thetae(pres, tmpk, dwpk);
    CHECK(theta_e == doctest::Approx(expected_thetae));
}

TEST_CASE("Testing Wetbulb Temperature") {
    static constexpr sharp::lifter_wobus wobf;
    static sharp::lifter_cm1 cm1;

    cm1.ma_type = sharp::adiabat::adiab_liq;

    static constexpr float pres = 90000.0f;
    static constexpr float tmpk = 20.0f + sharp::ZEROCNK;
    static constexpr float dwpk = 10.0f + sharp::ZEROCNK;

    float wblbk = sharp::wetbulb(wobf, pres, tmpk, dwpk);
    printf("TD: %f\tTW: %f\tTA: %f\n", dwpk, wblbk, tmpk);

    wblbk = sharp::wetbulb(cm1, pres, tmpk, dwpk);
    printf("TD: %f\tTW: %f\tTA: %f\n", dwpk, wblbk, tmpk);
}

constexpr std::ptrdiff_t LR_N = 3;
constexpr float lr_pres[LR_N] = {100000, 95000, 90000};
constexpr float lr_tmpk[LR_N] = {300, 297, 294};
constexpr float lr_hght_0[LR_N] = {0, 500, 1000};
constexpr float lr_hght_300[LR_N] = {300, 800, 1300};

TEST_CASE("Testing lapse_rate over layers outside the profile") {
    for (const float* hght : {lr_hght_0, lr_hght_300}) {
        CAPTURE(hght[0]);
        CHECK(sharp::lapse_rate(sharp::PressureLayer(85000, 80000), lr_pres,
                                hght, lr_tmpk, LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::PressureLayer(110000, 105000), lr_pres,
                                hght, lr_tmpk, LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::HeightLayer(1500, 2000), hght, lr_tmpk,
                                LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::HeightLayer(-500, -100), hght, lr_tmpk,
                                LR_N) == sharp::MISSING);
    }

    constexpr float pres[1] = {100000};
    constexpr float hght[1] = {0};
    constexpr float tmpk[1] = {300};
    CHECK(sharp::lapse_rate(sharp::PressureLayer(85000, 80000), pres, hght,
                            tmpk, 1) == sharp::MISSING);
    CHECK(sharp::lapse_rate(sharp::PressureLayer(110000, 105000), pres, hght,
                            tmpk, 1) == sharp::MISSING);
}

template <typename T>
static void check_lr_max(const sharp::PressureLayer search, const float depth,
                         const float pres[], const float hght[],
                         const float tmpk[], const std::ptrdiff_t N, const T lr,
                         const float bottom, const float top) {
    INFO("search ", search.bottom, " to ", search.top, ", depth ", depth);
    sharp::PressureLayer max_lyr = {0, 0};
    CHECK(sharp::lapse_rate_max(search, depth, pres, hght, tmpk, N, &max_lyr) ==
          lr);
    CHECK(max_lyr.bottom == bottom);
    CHECK(max_lyr.top == top);
}

template <typename T>
static void check_lr_max(const sharp::HeightLayer search, const float depth,
                         const float hght[], const float tmpk[],
                         const std::ptrdiff_t N, const T lr, const float bottom,
                         const float top) {
    INFO("search ", search.bottom, " to ", search.top, ", depth ", depth);
    sharp::HeightLayer max_lyr = {0, 0};
    CHECK(sharp::lapse_rate_max(search, depth, hght, tmpk, N, &max_lyr) == lr);
    CHECK(max_lyr.bottom == bottom);
    CHECK(max_lyr.top == top);
}

constexpr std::ptrdiff_t LRM_N = 5;
constexpr float top_pres[LRM_N] = {100000, 90000, 80000, 77500, 75000};
constexpr float top_hght_0[LRM_N] = {0, 1000, 2000, 2250, 2500};
constexpr float top_hght_300[LRM_N] = {300, 1300, 2300, 2550, 2800};
constexpr float top_tmpk[LRM_N] = {310, 298, 292, 290, 289};
constexpr float top_tmpk_steep[LRM_N] = {310, 298, 292, 290, 286};
constexpr float bot_pres[LRM_N] = {85000, 80000, 70000, 60000, 50000};
constexpr float bot_hght[LRM_N] = {1500, 2000, 3100, 4300, 5700};
constexpr float bot_tmpk[LRM_N] = {295, 292, 284, 276, 266};
constexpr float bot_tmpk_steep[LRM_N] = {295, 291, 284, 276, 266};

TEST_CASE("Testing lapse_rate_max over a layer that leaves the profile") {
    for (const float* hght : {top_hght_0, top_hght_300}) {
        CAPTURE(hght[0]);
        check_lr_max(sharp::PressureLayer(80000, 60000), 5000, top_pres, hght,
                     top_tmpk, LRM_N, doctest::Approx(6.0f), 80000, 75000);

        sharp::HeightLayer max_hlyr = {0, 0};
        CHECK(sharp::lapse_rate_max(sharp::HeightLayer(2000, 6000), 500, hght,
                                    top_tmpk, LRM_N,
                                    &max_hlyr) == doctest::Approx(6.0f));
    }

    check_lr_max(sharp::PressureLayer(100000, 50000), 10000, bot_pres, bot_hght,
                 bot_tmpk, LRM_N, doctest::Approx(80.0f / 11.0f), 80000, 70000);
    CHECK(sharp::lapse_rate_max(sharp::PressureLayer(100000, 90000), 5000,
                                bot_pres, bot_hght, bot_tmpk,
                                LRM_N) == sharp::MISSING);
}

TEST_CASE("Testing lapse_rate_max at the bottom of the profile") {
    check_lr_max(sharp::PressureLayer(100000, 50000), 10000, bot_pres, bot_hght,
                 bot_tmpk_steep, LRM_N, doctest::Approx(7.15671f), 85000,
                 75000);
    check_lr_max(sharp::PressureLayer(85000, 50000), 10000, bot_pres, bot_hght,
                 bot_tmpk_steep, LRM_N, doctest::Approx(7.15671f), 85000,
                 75000);
    check_lr_max(sharp::PressureLayer(100000, 70000, -2000), 10000, bot_pres,
                 bot_hght, bot_tmpk_steep, LRM_N, doctest::Approx(6.99398f),
                 84000, 74000);
    check_lr_max(sharp::PressureLayer(100000, 40000), 40000, bot_pres, bot_hght,
                 bot_tmpk_steep, LRM_N, sharp::MISSING, sharp::MISSING,
                 sharp::MISSING);

    check_lr_max(sharp::HeightLayer(-1500, 4000), 1000, bot_hght,
                 bot_tmpk_steep, LRM_N, doctest::Approx(7.18182f), 0, 1000);
    check_lr_max(sharp::HeightLayer(0, 4200), 1000, bot_hght, bot_tmpk_steep,
                 LRM_N, doctest::Approx(7.18182f), 0, 1000);
    check_lr_max(sharp::HeightLayer(-1450, 2000), 1000, bot_hght,
                 bot_tmpk_steep, LRM_N, doctest::Approx(7.1f), 50, 1050);
}

TEST_CASE("Testing lapse_rate_max at the top of the profile") {
    for (const float* hght : {top_hght_0, top_hght_300}) {
        CAPTURE(hght[0]);
        check_lr_max(sharp::PressureLayer(80000, 60000), 5000, top_pres, hght,
                     top_tmpk_steep, LRM_N, doctest::Approx(12.0f), 80000,
                     75000);
        check_lr_max(sharp::PressureLayer(80000, 60000), 10000, top_pres, hght,
                     top_tmpk_steep, LRM_N, sharp::MISSING, sharp::MISSING,
                     sharp::MISSING);

        check_lr_max(sharp::HeightLayer(2000, 6000), 500, hght, top_tmpk_steep,
                     LRM_N, doctest::Approx(12.0f), 2000, 2500);
        check_lr_max(sharp::HeightLayer(2000, 6000), 1000, hght, top_tmpk_steep,
                     LRM_N, sharp::MISSING, sharp::MISSING, sharp::MISSING);
        CHECK(sharp::lapse_rate_max(sharp::HeightLayer(0, 6000), 3000, hght,
                                    top_tmpk_steep, LRM_N) == sharp::MISSING);
    }
}

TEST_CASE("Testing lapse_rate_max with a delta that doesn't step upward") {
    for (const float delta : {0.0f, -100.0f, nanval}) {
        CAPTURE(delta);
        const sharp::HeightLayer search(0, 1000, delta);
        check_lr_max(search, 500, lr_hght_0, lr_tmpk, LR_N, sharp::MISSING,
                     sharp::MISSING, sharp::MISSING);
        CHECK(sharp::lapse_rate_max(search, 500, lr_hght_0, lr_tmpk, LR_N) ==
              sharp::MISSING);
    }

    for (const float delta : {0.0f, 1000.0f, nanval}) {
        CAPTURE(delta);
        const sharp::PressureLayer search(100000, 90000, delta);
        check_lr_max(search, 5000, lr_pres, lr_hght_0, lr_tmpk, LR_N,
                     sharp::MISSING, sharp::MISSING, sharp::MISSING);
        CHECK(sharp::lapse_rate_max(search, 5000, lr_pres, lr_hght_0, lr_tmpk,
                                    LR_N) == sharp::MISSING);
    }
}
