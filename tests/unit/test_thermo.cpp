#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>

#include "doctest.h"

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

// lapse_rate over layers that leave the profile (SHARPlib-103). A layer wholly
// outside the profile has no lapse rate and returns MISSING in every build.
// Each "was" comment is the output measured before this change, in both QC
// and NO_QC builds. The pressure overload clipped one end of such a layer and
// not the other, which inverted it, and the conversion to height then threw
// std::range_error or gave the lapse rate of a different layer.
constexpr std::ptrdiff_t LR_N = 3;
constexpr float lr_pres[LR_N] = {100000, 95000, 90000};
constexpr float lr_tmpk[LR_N] = {300, 297, 294};
// the same profile with the surface at 0 m and at 300 m
constexpr float lr_hght_0[LR_N] = {0, 500, 1000};
constexpr float lr_hght_300[LR_N] = {300, 800, 1300};

TEST_CASE("Testing lapse_rate over layers outside the profile") {
    for (const float* hght : {lr_hght_0, lr_hght_300}) {
        CAPTURE(hght[0]);
        // wholly above: was std::range_error (0 m), 6 K/km (300 m)
        CHECK(sharp::lapse_rate(sharp::PressureLayer(85000, 80000), lr_pres,
                                hght, lr_tmpk, LR_N) == sharp::MISSING);
        // wholly below: was std::range_error
        CHECK(sharp::lapse_rate(sharp::PressureLayer(110000, 105000), lr_pres,
                                hght, lr_tmpk, LR_N) == sharp::MISSING);
        // height layers, unchanged
        CHECK(sharp::lapse_rate(sharp::HeightLayer(1500, 2000), hght, lr_tmpk,
                                LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::HeightLayer(-500, -100), hght, lr_tmpk,
                                LR_N) == sharp::MISSING);
    }

    // one level at 0 m, as found by the SHARPlib-cld sanitizer sweep
    constexpr float pres[1] = {100000};
    constexpr float hght[1] = {0};
    constexpr float tmpk[1] = {300};
    // was std::range_error
    CHECK(sharp::lapse_rate(sharp::PressureLayer(85000, 80000), pres, hght,
                            tmpk, 1) == sharp::MISSING);
    // was std::range_error
    CHECK(sharp::lapse_rate(sharp::PressureLayer(110000, 105000), pres, hght,
                            tmpk, 1) == sharp::MISSING);
}

TEST_CASE("Testing lapse_rate over layers at the profile edge") {
    // Layers that touch the profile at one point have no depth and return
    // MISSING. Layers that partly overlap it are clipped to it. Both are
    // unchanged.
    for (const float* hght : {lr_hght_0, lr_hght_300}) {
        CAPTURE(hght[0]);
        // touching at one point
        CHECK(sharp::lapse_rate(sharp::PressureLayer(90000, 80000), lr_pres,
                                hght, lr_tmpk, LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::PressureLayer(105000, 100000), lr_pres,
                                hght, lr_tmpk, LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::HeightLayer(1000, 2000), hght, lr_tmpk,
                                LR_N) == sharp::MISSING);
        CHECK(sharp::lapse_rate(sharp::HeightLayer(-500, 0), hght, lr_tmpk,
                                LR_N) == sharp::MISSING);

        // partly above and partly below
        CHECK(sharp::lapse_rate(sharp::PressureLayer(95000, 80000), lr_pres,
                                hght, lr_tmpk, LR_N) == doctest::Approx(6.0f));
        CHECK(sharp::lapse_rate(sharp::PressureLayer(105000, 95000), lr_pres,
                                hght, lr_tmpk, LR_N) == doctest::Approx(6.0f));
        CHECK(sharp::lapse_rate(sharp::HeightLayer(500, 1500), hght, lr_tmpk,
                                LR_N) == doctest::Approx(6.0f));
        CHECK(sharp::lapse_rate(sharp::HeightLayer(-500, 500), hght, lr_tmpk,
                                LR_N) == doctest::Approx(6.0f));
    }
}

TEST_CASE("Testing lapse_rate_max over a layer that leaves the profile") {
    // The sounding ends at 750 hPa, inside the 800-600 hPa search layer. The
    // lowest kilometre is superadiabatic, so the lapse rate of the whole
    // profile, 8.4 K/km, is larger than any in the search layer. The largest
    // there is 6 K/km, from 800 to 750 hPa.
    constexpr std::ptrdiff_t N = 5;
    constexpr float pres[N] = {100000, 90000, 80000, 77500, 75000};
    constexpr float tmpk[N] = {310, 298, 292, 290, 289};
    constexpr float hght_0[N] = {0, 1000, 2000, 2250, 2500};
    constexpr float hght_300[N] = {300, 1300, 2300, 2550, 2800};

    for (const float* hght : {hght_0, hght_300}) {
        CAPTURE(hght[0]);
        // was std::range_error (0 m); 8.4 K/km over 74000-69000 Pa, a layer
        // wholly above the profile (300 m)
        sharp::PressureLayer max_plyr = {0, 0};
        CHECK(sharp::lapse_rate_max(sharp::PressureLayer(80000, 60000), 5000,
                                    pres, hght, tmpk, N,
                                    &max_plyr) == doctest::Approx(6.0f));
        CHECK(max_plyr.bottom == 80000);
        CHECK(max_plyr.top == 75000);

        // the height search was already right, unchanged
        sharp::HeightLayer max_hlyr = {0, 0};
        CHECK(sharp::lapse_rate_max(sharp::HeightLayer(2000, 6000), 500, hght,
                                    tmpk, N,
                                    &max_hlyr) == doctest::Approx(6.0f));
    }

    // A sounding whose surface, at 850 hPa, is above the 1000 hPa bottom of
    // the search layer. The largest lapse rate, 7.27 K/km, is from 800 to
    // 700 hPa.
    constexpr std::ptrdiff_t N_HI = 5;
    constexpr float pres_hi[N_HI] = {85000, 80000, 70000, 60000, 50000};
    constexpr float hght_hi[N_HI] = {1500, 2000, 3100, 4300, 5700};
    constexpr float tmpk_hi[N_HI] = {295, 292, 284, 276, 266};
    // was std::range_error
    sharp::PressureLayer max_plyr = {0, 0};
    CHECK(sharp::lapse_rate_max(sharp::PressureLayer(100000, 50000), 10000,
                                pres_hi, hght_hi, tmpk_hi, N_HI,
                                &max_plyr) == doctest::Approx(80.0f / 11.0f));
    CHECK(max_plyr.bottom == 80000);
    CHECK(max_plyr.top == 70000);
    // a search layer wholly below the surface; was std::range_error
    CHECK(sharp::lapse_rate_max(sharp::PressureLayer(100000, 90000), 5000,
                                pres_hi, hght_hi, tmpk_hi,
                                N_HI) == sharp::MISSING);
}
