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
