#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <SHARPlib/thermo.h>
#include <SHARPlib/constants.h>

TEST_CASE("Testing theta") {

    constexpr float tmpk = 10.0 + sharp::ZEROCNK;
    constexpr float pres = 900.0 * sharp::HPA_TO_PA;
    CHECK(sharp::theta(pres, tmpk, sharp::THETA_REF_PRESSURE) ==
          doctest::Approx(291.794));

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
    constexpr float tmpk = 10.0 + sharp::ZEROCNK;
    constexpr float theta = 30.0 + sharp::ZEROCNK;

#ifndef NO_QC
    CHECK(sharp::theta_level(sharp::MISSING, tmpk) == sharp::MISSING);
    CHECK(sharp::theta_level(theta, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta_level(sharp::MISSING, sharp::MISSING) == sharp::MISSING);
#endif

}

TEST_CASE("Testing temperature_at_mixratio") {

#ifndef NO_QC
    CHECK(sharp::temperature_at_mixratio(sharp::MISSING, 1000.0) 
                                              == sharp::MISSING);
    CHECK(sharp::temperature_at_mixratio(10.0, sharp::MISSING) 
                                            == sharp::MISSING);
    CHECK(sharp::temperature_at_mixratio(sharp::MISSING, sharp::MISSING) 
                                                      == sharp::MISSING);
#endif

}


TEST_CASE("Testing lcl_temperature") {

#ifndef NO_QC
    CHECK(sharp::lcl_temperature(sharp::MISSING, 10.0) == sharp::MISSING);
    CHECK(sharp::lcl_temperature(10.0, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::lcl_temperature(sharp::MISSING, sharp::MISSING) 
                                              == sharp::MISSING);
#endif

}

TEST_CASE("Testing vapor_pressure") {

#ifndef NO_QC
    CHECK(sharp::vapor_pressure(100000.0f, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::vapor_pressure(sharp::MISSING, sharp::ZEROCNK) ==
          sharp::MISSING);
#endif

}

TEST_CASE("Testing wobf") {

#ifndef NO_QC
    CHECK(sharp::wobf(sharp::MISSING) == sharp::MISSING);
#endif

}
