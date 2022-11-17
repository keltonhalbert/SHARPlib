#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "thermo.h"
#include "constants.h"

TEST_CASE("Testing theta") {

    CHECK(sharp::theta(900.0, 10.0, 1000.0) == doctest::Approx(18.6533));

#ifndef NO_QC
    CHECK(sharp::theta(sharp::MISSING, 10.0, 1000.0) == sharp::MISSING);
    CHECK(sharp::theta(900.0, sharp::MISSING, 1000.0) == sharp::MISSING);
    CHECK(sharp::theta(900.0, 10.0, sharp::MISSING) == sharp::MISSING);

    CHECK(sharp::theta(sharp::MISSING, sharp::MISSING, 1000.0) == sharp::MISSING);
    CHECK(sharp::theta(sharp::MISSING, 10.0, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta(900.0, sharp::MISSING, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta(sharp::MISSING, sharp::MISSING, sharp::MISSING) 
                                                    == sharp::MISSING);
#endif
}

TEST_CASE("Testing theta_level") {

#ifndef NO_QC
    CHECK(sharp::theta_level(sharp::MISSING, 10.0) == sharp::MISSING);
    CHECK(sharp::theta_level(30.0, sharp::MISSING) == sharp::MISSING);
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
    CHECK(sharp::vapor_pressure(sharp::MISSING) == sharp::MISSING);
#endif

}

TEST_CASE("Testing wobf") {

#ifndef NO_QC
    CHECK(sharp::wobf(sharp::MISSING) == sharp::MISSING);
#endif

}
