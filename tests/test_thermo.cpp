#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "thermo.h"
#include "constants.h"

TEST_CASE("Testing theta") {

    CHECK(sharp::theta(900.0, 10.0, 1000.0) == doctest::Approx(18.6444));
    CHECK(sharp::theta(sharp::MISSING, 10.0, 1000.0) == sharp::MISSING);
    CHECK(sharp::theta(900.0, sharp::MISSING, 1000.0) == sharp::MISSING);
    CHECK(sharp::theta(900.0, 10.0, sharp::MISSING) == sharp::MISSING);

    CHECK(sharp::theta(sharp::MISSING, sharp::MISSING, 1000.0) == sharp::MISSING);
    CHECK(sharp::theta(sharp::MISSING, 10.0, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta(900.0, sharp::MISSING, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta(sharp::MISSING, sharp::MISSING, sharp::MISSING) 
                                                    == sharp::MISSING);
}

TEST_CASE("Testing theta_level") {

    CHECK(sharp::theta_level(sharp::MISSING, 10.0) == sharp::MISSING);
    CHECK(sharp::theta_level(30.0, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::theta_level(sharp::MISSING, sharp::MISSING) == sharp::MISSING);

}

TEST_CASE("Testing temperature_at_mixratio") {
    CHECK(sharp::temperature_at_mixratio(sharp::MISSING, 1000.0) 
                                              == sharp::MISSING);
    CHECK(sharp::temperature_at_mixratio(10.0, sharp::MISSING) 
                                            == sharp::MISSING);
    CHECK(sharp::temperature_at_mixratio(sharp::MISSING, sharp::MISSING) 
                                                      == sharp::MISSING);
}


TEST_CASE("Testing lcl_temperature") {
    CHECK(sharp::lcl_temperature(sharp::MISSING, 10.0) == sharp::MISSING);
    CHECK(sharp::lcl_temperature(10.0, sharp::MISSING) == sharp::MISSING);
    CHECK(sharp::lcl_temperature(sharp::MISSING, sharp::MISSING) 
                                              == sharp::MISSING);
}

TEST_CASE("Testing vapor_pressure") {
    CHECK(sharp::vapor_pressure(sharp::MISSING) == sharp::MISSING);
}

TEST_CASE("Testing wobf") {
    CHECK(sharp::wobf(sharp::MISSING) == sharp::MISSING);
}
