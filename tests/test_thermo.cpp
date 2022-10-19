#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "thermo.h"

TEST_CASE("Testing the theta (Potential Temperature) function") {

    CHECK(sharp::theta(900.0, 10.0, 1000.0) == doctest::Approx(18.6444));

}

