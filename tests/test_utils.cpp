#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <stdexcept>
#include <cmath>

#include "doctest.h"
#include "utils.h"

TEST_CASE("Testing HeightLayer structs") {
    CHECK_THROWS_AS(
        sharp::HeightLayer layer2( 1000.0, 0.0 ), const std::range_error&
        );
}

TEST_CASE("Testing PressureLayer structs") {
    CHECK_THROWS_AS(
        sharp::PressureLayer layer2( 100.0, 1000.0 ), const std::range_error&
        );
}
