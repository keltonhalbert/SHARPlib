#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/qc.h>

#include <limits>

#include "doctest.h"

TEST_CASE("Testing is_missing") {
    CHECK(sharp::is_missing(sharp::MISSING));
    CHECK(sharp::is_missing(std::numeric_limits<float>::quiet_NaN()));

    CHECK_FALSE(sharp::is_missing(0.0f));
    CHECK_FALSE(sharp::is_missing(sharp::ZEROCNK));
    CHECK_FALSE(sharp::is_missing(sharp::MISSING + 1.0f));
}
