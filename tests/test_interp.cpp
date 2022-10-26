#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <limits>

#include "doctest.h"
#include "interp.h"
#include "constants.h"

TEST_CASE("Testing lerp (float)") {

    CHECK(sharp::lerp((float)10.0, (float)20.0, (float)0) == (float)10.0);
    CHECK(sharp::lerp((float)10.0, (float)20.0, (float)1) == (float)20.0);
    CHECK(sharp::lerp((float)10.0, (float)20.0, (float)0.5) == (float)15.0);

    // make sure reordering the operations 
    // doesnt change the result!
    CHECK(sharp::lerp((float)20.0, (float)10.0, (float)0) == (float)20.0);
    CHECK(sharp::lerp((float)20.0, (float)10.0, (float)1) == (float)10.0);
    CHECK(sharp::lerp((float)20.0, (float)10.0, (float)0.5) == (float)15.0);

    const float inf = std::numeric_limits<float>::infinity();
    // test the infinity bound - should return first arg
    CHECK(sharp::lerp((float)10.0, (float)10.0, inf) == (float)10.0);

    // this one should return infinity
    CHECK(std::isinf(sharp::lerp((float)10.0, (float)50.0, inf)));
}

TEST_CASE("Testing lerp (double)") {

    CHECK(sharp::lerp((double)10.0, (double)20.0, (double)0) == (double)10.0);
    CHECK(sharp::lerp((double)10.0, (double)20.0, (double)1) == (double)20.0);
    CHECK(sharp::lerp((double)10.0, (double)20.0, (double)0.5) == (double)15.0);

    // make sure reordering the operations 
    // doesnt change the result!
    CHECK(sharp::lerp((double)20.0, (double)10.0, (double)0) == (double)20.0);
    CHECK(sharp::lerp((double)20.0, (double)10.0, (double)1) == (double)10.0);
    CHECK(sharp::lerp((double)20.0, (double)10.0, (double)0.5) == (double)15.0);

    const double inf = std::numeric_limits<double>::infinity();
    // test the infinity bound - should return first arg
    CHECK(sharp::lerp((double)10.0, (double)10.0, inf) == (double)10.0);

    // this one should return infinity
    CHECK(std::isinf(sharp::lerp((double)10.0, (double)50.0, inf)));
}
