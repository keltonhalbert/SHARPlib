#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>

#include <cmath>
#include <limits>

#include "doctest.h"

TEST_CASE("Testing lerp (float)") {
    CHECK(sharp::lerp(10.0f, 20.0f, 0.f) == 10.0f);
    CHECK(sharp::lerp(10.0f, 20.0f, 1.f) == 20.0f);
    CHECK(sharp::lerp(10.0f, 20.0f, 0.5f) == 15.0f);

    // make sure reordering the operations
    // doesnt change the result!
    CHECK(sharp::lerp(20.0f, 10.0f, 0.f) == 20.0f);
    CHECK(sharp::lerp(20.0f, 10.0f, 1.f) == 10.0f);
    CHECK(sharp::lerp(20.0f, 10.0f, 0.5f) == 15.0f);

    const float inf = std::numeric_limits<float>::infinity();
    // test the infinity bound - should return first arg
    CHECK(sharp::lerp(10.0f, 10.0f, inf) == 10.0f);

    // this one should return infinity
    CHECK(std::isinf(sharp::lerp(10.0f, 50.0f, inf)));
}

TEST_CASE("Testing interp_height") {
    float height_arr[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    float data_arr[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int arr_len = 10;

#ifndef NO_QC
    // test out of bounds values
    CHECK(sharp::interp_height(0, height_arr, data_arr, arr_len) ==
          sharp::MISSING);
    CHECK(sharp::interp_height(1100, height_arr, data_arr, arr_len) ==
          sharp::MISSING);
#endif

    // test exact values along the edges of the arrays
    CHECK(sharp::interp_height(100, height_arr, data_arr, arr_len) == 1);
    CHECK(sharp::interp_height(1000, height_arr, data_arr, arr_len) == 10);

    // test an exact value in the middle
    CHECK(sharp::interp_height(500, height_arr, data_arr, arr_len) == 5);

    // test between levels
    CHECK(sharp::interp_height(550, height_arr, data_arr, arr_len) == 5.5);
    CHECK(sharp::interp_height(110, height_arr, data_arr, arr_len) ==
          doctest::Approx(1.1));
    CHECK(sharp::interp_height(391, height_arr, data_arr, arr_len) ==
          doctest::Approx(3.91));

#ifndef NO_QC
    // test missing values
    constexpr float inf = std::numeric_limits<float>::infinity();
    // constexpr float nan = std::numeric_limits<float>::quiet_NaN();
    CHECK(sharp::interp_height(sharp::MISSING, height_arr, data_arr, arr_len) ==
          sharp::MISSING);
    CHECK(sharp::interp_height(sharp::MISSING, height_arr, data_arr, arr_len) ==
          sharp::MISSING);
    // To-Do: NaN check doesn't work, so fix it
    CHECK(sharp::interp_height(inf, height_arr, data_arr, arr_len) ==
          sharp::MISSING);
    // CHECK(sharp::interp_height(nan, height_arr, data_arr, arr_len) ==
    // sharp::MISSING);
#endif
}

TEST_CASE("Testing interp_pressure") {
    float pres_arr[10] = {1000, 900, 800, 700, 600, 500, 400, 300, 200, 100};
    float data_arr[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int arr_len = 10;

#ifndef NO_QC
    // test out of bounds values
    CHECK(sharp::interp_pressure(0, pres_arr, data_arr, arr_len) ==
          sharp::MISSING);
    CHECK(sharp::interp_pressure(1100, pres_arr, data_arr, arr_len) ==
          sharp::MISSING);
#endif

    // test exact values along the edges of the arrays
    CHECK(sharp::interp_pressure(1000, pres_arr, data_arr, arr_len) == 1);
    CHECK(sharp::interp_pressure(100, pres_arr, data_arr, arr_len) == 10);

    // test an exact value in the middle
    CHECK(sharp::interp_pressure(500, pres_arr, data_arr, arr_len) == 6);

    // test between levels -- float values were generated using
    // numpy.interp in python to get a known value
    CHECK(sharp::interp_pressure(975, pres_arr, data_arr, arr_len) ==
          doctest::Approx(1.2402969255248580));
    CHECK(sharp::interp_pressure(950, pres_arr, data_arr, arr_len) ==
          doctest::Approx(1.4868360226532418));
    CHECK(sharp::interp_pressure(925, pres_arr, data_arr, arr_len) ==
          doctest::Approx(1.7399502648876880));

#ifndef NO_QC
    // test missing values
    constexpr float inf = std::numeric_limits<float>::infinity();
    // constexpr float nan = std::numeric_limits<float>::quiet_NaN();
    CHECK(sharp::interp_pressure(sharp::MISSING, pres_arr, data_arr, arr_len) ==
          sharp::MISSING);
    CHECK(sharp::interp_pressure(sharp::MISSING, pres_arr, data_arr, arr_len) ==
          sharp::MISSING);
    CHECK(sharp::interp_pressure(inf, pres_arr, data_arr, arr_len) ==
          sharp::MISSING);
#endif
    // To-Do: NaN check doesn't work, so fix it
    // CHECK(sharp::interp_pressure(nan, pres_arr, data_arr, arr_len) ==
    // sharp::MISSING);
}
