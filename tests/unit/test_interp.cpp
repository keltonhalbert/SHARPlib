#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>

#include <cmath>
#include <limits>

#include "doctest.h"

constexpr float infval = std::numeric_limits<float>::infinity();
constexpr float nanval = std::numeric_limits<float>::quiet_NaN();
TEST_CASE("Testing lerp (float)") {
    CHECK(sharp::lerp(10.0f, 20.0f, 0.f) == 10.0f);
    CHECK(sharp::lerp(10.0f, 20.0f, 1.f) == 20.0f);
    CHECK(sharp::lerp(10.0f, 20.0f, 0.5f) == 15.0f);

    // make sure reordering the operations
    // doesnt change the result!
    CHECK(sharp::lerp(20.0f, 10.0f, 0.f) == 20.0f);
    CHECK(sharp::lerp(20.0f, 10.0f, 1.f) == 10.0f);
    CHECK(sharp::lerp(20.0f, 10.0f, 0.5f) == 15.0f);

    constexpr float inf = std::numeric_limits<float>::infinity();
    // test the infinity bound - should return first arg
    CHECK(sharp::lerp(10.0f, 10.0f, inf) == 10.0f);

    // this one should return infinity
    CHECK(std::isinf(sharp::lerp(10.0f, 50.0f, inf)));
}

TEST_CASE("Testing interp_height") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float height_arr[N] = {100, 200, 300, 400, 500,
                                     600, 700, 800, 900, 1000};
    constexpr float data_arr[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

#ifndef NO_QC
    // test out of bounds values
    CHECK(sharp::interp_height(0, height_arr, data_arr, N) == sharp::MISSING);
    CHECK(sharp::interp_height(1100, height_arr, data_arr, N) ==
          sharp::MISSING);

    // test missing values
    CHECK(sharp::interp_height(sharp::MISSING, height_arr, data_arr, N) ==
          sharp::MISSING);
    CHECK(sharp::interp_height(sharp::MISSING, height_arr, data_arr, N) ==
          sharp::MISSING);
    CHECK(sharp::interp_height(infval, height_arr, data_arr, N) ==
          sharp::MISSING);
    CHECK(sharp::interp_height(nanval, height_arr, data_arr, N) ==
          sharp::MISSING);
#endif

    // test exact values along the edges of the arrays
    CHECK(sharp::interp_height(100, height_arr, data_arr, N) == 1);
    CHECK(sharp::interp_height(1000, height_arr, data_arr, N) == 10);

    // test an exact value in the middle
    CHECK(sharp::interp_height(500, height_arr, data_arr, N) == 5);

    // test between levels
    CHECK(sharp::interp_height(550, height_arr, data_arr, N) == 5.5);
    CHECK(sharp::interp_height(110, height_arr, data_arr, N) ==
          doctest::Approx(1.1));
    CHECK(sharp::interp_height(391, height_arr, data_arr, N) ==
          doctest::Approx(3.91));
}

TEST_CASE("Testing interp_pressure") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    constexpr float pres_arr[N] = {100000, 90000, 80000, 70000, 60000,
                                   50000,  40000, 30000, 20000, 10000};
    constexpr float data_arr[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

#ifndef NO_QC
    // test out of bounds values
    CHECK(sharp::interp_pressure(0, pres_arr, data_arr, N) == sharp::MISSING);
    CHECK(sharp::interp_pressure(110000, pres_arr, data_arr, N) ==
          sharp::MISSING);
    // test missing values
    CHECK(sharp::interp_pressure(sharp::MISSING, pres_arr, data_arr, N) ==
          sharp::MISSING);
    CHECK(sharp::interp_pressure(sharp::MISSING, pres_arr, data_arr, N) ==
          sharp::MISSING);
    CHECK(sharp::interp_pressure(infval, pres_arr, data_arr, N) ==
          sharp::MISSING);
    CHECK(sharp::interp_pressure(nanval, pres_arr, data_arr, N) ==
          sharp::MISSING);
#endif

    // test exact values along the edges of the arrays
    CHECK(sharp::interp_pressure(100000, pres_arr, data_arr, N) == 1);
    CHECK(sharp::interp_pressure(10000, pres_arr, data_arr, N) == 10);

    // test an exact value in the middle
    CHECK(sharp::interp_pressure(50000, pres_arr, data_arr, N) == 6);

    // test between levels -- float values were generated using
    // numpy.interp in python to get a known value
    CHECK(sharp::interp_pressure(97500, pres_arr, data_arr, N) ==
          doctest::Approx(1.2402969255248580));
    CHECK(sharp::interp_pressure(95000, pres_arr, data_arr, N) ==
          doctest::Approx(1.4868360226532418));
    CHECK(sharp::interp_pressure(92500, pres_arr, data_arr, N) ==
          doctest::Approx(1.7399502648876880));
}

TEST_CASE("Testing find_first_pressure") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    constexpr float pres_arr[N] = {100000, 90000, 80000, 70000, 60000,
                                   50000,  40000, 30000, 20000, 10000};
    constexpr float data_arr1[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    constexpr float data_arr2[N] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

    CHECK(sharp::find_first_pressure(5.0f, pres_arr, data_arr1, N) == 60000.0f);
    CHECK(sharp::find_first_pressure(5.5f, pres_arr, data_arr1, N) ==
          doctest::Approx(54772.3f));
    CHECK(sharp::find_first_pressure(5.0f, pres_arr, data_arr2, N) == 50000.0f);
    CHECK(sharp::find_first_pressure(5.5f, pres_arr, data_arr2, N) ==
          doctest::Approx(54772.3f));
}

TEST_CASE("Testing find_first_height") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float hght_arr[N] = {100, 200, 300, 400, 500,
                                   600, 700, 800, 900, 1000};
    constexpr float data_arr1[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    constexpr float data_arr2[N] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

    CHECK(sharp::find_first_height(5, hght_arr, data_arr1, N) == 500.0f);
    CHECK(sharp::find_first_height(5.5f, hght_arr, data_arr1, N) == 550.0f);
    CHECK(sharp::find_first_height(5, hght_arr, data_arr2, N) == 600.0f);
    CHECK(sharp::find_first_height(5.5f, hght_arr, data_arr2, N) == 550.0f);
}

#ifndef NO_QC
TEST_CASE("Testing find_first_height_with_missing") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float hght_arr[N] = {100, 200, 300, 400, 500,
                                   600, 700, 800, 900, 1000};
    constexpr float data_arr1[N] = {1, 2, sharp::MISSING, 4, 5, 6, 7, 8, 9, 10};
    constexpr float data_arr2[N] = {10, 9, sharp::MISSING, 7, 6, 5, 4, 3, 2, 1};

    CHECK(sharp::find_first_height(5, hght_arr, data_arr1, N) == 500.0f);
    CHECK(sharp::find_first_height(5.5f, hght_arr, data_arr1, N) == 550.0f);
    CHECK(sharp::find_first_height(5, hght_arr, data_arr2, N) == 600.0f);
    CHECK(sharp::find_first_height(5.5f, hght_arr, data_arr2, N) == 550.0f);

    CHECK(sharp::find_first_height(sharp::MISSING, hght_arr, data_arr1, N) ==
          sharp::MISSING);
    CHECK(sharp::find_first_height(3, hght_arr, data_arr1, N) == 300.0);
    CHECK(sharp::find_first_height(8, hght_arr, data_arr2, N) == 300.0);
}
#endif

#ifndef NO_QC
TEST_CASE("Testing find_first_pressure_with_missing") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    constexpr float pres_arr[N] = {100000, 90000, 80000, 70000, 60000,
                                   50000,  40000, 30000, 20000, 10000};
    constexpr float data_arr1[N] = {1, 2, sharp::MISSING, 4, 5, 6, 7, 8, 9, 10};
    constexpr float data_arr2[N] = {10, 9, sharp::MISSING, 7, 6, 5, 4, 3, 2, 1};

    CHECK(sharp::find_first_pressure(5.0f, pres_arr, data_arr1, N) == 60000.0f);
    CHECK(sharp::find_first_pressure(5.5f, pres_arr, data_arr1, N) ==
          doctest::Approx(54772.3f));
    CHECK(sharp::find_first_pressure(5.0f, pres_arr, data_arr2, N) == 50000.0f);
    CHECK(sharp::find_first_pressure(5.5f, pres_arr, data_arr2, N) ==
          doctest::Approx(54772.3f));

    CHECK(sharp::find_first_pressure(sharp::MISSING, pres_arr, data_arr2, N) ==
          sharp::MISSING);
    CHECK(sharp::find_first_pressure(3.0f, pres_arr, data_arr1, N) ==
          doctest::Approx(79372.6f));
    CHECK(sharp::find_first_pressure(8.0f, pres_arr, data_arr2, N) ==
          doctest::Approx(79372.6f));
}
#endif
