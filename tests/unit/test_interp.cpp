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

#ifndef NO_QC
constexpr float hght[4] = {0, 100, 200, 300};
constexpr float pres[4] = {100000, 90000, 80000, 70000};
constexpr float nan_mid[3] = {1, nanval, 9};

TEST_CASE("Testing interp with NaN and MISSING data") {
    CHECK(sharp::interp_height(150, hght, nan_mid, 3) == 7);
    CHECK(sharp::interp_height(50, hght, nan_mid, 3) == 3);
    CHECK(sharp::interp_pressure(85000, pres, nan_mid, 3) ==
          doctest::Approx(6.82651615));
    CHECK(sharp::interp_pressure(95000, pres, nan_mid, 3) ==
          doctest::Approx(2.83893538));

    constexpr float nan_top[2] = {1, nanval};
    constexpr float nan_below[3] = {nanval, nanval, 9};
    constexpr float nan_above[3] = {1, nanval, nanval};
    CHECK(sharp::interp_height(50, hght, nan_top, 2) == sharp::MISSING);
    CHECK(sharp::interp_height(50, hght, nan_below, 3) == sharp::MISSING);
    CHECK(sharp::interp_height(150, hght, nan_above, 3) == sharp::MISSING);
    CHECK(sharp::interp_pressure(95000, pres, nan_top, 2) == sharp::MISSING);
    CHECK(sharp::interp_pressure(95000, pres, nan_below, 3) == sharp::MISSING);
    CHECK(sharp::interp_pressure(85000, pres, nan_above, 3) == sharp::MISSING);

    constexpr float hght_td[3] = {2950, 3000, 3050};
    for (const float bad : {nanval, sharp::MISSING}) {
        CAPTURE(bad);
        const float top[2] = {bad, 280};
        CHECK(sharp::interp_height(100, hght, top, 2) == 280);
        CHECK(sharp::interp_pressure(90000, pres, top, 2) == 280);

        const float bottom[2] = {280, bad};
        CHECK(sharp::interp_height(0, hght, bottom, 2) == 280);
        CHECK(sharp::interp_pressure(100000, pres, bottom, 2) == 280);

        const float td[3] = {270.0f, 269.5f, bad};
        CHECK(sharp::interp_height(3000, hght_td, td, 3) == 269.5f);
        CHECK(sharp::interp_pressure(90000, pres, td, 3) == 269.5f);
    }
    constexpr float td_mis[3] = {270.0f, 269.5f, sharp::MISSING};
    CHECK(sharp::interp_height(2999.9f, hght_td, td_mis, 3) ==
          doctest::Approx(269.501));
}

TEST_CASE("Testing interp at exact levels with a complete bracket") {
    constexpr float data4[4] = {sharp::MISSING, 2.25f, 7.75f, nanval};
    CHECK(sharp::interp_height(100, hght, data4, 4) == 2.25f);
    CHECK(sharp::interp_height(150, hght, data4, 4) == 5.0f);
    CHECK(sharp::interp_pressure(90000, pres, data4, 4) == 2.25f);
    CHECK(sharp::interp_pressure(85000, pres, data4, 4) == 4.91907024f);
}

TEST_CASE("Testing find_first with NaN and MISSING data") {
    constexpr float data3[3] = {1, 5, 9};
    CHECK(sharp::find_first_height(nanval, hght, data3, 3) == sharp::MISSING);
    CHECK(sharp::find_first_pressure(nanval, pres, data3, 3) == sharp::MISSING);

    CHECK(sharp::find_first_height(5, hght, nan_mid, 3) == 100);
    CHECK(sharp::find_first_pressure(5, pres, nan_mid, 3) ==
          doctest::Approx(89442.7));

    for (const float bad : {nanval, sharp::MISSING}) {
        CAPTURE(bad);
        const float v5_bad[2] = {5, bad};
        const float bad_v5[2] = {bad, 5};
        const float bad_v5_bad[3] = {bad, 5, bad};
        CHECK(sharp::find_first_height(5, hght, v5_bad, 2) == 0);
        CHECK(sharp::find_first_height(5, hght, bad_v5, 2) == 100);
        CHECK(sharp::find_first_height(5, hght, bad_v5_bad, 3) == 100);
        CHECK(sharp::find_first_pressure(5, pres, v5_bad, 2) == 100000);
        CHECK(sharp::find_first_pressure(5, pres, bad_v5, 2) == 90000);
        CHECK(sharp::find_first_pressure(5, pres, bad_v5_bad, 3) == 90000);

        CHECK(sharp::find_first_height(6, hght, v5_bad, 2) == sharp::MISSING);
        CHECK(sharp::find_first_pressure(6, pres, v5_bad, 2) == sharp::MISSING);
    }
}
#endif

TEST_CASE("Testing interp on an empty profile") {
    CHECK(sharp::interp_height(0, nullptr, nullptr, 0) == sharp::MISSING);
    CHECK(sharp::interp_pressure(100000, nullptr, nullptr, 0) ==
          sharp::MISSING);
    CHECK(sharp::interp_height(0, nullptr, nullptr, -1) == sharp::MISSING);
    CHECK(sharp::interp_pressure(100000, nullptr, nullptr, -1) ==
          sharp::MISSING);
}

TEST_CASE("Testing interp on a single-level profile") {
    constexpr float hght1[1] = {100};
    constexpr float pres1[1] = {85000};
    constexpr float data1[1] = {280.5f};

    CHECK(sharp::interp_height(100, hght1, data1, 1) == 280.5f);
    CHECK(sharp::interp_pressure(85000, pres1, data1, 1) == 280.5f);

    CHECK(sharp::interp_height(99, hght1, data1, 1) == sharp::MISSING);
    CHECK(sharp::interp_height(101, hght1, data1, 1) == sharp::MISSING);
    CHECK(sharp::interp_pressure(85001, pres1, data1, 1) == sharp::MISSING);
    CHECK(sharp::interp_pressure(84999, pres1, data1, 1) == sharp::MISSING);

#ifndef NO_QC
    for (const float bad : {sharp::MISSING, nanval}) {
        CAPTURE(bad);
        const float data[1] = {bad};
        CHECK(sharp::interp_height(100, hght1, data, 1) == sharp::MISSING);
        CHECK(sharp::interp_pressure(85000, pres1, data, 1) == sharp::MISSING);
    }
#endif
}

TEST_CASE("Testing find_first on a single-level profile") {
    constexpr float hght1[1] = {100};
    constexpr float pres1[1] = {85000};
    constexpr float v5[1] = {5};
    CHECK(sharp::find_first_height(5, hght1, v5, 1) == 100);
    CHECK(sharp::find_first_height(6, hght1, v5, 1) == sharp::MISSING);
    CHECK(sharp::find_first_pressure(5, pres1, v5, 1) == 85000);
    CHECK(sharp::find_first_pressure(6, pres1, v5, 1) == sharp::MISSING);
}
