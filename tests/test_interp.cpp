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

TEST_CASE("Testing interp_height") {

    float height_arr[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    float data_arr[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int arr_len = 10;

    // test out of bounds values
    CHECK(sharp::interp_height(0, height_arr, data_arr, arr_len) 
                                              == sharp::MISSING);
    CHECK(sharp::interp_height(1100, height_arr, data_arr, arr_len) 
                                                 == sharp::MISSING);

    // test exact values along the edges of the arrays
    CHECK(sharp::interp_height(100, height_arr, data_arr, arr_len) == 1);
    CHECK(sharp::interp_height(1000, height_arr, data_arr, arr_len) == 10);

    // test an exact value in the middle
    CHECK(sharp::interp_height(500, height_arr, data_arr, arr_len) == 5);

    // test between levels
    CHECK(sharp::interp_height(550, height_arr, data_arr, arr_len) == 5.5);
    CHECK(sharp::interp_height(110, height_arr, data_arr, arr_len) 
                                          == doctest::Approx(1.1));
    CHECK(sharp::interp_height(391, height_arr, data_arr, arr_len) 
                                         == doctest::Approx(3.91));

}



