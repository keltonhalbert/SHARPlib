#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <algorithm>
#include <functional>
#include <iostream>

#include "doctest.h"
#include <SHARPlib/algorithms.h>
#include <SHARPlib/layer.h>



TEST_CASE("lower bound") {

    float array[10] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    float val = 13;
    int N = 10;

    constexpr auto comp = std::less<float>();
    int sharp_lower = sharp::lower_bound(array, N, val, comp);
    auto std_lower = std::lower_bound(&array[0], &array[N-1], val); 

    int std_lower_idx = std_lower - &array[0];

    CHECK(sharp_lower == std_lower_idx);  
}

TEST_CASE("upper bound") {

    float array[10] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    float val = 13;
    int N = 10;

    constexpr auto comp = std::less<float>();
    int sharp_upper = sharp::upper_bound(array, N, val, comp);
    auto std_upper = std::upper_bound(&array[0], &array[N-1], val); 

    int std_upper_idx = std_upper - &array[0];

    CHECK(sharp_upper == std_upper_idx);  

}
