#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <SHARPlib/algorithms.h>
#include <SHARPlib/layer.h>

#include <algorithm>
#include <functional>

#include "doctest.h"

TEST_CASE("lower bound (in bounds)") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float array[N] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    constexpr float val = 13;

    constexpr auto comp = std::less<float>();
    const std::ptrdiff_t sharp_lower = sharp::lower_bound(array, N, val, comp);
    const auto std_lower = std::lower_bound(&array[0], &array[N - 1], val);
    const std::ptrdiff_t std_lower_idx = std_lower - &array[0];

    CHECK(sharp_lower == std_lower_idx);
    CHECK(array[sharp_lower] == array[std_lower_idx]);
}

TEST_CASE("lower bound (out of bounds)") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float array[N] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    float val = 7;

    constexpr auto comp = std::less<float>();
    std::ptrdiff_t sharp_lower = sharp::lower_bound(array, N, val, comp);
    auto std_lower = std::lower_bound(&array[0], &array[N - 1], val);
    std::ptrdiff_t std_lower_idx = std_lower - &array[0];

    CHECK(sharp_lower == 0);
    CHECK(sharp_lower == std_lower_idx);
    CHECK(array[sharp_lower] == array[std_lower_idx]);

    val = 27;
    sharp_lower = sharp::lower_bound(array, N, val, comp);
    std_lower = std::lower_bound(&array[0], &array[N - 1], val);
    std_lower_idx = std_lower - &array[0];

    CHECK(sharp_lower == N - 1);
    CHECK(sharp_lower == std_lower_idx);
    CHECK(array[sharp_lower] == array[std_lower_idx]);
}

TEST_CASE("upper bound (in bounds)") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float array[N] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    constexpr float val = 13;

    constexpr auto comp = std::less<float>();
    const std::ptrdiff_t sharp_upper = sharp::upper_bound(array, N, val, comp);
    const auto std_upper = std::upper_bound(&array[0], &array[N - 1], val);
    const std::ptrdiff_t std_upper_idx = std_upper - &array[0];

    CHECK(sharp_upper == std_upper_idx);
    CHECK(array[sharp_upper] == array[std_upper_idx]);
}

TEST_CASE("upper bound (out of bounds)") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float array[N] = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    float val = 7;

    constexpr auto comp = std::less<float>();
    std::ptrdiff_t sharp_upper = sharp::upper_bound(array, N, val, comp);
    auto std_upper = std::upper_bound(&array[0], &array[N - 1], val);
    std::ptrdiff_t std_upper_idx = std_upper - &array[0];

    CHECK(sharp_upper == 0);
    CHECK(sharp_upper == std_upper_idx);
    CHECK(array[sharp_upper] == array[std_upper_idx]);

    val = 27;
    sharp_upper = sharp::upper_bound(array, N, val, comp);
    std_upper = std::upper_bound(&array[0], &array[N - 1], val);
    std_upper_idx = std_upper - &array[0];

    CHECK(sharp_upper == N - 1);
    CHECK(sharp_upper == std_upper_idx);
    CHECK(array[sharp_upper] == array[std_upper_idx]);
}

TEST_CASE("Reverse ordered (upper|lower)bound") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float array[N] = {19, 18, 17, 16, 15, 14, 13, 12, 11, 10};
    constexpr float val = 16;

    constexpr auto comp = std::greater<float>();
    const std::ptrdiff_t sharp_lower = sharp::lower_bound(array, N, val, comp);
    const std::ptrdiff_t sharp_upper = sharp::upper_bound(array, N, val, comp);

    CHECK(sharp_lower == 3);
    CHECK(sharp_upper == 4);
}
