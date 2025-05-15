#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/layer.h>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "doctest.h"
constexpr float nanval = std::numeric_limits<float>::quiet_NaN();
constexpr float infval = std::numeric_limits<float>::infinity();

TEST_CASE("Testing HeightLayer structs") {
    CHECK_THROWS_AS(sharp::HeightLayer layer1(100000.0f, 0.0f),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::HeightLayer layer2(nanval, nanval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::HeightLayer layer3(100000.0f, nanval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::HeightLayer layer4(nanval, 50000.0f),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::HeightLayer layer5(infval, infval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::HeightLayer layer6(100000.0f, infval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::HeightLayer layer7(infval, 50000.0f),
                    const std::range_error&);
}

TEST_CASE("Testing PressureLayer structs") {
    CHECK_THROWS_AS(sharp::PressureLayer layer1(10000.0, 100000.0),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::PressureLayer layer2(nanval, nanval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::PressureLayer layer3(100000.0f, nanval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::PressureLayer layer4(nanval, 50000.0f),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::PressureLayer layer5(infval, infval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::PressureLayer layer6(100000.0f, infval),
                    const std::range_error&);
    CHECK_THROWS_AS(sharp::PressureLayer layer7(infval, 50000.0f),
                    const std::range_error&);
}

TEST_CASE("Testing conversion between layers") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    constexpr float pres[N] = {100000.0, 90000.0, 80000.0, 70000.0, 60000.0,
                               50000.0,  40000.0, 30000.0, 20000.0, 10000.0};
    constexpr float hght[N] = {0.0,    500.0,  1500.0, 2500.0,  4000.0,
                               5500.0, 7500.0, 8500.0, 10500.0, 12500.0};

    // check within bounds
    sharp::HeightLayer hlyr = {0.0, 3000.0};
    sharp::PressureLayer plyr = {100000.0, 75000.0};

    auto out_plyr = sharp::height_layer_to_pressure(hlyr, pres, hght, N);
    auto out_hlyr = sharp::pressure_layer_to_height(plyr, pres, hght, N);

    CHECK(out_plyr.bottom == 100000.0);
    CHECK(out_hlyr.bottom == 0.0);
    CHECK(out_plyr.top == doctest::Approx(66666.7));
    CHECK(out_hlyr.top == doctest::Approx(1983.32));

    // check ot of bounds
    sharp::HeightLayer h_oob1 = {-100, 250.0};
    sharp::HeightLayer h_oob2 = {11500.0, 14000.0};
    sharp::PressureLayer p_oob1 = {115000.0, 90000.0};
    sharp::PressureLayer p_oob2 = {15000.0, 5000.0};

    auto oob1 = sharp::height_layer_to_pressure(h_oob1, pres, hght, N);
    auto oob2 = sharp::height_layer_to_pressure(h_oob2, pres, hght, N);
    auto oob3 = sharp::pressure_layer_to_height(p_oob1, pres, hght, N);
    auto oob4 = sharp::pressure_layer_to_height(p_oob2, pres, hght, N);

    CHECK(oob1.bottom == sharp::MISSING);
    CHECK(oob2.bottom == sharp::MISSING);
    CHECK(oob3.bottom == sharp::MISSING);
    CHECK(oob4.bottom == sharp::MISSING);
    CHECK(oob1.top == sharp::MISSING);
    CHECK(oob2.top == sharp::MISSING);
    CHECK(oob3.top == sharp::MISSING);
    CHECK(oob4.top == sharp::MISSING);
}

TEST_CASE("Testing layer bounds checking and searching") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float pres[N] = {100000.0, 90000.0, 80000.0, 70000.0, 60000.0,
                               50000.0,  40000.0, 30000.0, 20000.0, 10000.0};

    sharp::PressureLayer out_of_bounds_1 = {100000.0, 5000.0};
    sharp::PressureLayer out_of_bounds_2 = {110000.0, 10000.0};
    sharp::PressureLayer out_of_bounds_3 = {110000.0, 5000.0};
    sharp::PressureLayer out_of_bounds_4 = {5000.0, 2500.0};
    sharp::PressureLayer out_of_bounds_5 = {115000.0, 110000.0};

    const sharp::LayerIndex idx_1 =
        sharp::get_layer_index(out_of_bounds_1, pres, N);
    const sharp::LayerIndex idx_2 =
        sharp::get_layer_index(out_of_bounds_2, pres, N);
    const sharp::LayerIndex idx_3 =
        sharp::get_layer_index(out_of_bounds_3, pres, N);
    const sharp::LayerIndex idx_4 =
        sharp::get_layer_index(out_of_bounds_4, pres, N);
    const sharp::LayerIndex idx_5 =
        sharp::get_layer_index(out_of_bounds_5, pres, N);

    CHECK(idx_1.kbot == 0);
    CHECK(idx_1.ktop == 9);
    CHECK(idx_2.kbot == 0);
    CHECK(idx_2.ktop == 9);
    CHECK(idx_3.kbot == 0);
    CHECK(idx_3.ktop == 9);
    CHECK(idx_4.kbot == 9);
    CHECK(idx_4.ktop == 9);
    CHECK(idx_5.kbot == 0);
    CHECK(idx_5.ktop == 0);
}

TEST_CASE("Testing layer_max over pressure layer") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    constexpr float pres[N] = {100000.0, 90000.0, 80000.0, 70000.0, 60000.0,
                               50000.0,  40000.0, 30000.0, 20000.0, 10000.0};
    constexpr float data[N] = {0, 0, 0, 0, 10, 0, 0, 0, 0, 0};
    sharp::PressureLayer layer1(100000.0, 60000.0);  // max at top of layer
    sharp::PressureLayer layer2(110000.0, 5000.0);   // out of bounds recovery
    sharp::PressureLayer layer3(100000.0, 85000.0);  // max is not in layer
    sharp::PressureLayer layer4(80000.0, 65000.0);   // max is just above layer
    sharp::PressureLayer layer5(55000.0, 40000.0);   // max is just below layer
    float pmax = -9999.0;

    CHECK(sharp::layer_max(layer1, pres, data, N, &pmax) == 10.0);
    CHECK(pmax == 60000.0);
    CHECK(sharp::layer_max(layer2, pres, data, N, nullptr) == 10.0);
    CHECK(sharp::layer_max(layer3, pres, data, N, nullptr) == 0.0);
    CHECK(sharp::layer_max(layer4, pres, data, N, &pmax) ==
          doctest::Approx(4.80749));
    CHECK(pmax == 65000.0);
    CHECK(sharp::layer_max(layer5, pres, data, N, &pmax) ==
          doctest::Approx(5.22758));
}

TEST_CASE("Testing layer_min over pressure layer") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    float pres[N] = {100000.0, 90000.0, 80000.0, 70000.0, 60000.0,
                     50000.0,  40000.0, 30000.0, 20000.0, 10000.0};
    float data[N] = {0, 0, 0, 0, -10, 0, 0, 0, 0, 0};
    sharp::PressureLayer layer1(100000.0, 60000.0);  // min at top of layer
    sharp::PressureLayer layer2(110000.0, 5000.0);   // out of bounds recovery
    sharp::PressureLayer layer3(100000.0, 85000.0);  // min is not in layer
    sharp::PressureLayer layer4(80000.0, 65000.0);   // min is just above layer
    sharp::PressureLayer layer5(55000.0, 40000.0);   // min is just below layer

    CHECK(sharp::layer_min(layer1, pres, data, N) == -10.0);
    CHECK(sharp::layer_min(layer2, pres, data, N) == -10.0);
    CHECK(sharp::layer_min(layer3, pres, data, N) == 0.0);
    CHECK(sharp::layer_min(layer4, pres, data, N) == doctest::Approx(-4.80749));
    CHECK(sharp::layer_min(layer5, pres, data, N) == doctest::Approx(-5.22758));
}

TEST_CASE("Testing layer_max over height layer") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float hght[N] = {1000.0, 2000.0, 3000.0, 4000.0, 5000.0,
                               6000.0, 7000.0, 8000.0, 9000.0, 10000.0};
    constexpr float data[N] = {0, 0, 0, 0, 10, 0, 0, 0, 0, 0};
    sharp::HeightLayer layer1(1000.0, 6000.0);   // max at top of layer
    sharp::HeightLayer layer2(-100.0, 20000.0);  // out of bounds recovery
    sharp::HeightLayer layer3(6000.0, 10000.0);  // max is not in layer
    sharp::HeightLayer layer4(1000.0, 4750.0);   // max is just above layer
    sharp::HeightLayer layer5(5250.0, 9000.0);   // max is just below layer

    CHECK(sharp::layer_max(layer1, hght, data, N, nullptr) == 10.0);
    CHECK(sharp::layer_max(layer2, hght, data, N, nullptr) == 10.0);
    CHECK(sharp::layer_max(layer3, hght, data, N, nullptr) == 0.0);
    CHECK(sharp::layer_max(layer4, hght, data, N, nullptr) == 7.5);
    CHECK(sharp::layer_max(layer5, hght, data, N, nullptr) == 7.5);
}

TEST_CASE("Testing layer_min over height layer") {
    constexpr std::ptrdiff_t N = 10;
    constexpr float hght[N] = {1000.0, 2000.0, 3000.0, 4000.0, 5000.0,
                               6000.0, 7000.0, 8000.0, 9000.0, 10000.0};
    constexpr float data[N] = {0, 0, 0, 0, -10, 0, 0, 0, 0, 0};
    sharp::HeightLayer layer1(1000.0, 6000.0);   // min at top of layer
    sharp::HeightLayer layer2(-100.0, 20000.0);  // out of bounds recovery
    sharp::HeightLayer layer3(6000.0, 10000.0);  // min is not in layer
    sharp::HeightLayer layer4(1000.0, 4750.0);   // min is just above layer
    sharp::HeightLayer layer5(5250.0, 9000.0);   // min is just below layer

    CHECK(sharp::layer_min(layer1, hght, data, 10) == -10.0);
    CHECK(sharp::layer_min(layer2, hght, data, 10) == -10.0);
    CHECK(sharp::layer_min(layer3, hght, data, 10) == 0.0);
    CHECK(sharp::layer_min(layer4, hght, data, 10) == -7.5);
    CHECK(sharp::layer_min(layer5, hght, data, 10) == -7.5);
}

TEST_CASE("Testing layer_mean over a pressure layer") {
    constexpr std::ptrdiff_t N = 10;
    // pressure is always in Pa
    float pres[N] = {100000.0, 90000.0, 80000.0, 70000.0, 60000.0,
                     50000.0,  40000.0, 30000.0, 20000.0, 10000.0};

    float data[N] = {0, 0, 0, 0, 10, 0, 0, 0, 0, 0};

    sharp::PressureLayer layer1(100000.0, 60000.0);  // max at top of layer
    sharp::PressureLayer layer2(110000.0, 5000.0);   // out of bounds recovery
    sharp::PressureLayer layer3(100000.0, 10000.0);  // full layer

    CHECK(sharp::layer_mean(layer1, pres, data, N) == 1.25);
    CHECK(sharp::layer_mean(layer2, pres, data, N) == doctest::Approx(1.1111));
    CHECK(sharp::layer_mean(layer3, pres, data, N) == doctest::Approx(1.1111));
}
