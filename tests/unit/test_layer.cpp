#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/layer.h>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <new>
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

// A layer_min or layer_max result: the value and the level reported with it.
struct Extreme {
    float value;
    float level;
};

// Checks layer_min and layer_max over one layer.
template <typename L>
static void check_min_max(const L layer, const float coord[],
                          const float data[], const std::ptrdiff_t N,
                          const Extreme min, const Extreme max) {
    INFO("layer ", layer.bottom, " to ", layer.top);
    float min_lvl = 0.0f;
    float max_lvl = 0.0f;
    const float min_val = sharp::layer_min(layer, coord, data, N, &min_lvl);
    const float max_val = sharp::layer_max(layer, coord, data, N, &max_lvl);
    CHECK(min_val == doctest::Approx(min.value));
    CHECK(min_lvl == min.level);
    CHECK(max_val == doctest::Approx(max.value));
    CHECK(max_lvl == max.level);
}

// QC builds skip MISSING and NaN data in layer_min and layer_max, and return
// MISSING for a layer wholly outside the profile (SHARPlib-5yf.2). Each
// "was" comment gives the min and max measured before P1, which is before
// this change and the interp change of SHARPlib-5yf.1. Where the interp
// change alone gave something else, that output follows in brackets.
#ifndef NO_QC
constexpr float MISSING = sharp::MISSING;

TEST_CASE("Testing layer_min and layer_max with missing data at a boundary") {
    constexpr std::ptrdiff_t N = 6;
    constexpr float hght[N] = {0, 500, 1000, 1500, 2000, 2500};
    constexpr float pres[N] = {100000, 95000, 90000, 85000, 80000, 75000};
    constexpr float tmpk_nan[N] = {258, 258, nanval, 268, 268, 268};
    constexpr float tmpk_mis[N] = {258, 258, MISSING, 268, 268, 268};

    // NaN next to the bottom of the layer
    // was (NaN, 750), (NaN, 750) [(260.5, 750), (268, 1500)]
    check_min_max(sharp::HeightLayer(750, 2500), hght, tmpk_nan, N,
                  {260.5f, 750}, {268, 1500});
    // was (NaN, 92500), (NaN, 92500) [(260.398, 92500), (268, 85000)]
    check_min_max(sharp::PressureLayer(92500, 75000), pres, tmpk_nan, N,
                  {260.397675f, 92500}, {268, 85000});

    // NaN next to the top of the layer
    // was (258, 0), (268, 1250) [(258, 0), (265.5, 1250)]
    check_min_max(sharp::HeightLayer(0, 1250), hght, tmpk_nan, N, {258, 0},
                  {265.5f, 1250});
    // was (258, 0), (258, 0) [(258, 0), (260.5, 750)]
    check_min_max(sharp::HeightLayer(0, 750), hght, tmpk_nan, N, {258, 0},
                  {260.5f, 750});
    // was (258, 100000), (268, 87500) [(258, 100000), (265.394, 87500)]
    check_min_max(sharp::PressureLayer(100000, 87500), pres, tmpk_nan, N,
                  {258, 100000}, {265.393829f, 87500});

    // MISSING next to the bottom and top of the layer
    // was (MISSING, 1000), (268, 1500)
    check_min_max(sharp::HeightLayer(750, 2500), hght, tmpk_mis, N,
                  {260.5f, 750}, {268, 1500});
    // was (MISSING, 1000), (265.5, 1250)
    check_min_max(sharp::HeightLayer(0, 1250), hght, tmpk_mis, N, {258, 0},
                  {265.5f, 1250});
}

TEST_CASE("Testing layer_min and layer_max skip MISSING and NaN levels") {
    constexpr std::ptrdiff_t N = 5;
    constexpr float hght[N] = {0, 100, 200, 300, 400};
    constexpr float pres[N] = {100000, 90000, 80000, 70000, 60000};
    constexpr float data_mis[N] = {3, 1, MISSING, 6, 4};
    constexpr float data_nan[N] = {3, 1, nanval, 6, 4};

    // was (MISSING, 200), (6, 300)
    check_min_max(sharp::HeightLayer(0, 400), hght, data_mis, N, {1, 100},
                  {6, 300});
    // was (MISSING, 80000), (6, 70000)
    check_min_max(sharp::PressureLayer(100000, 60000), pres, data_mis, N,
                  {1, 90000}, {6, 70000});
    // was (1, 100), (6, 300)
    check_min_max(sharp::HeightLayer(0, 400), hght, data_nan, N, {1, 100},
                  {6, 300});
    // was (1, 90000), (6, 70000)
    check_min_max(sharp::PressureLayer(100000, 60000), pres, data_nan, N,
                  {1, 90000}, {6, 70000});

    // an isolated NaN away from the layer bottom and top
    constexpr float data_iso[N] = {1, 2, nanval, 4, 5};
    // was (1, 0), (5, 400)
    check_min_max(sharp::HeightLayer(0, 400), hght, data_iso, N, {1, 0},
                  {5, 400});
}

TEST_CASE("Testing layer_min and layer_max with a MISSING endpoint") {
    // The endpoint has no valid level on one side of it, so it interpolates
    // to MISSING and the result comes from the valid levels.
    constexpr std::ptrdiff_t N = 5;
    constexpr float hght[N] = {0, 100, 200, 300, 400};
    constexpr float pres[N] = {100000, 90000, 80000, 70000, 60000};
    constexpr float bot_mis[N] = {MISSING, MISSING, 3, 1, 6};
    constexpr float bot_nan[N] = {nanval, nanval, 3, 1, 6};
    constexpr float top_mis[N] = {3, 1, 6, MISSING, MISSING};
    constexpr float top_nan[N] = {3, 1, 6, nanval, nanval};

    // bottom endpoint
    // was (MISSING, 50), (6, 400)
    check_min_max(sharp::HeightLayer(50, 400), hght, bot_mis, N, {1, 300},
                  {6, 400});
    // was (NaN, 50), (NaN, 50) [(MISSING, 50), (6, 400)]
    check_min_max(sharp::HeightLayer(50, 400), hght, bot_nan, N, {1, 300},
                  {6, 400});
    // was (MISSING, 95000), (6, 60000)
    check_min_max(sharp::PressureLayer(95000, 60000), pres, bot_mis, N,
                  {1, 70000}, {6, 60000});
    // was (NaN, 95000), (NaN, 95000) [(MISSING, 95000), (6, 60000)]
    check_min_max(sharp::PressureLayer(95000, 60000), pres, bot_nan, N,
                  {1, 70000}, {6, 60000});

    // top endpoint
    // was (MISSING, 300), (6, 200)
    check_min_max(sharp::HeightLayer(0, 350), hght, top_mis, N, {1, 100},
                  {6, 200});
    // was (1, 100), (6, 200) [(MISSING, 350), (6, 200)]
    check_min_max(sharp::HeightLayer(0, 350), hght, top_nan, N, {1, 100},
                  {6, 200});
    // was (MISSING, 70000), (6, 80000)
    check_min_max(sharp::PressureLayer(100000, 65000), pres, top_mis, N,
                  {1, 90000}, {6, 80000});
    // was (1, 90000), (6, 80000) [(MISSING, 65000), (6, 80000)]
    check_min_max(sharp::PressureLayer(100000, 65000), pres, top_nan, N,
                  {1, 90000}, {6, 80000});
}

TEST_CASE("Testing layer_min and layer_max over a layer with no valid data") {
    constexpr std::ptrdiff_t N = 5;
    constexpr float hght[N] = {0, 100, 200, 300, 400};
    constexpr float pres[N] = {100000, 90000, 80000, 70000, 60000};
    constexpr float data_mis[N] = {3, 1, MISSING, MISSING, MISSING};
    constexpr float data_nan[N] = {3, 1, nanval, nanval, nanval};
    constexpr float all_mis[N] = {MISSING, MISSING, MISSING, MISSING, MISSING};

    // was (MISSING, 150), (MISSING, 150)
    check_min_max(sharp::HeightLayer(150, 400), hght, data_mis, N,
                  {MISSING, 150}, {MISSING, 150});
    // was (NaN, 150), (NaN, 150) [(MISSING, 150), (MISSING, 150)]
    check_min_max(sharp::HeightLayer(150, 400), hght, data_nan, N,
                  {MISSING, 150}, {MISSING, 150});
    // was (MISSING, 85000), (MISSING, 85000)
    check_min_max(sharp::PressureLayer(85000, 60000), pres, data_mis, N,
                  {MISSING, 85000}, {MISSING, 85000});
    // was (MISSING, 0), (MISSING, 0)
    check_min_max(sharp::HeightLayer(0, 400), hght, all_mis, N, {MISSING, 0},
                  {MISSING, 0});
}

TEST_CASE("Testing layer_min and layer_max over layers outside the profile") {
    // Complete data. layer_min is unchanged, and layer_max used to return a
    // value from outside the requested layer.
    constexpr std::ptrdiff_t N = 3;
    constexpr float hght[N] = {0, 500, 1000};
    constexpr float tmpk_hght[N] = {258, 268, 278};
    constexpr float pres[N] = {100000, 90000, 80000};
    constexpr float tmpk_pres[N] = {278, 268, 258};

    // was (MISSING, 1500), (278, 1000)
    check_min_max(sharp::HeightLayer(1500, 2000), hght, tmpk_hght, N,
                  {MISSING, 1500}, {MISSING, 1500});
    // was (MISSING, -100), (258, 0)
    check_min_max(sharp::HeightLayer(-500, -100), hght, tmpk_hght, N,
                  {MISSING, -100}, {MISSING, -100});
    // was (MISSING, 70000), (258, 80000)
    check_min_max(sharp::PressureLayer(70000, 60000), pres, tmpk_pres, N,
                  {MISSING, 70000}, {MISSING, 70000});
    // was (MISSING, 105000), (278, 100000)
    check_min_max(sharp::PressureLayer(110000, 105000), pres, tmpk_pres, N,
                  {MISSING, 105000}, {MISSING, 105000});

    // without a level pointer
    CHECK(sharp::layer_max(sharp::HeightLayer(1500, 2000), hght, tmpk_hght,
                           N) == MISSING);
    CHECK(sharp::layer_max(sharp::PressureLayer(110000, 105000), pres,
                           tmpk_pres, N) == MISSING);
}
#endif

TEST_CASE("Testing layer_min and layer_max over layers at the profile edge") {
    // Layers that touch the profile at one point or partly overlap it give
    // the same results as before, in both builds.
    constexpr std::ptrdiff_t N = 3;
    constexpr float hght[N] = {0, 500, 1000};
    constexpr float tmpk_hght[N] = {258, 268, 278};
    constexpr float pres[N] = {100000, 90000, 80000};
    constexpr float tmpk_pres[N] = {278, 268, 258};

    // touching at one point
    check_min_max(sharp::HeightLayer(1000, 2000), hght, tmpk_hght, N,
                  {278, 1000}, {278, 1000});
    check_min_max(sharp::HeightLayer(-500, 0), hght, tmpk_hght, N, {258, 0},
                  {258, 0});
    check_min_max(sharp::PressureLayer(80000, 70000), pres, tmpk_pres, N,
                  {258, 80000}, {258, 80000});
    check_min_max(sharp::PressureLayer(105000, 100000), pres, tmpk_pres, N,
                  {278, 100000}, {278, 100000});

    // partly above and partly below
    check_min_max(sharp::HeightLayer(500, 1500), hght, tmpk_hght, N, {268, 500},
                  {278, 1000});
    check_min_max(sharp::HeightLayer(-500, 500), hght, tmpk_hght, N, {258, 0},
                  {268, 500});
    check_min_max(sharp::PressureLayer(90000, 70000), pres, tmpk_pres, N,
                  {258, 80000}, {268, 90000});
    check_min_max(sharp::PressureLayer(110000, 90000), pres, tmpk_pres, N,
                  {268, 90000}, {278, 100000});
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

// layer_mean over layers that leave the profile (SHARPlib-4v5). A layer wholly
// outside the profile has no mean and returns MISSING in every build. Each
// "was" comment is the output measured before this change, in both QC and
// NO_QC builds. The height overload clipped one end of such a layer and not
// the other, which inverted it, and the conversion to pressure then threw
// std::range_error.
constexpr std::ptrdiff_t LM_N = 3;
constexpr float lm_pres[LM_N] = {100000, 95000, 90000};
constexpr float lm_data[LM_N] = {300, 297, 294};
// the same profile with the surface at 0 m and at 300 m
constexpr float lm_hght_0[LM_N] = {0, 500, 1000};
constexpr float lm_hght_300[LM_N] = {300, 800, 1300};

TEST_CASE("Testing layer_mean over layers outside the profile") {
    for (const float* hght : {lm_hght_0, lm_hght_300}) {
        CAPTURE(hght[0]);
        for (const bool agl : {false, true}) {
            CAPTURE(agl);
            // wholly above: was std::range_error
            CHECK(sharp::layer_mean(sharp::HeightLayer(1500, 2000), hght,
                                    lm_pres, lm_data, LM_N,
                                    agl) == sharp::MISSING);
            // wholly below: was std::range_error
            CHECK(sharp::layer_mean(sharp::HeightLayer(-500, -100), hght,
                                    lm_pres, lm_data, LM_N,
                                    agl) == sharp::MISSING);
        }
    }
    // below a 300 m surface in meters MSL: was std::range_error
    CHECK(sharp::layer_mean(sharp::HeightLayer(0, 200), lm_hght_300, lm_pres,
                            lm_data, LM_N, false) == sharp::MISSING);

    // pressure layers, unchanged
    CHECK(sharp::layer_mean(sharp::PressureLayer(85000, 80000), lm_pres,
                            lm_data, LM_N) == sharp::MISSING);
    CHECK(sharp::layer_mean(sharp::PressureLayer(110000, 105000), lm_pres,
                            lm_data, LM_N) == sharp::MISSING);

    // one level at 0 m
    constexpr float pres[1] = {100000};
    constexpr float hght[1] = {0};
    constexpr float data[1] = {300};
    // was std::range_error
    CHECK(sharp::layer_mean(sharp::HeightLayer(100, 200), hght, pres, data,
                            1) == sharp::MISSING);
    // was std::range_error
    CHECK(sharp::layer_mean(sharp::HeightLayer(-200, -100), hght, pres, data,
                            1) == sharp::MISSING);
}

TEST_CASE("Testing layer_mean over layers at the profile edge") {
    // Layers that touch the profile at one point have no depth and return
    // MISSING. Layers that partly overlap it are clipped to it. Both are
    // unchanged.
    for (const float* hght : {lm_hght_0, lm_hght_300}) {
        CAPTURE(hght[0]);
        // touching at one point, meters AGL
        CHECK(sharp::layer_mean(sharp::HeightLayer(1000, 2000), hght, lm_pres,
                                lm_data, LM_N, true) == sharp::MISSING);
        CHECK(sharp::layer_mean(sharp::HeightLayer(-500, 0), hght, lm_pres,
                                lm_data, LM_N, true) == sharp::MISSING);

        // partly above and partly below, meters AGL
        CHECK(sharp::layer_mean(sharp::HeightLayer(500, 1500), hght, lm_pres,
                                lm_data, LM_N,
                                true) == doctest::Approx(295.5f));
        CHECK(sharp::layer_mean(sharp::HeightLayer(-500, 500), hght, lm_pres,
                                lm_data, LM_N,
                                true) == doctest::Approx(298.5f));
    }

    // a 300 m surface in meters MSL
    CHECK(sharp::layer_mean(sharp::HeightLayer(1300, 2000), lm_hght_300,
                            lm_pres, lm_data, LM_N) == sharp::MISSING);
    CHECK(sharp::layer_mean(sharp::HeightLayer(-500, 300), lm_hght_300, lm_pres,
                            lm_data, LM_N) == sharp::MISSING);
    CHECK(sharp::layer_mean(sharp::HeightLayer(1000, 2000), lm_hght_300,
                            lm_pres, lm_data,
                            LM_N) == doctest::Approx(294.909698f));
    CHECK(sharp::layer_mean(sharp::HeightLayer(-500, 500), lm_hght_300, lm_pres,
                            lm_data, LM_N) == doctest::Approx(299.40921f));

    // pressure layers
    CHECK(sharp::layer_mean(sharp::PressureLayer(90000, 80000), lm_pres,
                            lm_data, LM_N) == sharp::MISSING);
    CHECK(sharp::layer_mean(sharp::PressureLayer(105000, 100000), lm_pres,
                            lm_data, LM_N) == sharp::MISSING);
    CHECK(sharp::layer_mean(sharp::PressureLayer(95000, 80000), lm_pres,
                            lm_data, LM_N) == doctest::Approx(295.5f));
    CHECK(sharp::layer_mean(sharp::PressureLayer(105000, 95000), lm_pres,
                            lm_data, LM_N) == doctest::Approx(298.5f));
}

// Counts heap allocations, so the threshold-layer tests can show that the
// walker never allocates.
static std::size_t heap_allocations = 0;

void* operator new(std::size_t size) {
    ++heap_allocations;
    if (void* ptr = std::malloc(size ? size : 1)) return ptr;
    throw std::bad_alloc();
}
void operator delete(void* ptr) noexcept { std::free(ptr); }
void operator delete(void* ptr, std::size_t) noexcept { std::free(ptr); }

// What one sharp::for_each_threshold_layer walk reported. Fixed-size
// storage, so recording it doesn't allocate either.
struct ThresholdWalk {
    struct Layer {
        float bottom;
        float top;
        bool above;
        float pos_area;
        float neg_area;
    };
    static constexpr std::ptrdiff_t capacity = 32;
    Layer layers[capacity];
    std::ptrdiff_t num_layers = 0;
    std::ptrdiff_t reads[capacity];  // levels passed to the accessor, in order
    std::ptrdiff_t num_reads = 0;
    std::ptrdiff_t returned = 0;
};

// Walks data[] against the threshold and records every layer and level
// read. With stop_after > 0, the callback stops the walk at that layer.
static ThresholdWalk walk_threshold(const float height[], const float data[],
                                    const std::ptrdiff_t N,
                                    const float threshold,
                                    const float min_depth = 0.0f,
                                    const float min_area = 0.0f,
                                    const std::ptrdiff_t stop_after = 0) {
    ThresholdWalk walk;
    auto data_at = [&](const std::ptrdiff_t k) {
        if (walk.num_reads < walk.capacity) walk.reads[walk.num_reads] = k;
        ++walk.num_reads;
        return data[k];
    };
    auto on_layer = [&](const sharp::HeightLayer& layer, const bool above,
                        const float pos_area, const float neg_area) {
        if (walk.num_layers < walk.capacity) {
            walk.layers[walk.num_layers] = {layer.bottom, layer.top, above,
                                            pos_area, neg_area};
        }
        ++walk.num_layers;
        return walk.num_layers != stop_after;
    };
    const std::size_t heap_before = heap_allocations;
    walk.returned = sharp::for_each_threshold_layer(
        height, data_at, N, threshold, min_depth, min_area, on_layer);
    CHECK(heap_allocations == heap_before);
    CHECK(walk.returned == walk.num_layers);
    REQUIRE(walk.num_layers <= walk.capacity);
    return walk;
}

// Checks layer i of a walk against its expected bounds, side, and areas.
static void check_layer(const ThresholdWalk& walk, const std::ptrdiff_t i,
                        const float bottom, const float top, const bool above,
                        const float pos_area, const float neg_area) {
    INFO("layer ", i);
    REQUIRE(i < walk.num_layers);
    const ThresholdWalk::Layer& layer = walk.layers[i];
    CHECK(layer.bottom == doctest::Approx(bottom));
    CHECK(layer.top == doctest::Approx(top));
    CHECK(layer.above == above);
    CHECK(layer.pos_area == doctest::Approx(pos_area));
    CHECK(layer.neg_area == doctest::Approx(neg_area));
}

// Checks that walking -data against 0 gives the same layers as data, with
// the sides and areas swapped.
static void check_mirror(const float height[], const float data[],
                         const std::ptrdiff_t N, const float min_depth,
                         const float min_area) {
    float mirrored[ThresholdWalk::capacity];
    REQUIRE(N <= ThresholdWalk::capacity);
    for (std::ptrdiff_t k = 0; k < N; ++k) mirrored[k] = -data[k];
    const ThresholdWalk walk =
        walk_threshold(height, data, N, 0.0f, min_depth, min_area);
    const ThresholdWalk flip =
        walk_threshold(height, mirrored, N, 0.0f, min_depth, min_area);
    REQUIRE(flip.num_layers == walk.num_layers);
    for (std::ptrdiff_t i = 0; i < walk.num_layers; ++i) {
        INFO("layer ", i);
        CHECK(flip.layers[i].bottom == walk.layers[i].bottom);
        CHECK(flip.layers[i].top == walk.layers[i].top);
        CHECK(flip.layers[i].above == !walk.layers[i].above);
        CHECK(flip.layers[i].pos_area == -walk.layers[i].neg_area);
        CHECK(flip.layers[i].neg_area == -walk.layers[i].pos_area);
    }
}

TEST_CASE("Testing for_each_threshold_layer crossings and areas") {
    constexpr std::ptrdiff_t N = 4;
    constexpr float hght[N] = {0.0f, 1000.0f, 2000.0f, 3000.0f};

    // Crossings at 500 m and 2000 + 1000/3 m. Each area is a sum of
    // trapezoids split at the crossings.
    constexpr float data[N] = {2.0f, -2.0f, -2.0f, 4.0f};
    const float z_cross = 2000.0f + 1000.0f / 3.0f;
    const float neg_area = -500.0f - 2000.0f - 1000.0f / 3.0f;
    const float pos_area = 2.0f * (3000.0f - z_cross);

    const ThresholdWalk walk = walk_threshold(hght, data, N, 0.0f);
    CHECK(walk.num_layers == 3);
    check_layer(walk, 0, 0.0f, 500.0f, true, 500.0f, 0.0f);
    check_layer(walk, 1, 500.0f, z_cross, false, 0.0f, neg_area);
    check_layer(walk, 2, z_cross, 3000.0f, true, pos_area, 0.0f);

    // The same profile as relative humidity against 0.75: the same
    // crossings, and areas scaled by 0.1.
    constexpr float relh[N] = {0.95f, 0.55f, 0.55f, 1.15f};
    const ThresholdWalk rh_walk = walk_threshold(hght, relh, N, 0.75f);
    CHECK(rh_walk.num_layers == 3);
    check_layer(rh_walk, 0, 0.0f, 500.0f, true, 50.0f, 0.0f);
    check_layer(rh_walk, 1, 500.0f, z_cross, false, 0.0f, 0.1f * neg_area);
    check_layer(rh_walk, 2, z_cross, 3000.0f, true, 0.1f * pos_area, 0.0f);
}

TEST_CASE("Testing for_each_threshold_layer levels at the threshold") {
    constexpr std::ptrdiff_t N = 4;
    constexpr float hght[N] = {0.0f, 100.0f, 200.0f, 300.0f};

    // The first run takes the side of the first level off the threshold.
    {
        INFO("first level at the threshold");
        constexpr float data[N] = {0.0f, 1.0f, -1.0f, -1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, N, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 150.0f, true, 75.0f, 0.0f);
        check_layer(walk, 1, 150.0f, 300.0f, false, 0.0f, -125.0f);
    }
    {
        INFO("first two levels at the threshold");
        constexpr float data[N] = {0.0f, 0.0f, -2.0f, 2.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, N, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 250.0f, false, 0.0f, -150.0f);
        check_layer(walk, 1, 250.0f, 300.0f, true, 50.0f, 0.0f);
    }
    {
        INFO("one level at the threshold between opposite sides");
        constexpr float data[N] = {1.0f, 0.0f, -1.0f, -1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, N, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 100.0f, true, 50.0f, 0.0f);
        check_layer(walk, 1, 100.0f, 300.0f, false, 0.0f, -150.0f);
    }
    {
        INFO("every level at the threshold");
        constexpr float relh[3] = {0.75f, 0.75f, 0.75f};
        const ThresholdWalk walk = walk_threshold(hght, relh, 3, 0.75f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 200.0f, false, 0.0f, 0.0f);

        const ThresholdWalk noisy =
            walk_threshold(hght, relh, 3, 0.75f, 5000.0f, 1.0f);
        CHECK(noisy.num_layers == 1);
        check_layer(noisy, 0, 0.0f, 200.0f, false, 0.0f, 0.0f);
    }
}

TEST_CASE("Testing for_each_threshold_layer plateaus at the threshold") {
    constexpr std::ptrdiff_t N = 4;
    constexpr float hght[N] = {0.0f, 100.0f, 2000.0f, 2100.0f};

    // Between runs on the same side, a plateau continues the run: one
    // 2100 m moist layer, and one cold layer for a 0 C plateau.
    {
        INFO("moist, 75%, 75%, moist");
        constexpr float relh[N] = {0.8f, 0.75f, 0.75f, 0.8f};
        const ThresholdWalk walk = walk_threshold(hght, relh, N, 0.75f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2100.0f, true, 5.0f, 0.0f);
    }
    {
        INFO("cold, 0 C, 0 C, cold");
        constexpr float tmpk[N] = {272.15f, 273.15f, 273.15f, 272.15f};
        const ThresholdWalk walk = walk_threshold(hght, tmpk, N, 273.15f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2100.0f, false, 0.0f, -100.0f);
    }

    // Between opposite sides, the plateau stays with the run it continues,
    // and the crossing is where the data leave the threshold.
    {
        INFO("above, plateau, below");
        constexpr float data[N] = {1.0f, 0.0f, 0.0f, -1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, N, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 2000.0f, true, 50.0f, 0.0f);
        check_layer(walk, 1, 2000.0f, 2100.0f, false, 0.0f, -50.0f);
    }
    {
        INFO("below, plateau, above");
        constexpr float data[N] = {-1.0f, 0.0f, 0.0f, 1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, N, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 2000.0f, false, 0.0f, -50.0f);
        check_layer(walk, 1, 2000.0f, 2100.0f, true, 50.0f, 0.0f);
    }

    // The plateau's depth counts toward its run: 100 m above plus a 1900 m
    // plateau makes a 2000 m run, deep enough to absorb the 500 m below
    // run. Without the plateau, the deeper below run would win the column.
    {
        INFO("plateau depth counts toward min_depth");
        constexpr float hght_2[N] = {0.0f, 100.0f, 2000.0f, 2500.0f};
        constexpr float data[N] = {1.0f, 0.0f, 0.0f, -1.0f};
        const ThresholdWalk walk =
            walk_threshold(hght_2, data, N, 0.0f, 1000.0f, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2500.0f, true, 50.0f, -250.0f);
    }
}

TEST_CASE("Testing for_each_threshold_layer with 0, 1, and 2 levels") {
    constexpr float hght[2] = {500.0f, 600.0f};
    {
        INFO("N = 0");
        constexpr float data[2] = {1.0f, -1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 0, 0.0f);
        CHECK(walk.num_layers == 0);
        CHECK(walk.num_reads == 0);
    }
    // A lone level has no depth and no side, whatever its value.
    for (const float value : {3.0f, -3.0f, 0.0f}) {
        INFO("N = 1, value ", value);
        const float data[1] = {value};
        const ThresholdWalk walk = walk_threshold(hght, data, 1, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 500.0f, 500.0f, false, 0.0f, 0.0f);
    }
    {
        INFO("N = 2, one crossing");
        constexpr float data[2] = {1.0f, -1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 2, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 500.0f, 550.0f, true, 25.0f, 0.0f);
        check_layer(walk, 1, 550.0f, 600.0f, false, 0.0f, -25.0f);
    }
    {
        INFO("N = 2, no crossing");
        constexpr float data[2] = {1.0f, 3.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 2, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 500.0f, 600.0f, true, 200.0f, 0.0f);
    }
    {
        INFO("N = 2, first level at the threshold");
        constexpr float data[2] = {0.0f, -2.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 2, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 500.0f, 600.0f, false, 0.0f, -100.0f);
    }
    {
        INFO("N = 2, both levels at the threshold");
        constexpr float data[2] = {0.0f, 0.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 2, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 500.0f, 600.0f, false, 0.0f, 0.0f);
    }
}

#ifndef NO_QC
TEST_CASE("Testing for_each_threshold_layer bridges MISSING and NaN") {
    constexpr float MISSING = sharp::MISSING;

    // Each profile reduces to 2 at 0 m and -2 at 1000 m.
    {
        INFO("missing and NaN values");
        constexpr float hght[5] = {0.0f, 250.0f, 500.0f, 750.0f, 1000.0f};
        constexpr float data[5] = {2.0f, MISSING, nanval, MISSING, -2.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 5, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 500.0f, true, 500.0f, 0.0f);
        check_layer(walk, 1, 500.0f, 1000.0f, false, 0.0f, -500.0f);
    }
    {
        INFO("missing and NaN heights");
        constexpr float hght[4] = {0.0f, MISSING, nanval, 1000.0f};
        constexpr float data[4] = {2.0f, 5.0f, 5.0f, -2.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 4, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 500.0f, true, 500.0f, 0.0f);
        check_layer(walk, 1, 500.0f, 1000.0f, false, 0.0f, -500.0f);
    }
    {
        INFO("missing levels at the bottom and top");
        constexpr float hght[5] = {nanval, -500.0f, 0.0f, 1000.0f, MISSING};
        constexpr float data[5] = {5.0f, MISSING, 2.0f, -2.0f, 5.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 5, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 500.0f, true, 500.0f, 0.0f);
        check_layer(walk, 1, 500.0f, 1000.0f, false, 0.0f, -500.0f);

        constexpr float data_2[5] = {5.0f, nanval, 2.0f, -2.0f, nanval};
        const ThresholdWalk walk_2 = walk_threshold(hght, data_2, 5, 0.0f);
        CHECK(walk_2.num_layers == 2);
        check_layer(walk_2, 0, 0.0f, 500.0f, true, 500.0f, 0.0f);
        check_layer(walk_2, 1, 500.0f, 1000.0f, false, 0.0f, -500.0f);
    }
    {
        INFO("no valid level");
        constexpr float hght[3] = {0.0f, 100.0f, MISSING};
        constexpr float data[3] = {MISSING, nanval, 1.0f};
        const ThresholdWalk walk = walk_threshold(hght, data, 3, 0.0f);
        CHECK(walk.num_layers == 0);
    }
    {
        INFO("one valid level");
        constexpr float hght[3] = {0.0f, 100.0f, 200.0f};
        constexpr float data[3] = {MISSING, 5.0f, nanval};
        const ThresholdWalk walk = walk_threshold(hght, data, 3, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 100.0f, 100.0f, false, 0.0f, 0.0f);
    }
}
#endif

TEST_CASE("Testing for_each_threshold_layer reads each level once") {
    constexpr std::ptrdiff_t N = 8;
    constexpr float hght[N] = {0.0f,   100.0f, 200.0f, 300.0f,
                               400.0f, 500.0f, 600.0f, 700.0f};
    constexpr float data[N] = {1.0f, -1.0f, -2.0f, 1.0f,
                               0.0f, 0.0f,  -1.0f, 2.0f};
    for (const float min_depth : {0.0f, 150.0f, 1000.0f}) {
        INFO("min_depth ", min_depth);
        const ThresholdWalk walk =
            walk_threshold(hght, data, N, 0.0f, min_depth, 0.0f);
        REQUIRE(walk.num_reads == N);
        for (std::ptrdiff_t k = 0; k < N; ++k) CHECK(walk.reads[k] == k);
    }
}

TEST_CASE("Testing for_each_threshold_layer early stop") {
    // Crossings at 50, 150, 250, and 350 m.
    constexpr std::ptrdiff_t N = 5;
    constexpr float hght[N] = {0.0f, 100.0f, 200.0f, 300.0f, 400.0f};
    constexpr float data[N] = {1.0f, -1.0f, 1.0f, -1.0f, 1.0f};

    const ThresholdWalk full = walk_threshold(hght, data, N, 0.0f);
    CHECK(full.num_layers == 5);

    // With no noise thresholds, a layer is reported at the first level past
    // its top crossing, and the walk reads nothing after a stop.
    const ThresholdWalk first =
        walk_threshold(hght, data, N, 0.0f, 0.0f, 0.0f, 1);
    CHECK(first.returned == 1);
    check_layer(first, 0, 0.0f, 50.0f, true, 25.0f, 0.0f);
    CHECK(first.num_reads == 2);

    const ThresholdWalk third =
        walk_threshold(hght, data, N, 0.0f, 0.0f, 0.0f, 3);
    CHECK(third.returned == 3);
    check_layer(third, 2, 150.0f, 250.0f, true, 50.0f, 0.0f);
    CHECK(third.num_reads == 4);

    // With min_depth = 100, the 50 m runs at the ends join their neighbors.
    // The first layer is reported once the 150-250 m run reaches 100 m,
    // at its top crossing, which is found at the 300 m level.
    const ThresholdWalk deep = walk_threshold(hght, data, N, 0.0f, 100.0f);
    CHECK(deep.num_layers == 3);
    check_layer(deep, 0, 0.0f, 150.0f, false, 25.0f, -50.0f);
    check_layer(deep, 1, 150.0f, 250.0f, true, 50.0f, 0.0f);
    check_layer(deep, 2, 250.0f, 400.0f, false, 25.0f, -50.0f);

    const ThresholdWalk deep_first =
        walk_threshold(hght, data, N, 0.0f, 100.0f, 0.0f, 1);
    CHECK(deep_first.returned == 1);
    check_layer(deep_first, 0, 0.0f, 150.0f, false, 25.0f, -50.0f);
    CHECK(deep_first.num_reads == 4);
}

// The merge-rule profiles use values of +/-1 on either side of each
// crossing, 50 m from it, so every crossing lands on a round height.
TEST_CASE("Testing for_each_threshold_layer merge rule") {
    {
        // Runs: above 0-1000, then a zone of below 1000-1100, above
        // 1100-1200, and below 1200-1300, then above 1300-2300. The zone is
        // mostly below, but it lies between two above runs.
        INFO("same-side zone is absorbed");
        constexpr std::ptrdiff_t N = 7;
        constexpr float hght[N] = {0.0f,    950.0f,  1050.0f, 1150.0f,
                                   1250.0f, 1350.0f, 2300.0f};
        constexpr float data[N] = {1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 1.0f};
        const ThresholdWalk walk =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2300.0f, true, 2000.0f, -100.0f);
        check_mirror(hght, data, N, 500.0f, 0.0f);

        // With both thresholds at 0, the layers are the raw crossings.
        const ThresholdWalk raw = walk_threshold(hght, data, N, 0.0f);
        CHECK(raw.num_layers == 5);
        check_layer(raw, 0, 0.0f, 1000.0f, true, 975.0f, 0.0f);
        check_layer(raw, 1, 1000.0f, 1100.0f, false, 0.0f, -50.0f);
        check_layer(raw, 2, 1100.0f, 1200.0f, true, 50.0f, 0.0f);
        check_layer(raw, 3, 1200.0f, 1300.0f, false, 0.0f, -50.0f);
        check_layer(raw, 4, 1300.0f, 2300.0f, true, 975.0f, 0.0f);
    }
    {
        // Runs: above 0-1000, below 1000-1100, above 1100-1400, below
        // 1400-2400. The zone is 300 m above, 100 m below, so it joins the
        // lower layer.
        INFO("zone majority to the lower layer");
        constexpr std::ptrdiff_t N = 7;
        constexpr float hght[N] = {0.0f,    950.0f,  1050.0f, 1150.0f,
                                   1350.0f, 1450.0f, 2400.0f};
        constexpr float data[N] = {1.0f, 1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f};
        const ThresholdWalk walk =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 1400.0f, true, 1225.0f, -50.0f);
        check_layer(walk, 1, 1400.0f, 2400.0f, false, 0.0f, -975.0f);
        check_mirror(hght, data, N, 500.0f, 0.0f);
    }
    {
        // Runs: above 0-1000, below 1000-1300, above 1300-1400, below
        // 1400-2400. The zone is 100 m above, 300 m below, so it joins the
        // upper layer.
        INFO("zone majority to the upper layer");
        constexpr std::ptrdiff_t N = 7;
        constexpr float hght[N] = {0.0f,    950.0f,  1050.0f, 1250.0f,
                                   1350.0f, 1450.0f, 2400.0f};
        constexpr float data[N] = {1.0f, 1.0f,  -1.0f, -1.0f,
                                   1.0f, -1.0f, -1.0f};
        const ThresholdWalk walk =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 1000.0f, true, 975.0f, 0.0f);
        check_layer(walk, 1, 1000.0f, 2400.0f, false, 50.0f, -1225.0f);
        check_mirror(hght, data, N, 500.0f, 0.0f);
    }
    {
        // Runs: above 0-1000, below 1000-1200, above 1200-1400, below
        // 1400-2400. A 200 m / 200 m tie goes to the lower layer, and
        // still does with the signs swapped.
        INFO("zone tie to the lower layer");
        constexpr std::ptrdiff_t N = 8;
        constexpr float hght[N] = {0.0f,    950.0f,  1050.0f, 1150.0f,
                                   1250.0f, 1350.0f, 1450.0f, 2400.0f};
        constexpr float data[N] = {1.0f, 1.0f, -1.0f, -1.0f,
                                   1.0f, 1.0f, -1.0f, -1.0f};
        const ThresholdWalk walk =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(walk.num_layers == 2);
        check_layer(walk, 0, 0.0f, 1400.0f, true, 1125.0f, -150.0f);
        check_layer(walk, 1, 1400.0f, 2400.0f, false, 0.0f, -975.0f);
        check_mirror(hght, data, N, 500.0f, 0.0f);
    }
    {
        // Runs: below 0-300, above 300-1300, below 1300-1700. The zones at
        // the bottom and top join the one significant run, even though
        // they are on the other side.
        INFO("zones at the bottom and top");
        constexpr std::ptrdiff_t N = 6;
        constexpr float hght[N] = {0.0f,    250.0f,  350.0f,
                                   1250.0f, 1350.0f, 1700.0f};
        constexpr float data[N] = {-1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f};
        const ThresholdWalk walk =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 1700.0f, true, 950.0f, -650.0f);
        check_mirror(hght, data, N, 500.0f, 0.0f);
    }
    {
        // Runs: below 0-300, above 300-400, below 400-700, above 700-1700.
        // The three runs under 700 m are 300, 100, and 300 m deep, with
        // |area| 275, 50, and 250. Merged into one below stretch they would
        // pass min_depth = 500 or min_area = 400, but only a run's own size
        // counts, so all three join the one significant run above them.
        INFO("a run's size ignores what it absorbed");
        constexpr std::ptrdiff_t N = 7;
        constexpr float hght[N] = {0.0f,   250.0f, 350.0f, 450.0f,
                                   650.0f, 750.0f, 1700.0f};
        constexpr float data[N] = {-1.0f, -1.0f, 1.0f, -1.0f,
                                   -1.0f, 1.0f,  1.0f};
        const ThresholdWalk by_depth =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(by_depth.num_layers == 1);
        check_layer(by_depth, 0, 0.0f, 1700.0f, true, 1025.0f, -525.0f);

        const ThresholdWalk by_area =
            walk_threshold(hght, data, N, 0.0f, 0.0f, 400.0f);
        CHECK(by_area.num_layers == 1);
        check_layer(by_area, 0, 0.0f, 1700.0f, true, 1025.0f, -525.0f);
    }
}

TEST_CASE("Testing for_each_threshold_layer with no significant run") {
    // Moist 0-1000 m and dry 1000-2000 m: an exact 1000 m / 1000 m tie.
    // The column takes the side of its lowest run.
    constexpr std::ptrdiff_t N = 4;
    constexpr float hght[N] = {0.0f, 500.0f, 1500.0f, 2000.0f};
    {
        INFO("moist below dry");
        constexpr float relh[N] = {0.85f, 0.85f, 0.65f, 0.65f};
        const ThresholdWalk walk =
            walk_threshold(hght, relh, N, 0.75f, 1500.0f, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2000.0f, true, 75.0f, -75.0f);

        // min_depth is inclusive: at 1000 m both runs are significant.
        const ThresholdWalk at_depth =
            walk_threshold(hght, relh, N, 0.75f, 1000.0f, 0.0f);
        CHECK(at_depth.num_layers == 2);
        check_layer(at_depth, 0, 0.0f, 1000.0f, true, 75.0f, 0.0f);
        check_layer(at_depth, 1, 1000.0f, 2000.0f, false, 0.0f, -75.0f);
    }
    {
        INFO("dry below moist");
        constexpr float relh[N] = {0.65f, 0.65f, 0.85f, 0.85f};
        const ThresholdWalk walk =
            walk_threshold(hght, relh, N, 0.75f, 1500.0f, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2000.0f, false, 75.0f, -75.0f);
    }
    {
        // Moist 0-1000 m, dry 1000-2500 m: the majority side wins.
        INFO("dry majority");
        constexpr float hght_2[N] = {0.0f, 500.0f, 1500.0f, 2500.0f};
        constexpr float relh[N] = {0.85f, 0.85f, 0.65f, 0.65f};
        const ThresholdWalk walk =
            walk_threshold(hght_2, relh, N, 0.75f, 2000.0f, 0.0f);
        CHECK(walk.num_layers == 1);
        check_layer(walk, 0, 0.0f, 2500.0f, false, 75.0f, -125.0f);
    }
}

TEST_CASE("Testing for_each_threshold_layer depth and area thresholds") {
    {
        // Above 0-1000, a deep but weak below run 1000-2000 (area -95),
        // above 2000-3000. It passes min_depth but not min_area.
        INFO("deep, weak run");
        constexpr std::ptrdiff_t N = 8;
        constexpr float hght[N] = {0.0f,    900.0f,  950.0f,  1050.0f,
                                   1950.0f, 2050.0f, 2100.0f, 3000.0f};
        constexpr float data[N] = {1.0f,  1.0f, 0.1f, -0.1f,
                                   -0.1f, 0.1f, 1.0f, 1.0f};
        const ThresholdWalk both =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 200.0f);
        CHECK(both.num_layers == 1);
        check_layer(both, 0, 0.0f, 3000.0f, true, 1860.0f, -95.0f);

        const ThresholdWalk depth_only =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 0.0f);
        CHECK(depth_only.num_layers == 3);
        check_layer(depth_only, 0, 0.0f, 1000.0f, true, 930.0f, 0.0f);
        check_layer(depth_only, 1, 1000.0f, 2000.0f, false, 0.0f, -95.0f);
        check_layer(depth_only, 2, 2000.0f, 3000.0f, true, 930.0f, 0.0f);
    }
    {
        // Above 0-1000, a shallow but strong below run 1000-1100 (area
        // -500), above 1100-2100. It passes min_area but not min_depth.
        INFO("shallow, strong run");
        constexpr std::ptrdiff_t N = 5;
        constexpr float hght[N] = {0.0f, 950.0f, 1050.0f, 1150.0f, 2100.0f};
        constexpr float data[N] = {10.0f, 10.0f, -10.0f, 10.0f, 10.0f};
        const ThresholdWalk both =
            walk_threshold(hght, data, N, 0.0f, 500.0f, 200.0f);
        CHECK(both.num_layers == 1);
        check_layer(both, 0, 0.0f, 2100.0f, true, 19500.0f, -500.0f);

        const ThresholdWalk area_only =
            walk_threshold(hght, data, N, 0.0f, 0.0f, 200.0f);
        CHECK(area_only.num_layers == 3);
        check_layer(area_only, 0, 0.0f, 1000.0f, true, 9750.0f, 0.0f);
        check_layer(area_only, 1, 1000.0f, 1100.0f, false, 0.0f, -500.0f);
        check_layer(area_only, 2, 1100.0f, 2100.0f, true, 9750.0f, 0.0f);
    }
}
