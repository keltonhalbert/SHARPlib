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

// Layer conversions with MISSING data at an end level (SHARPlib-mut). An
// endpoint with no valid data level on its open side interpolates to
// MISSING. The conversion now returns the {MISSING, MISSING} layer it returns
// for a layer outside the profile, and an AGL conversion without a valid
// height[0] returns it too. Each "was" comment is the output measured before
// this change. NO_QC builds don't change.
#ifndef NO_QC
template <typename L>
static void check_missing_layer(const L layer) {
    CHECK(layer.bottom == MISSING);
    CHECK(layer.top == MISSING);
}

template <typename L>
static void check_layer_bounds(const L layer, const float bottom,
                               const float top) {
    CHECK(layer.bottom == doctest::Approx(bottom));
    CHECK(layer.top == doctest::Approx(top));
}

// The bead's profile, with MISSING pressure or height at the first, last, or
// an interior level.
constexpr std::ptrdiff_t ME_N = 5;
constexpr float me_hght[ME_N] = {0, 500, 1000, 1500, 2000};
constexpr float me_data[ME_N] = {300, 297, 294, 291, 288};
constexpr float me_pres[ME_N] = {100000, 95000, 90000, 85000, 80000};
constexpr float me_pres_bot[ME_N] = {MISSING, 95000, 90000, 85000, 80000};
constexpr float me_pres_top[ME_N] = {100000, 95000, 90000, 85000, MISSING};
constexpr float me_pres_mid[ME_N] = {100000, 95000, MISSING, 85000, 80000};
constexpr float me_hght_bot[ME_N] = {MISSING, 500, 1000, 1500, 2000};
constexpr float me_hght_top[ME_N] = {0, 500, 1000, 1500, MISSING};
constexpr float me_hght_mid[ME_N] = {0, 500, MISSING, 1500, 2000};

TEST_CASE("Testing height_layer_to_pressure with MISSING end pressure") {
    for (const float* pres : {me_pres_bot, me_pres_top}) {
        CAPTURE(pres[0]);
        for (const bool agl : {false, true}) {
            CAPTURE(agl);
            // was std::range_error
            check_missing_layer(sharp::height_layer_to_pressure(
                {0, 2000}, pres, me_hght, ME_N, agl));
            // was std::range_error
            check_missing_layer(sharp::height_layer_to_pressure(
                {250, 1750}, pres, me_hght, ME_N, agl));
            // a layer that doesn't need the MISSING level, unchanged
            check_layer_bounds(sharp::height_layer_to_pressure(
                                   {500, 1500}, pres, me_hght, ME_N, agl),
                               95000, 85000);
        }
    }
    // an interior MISSING level is bridged, unchanged
    check_layer_bounds(
        sharp::height_layer_to_pressure({0, 2000}, me_pres_mid, me_hght, ME_N),
        100000, 80000);
    check_layer_bounds(sharp::height_layer_to_pressure({250, 1750}, me_pres_mid,
                                                       me_hght, ME_N),
                       97500, 82500);
}

TEST_CASE("Testing pressure_layer_to_height with MISSING end height") {
    // was std::range_error, both AGL flags
    check_missing_layer(sharp::pressure_layer_to_height({97500, 82500}, me_pres,
                                                        me_hght_top, ME_N));
    check_missing_layer(sharp::pressure_layer_to_height(
        {97500, 82500}, me_pres, me_hght_top, ME_N, true));
    check_missing_layer(sharp::pressure_layer_to_height(
        {100000, 80000}, me_pres, me_hght_top, ME_N, true));
    check_missing_layer(sharp::pressure_layer_to_height(
        {100000, 80000}, me_pres, me_hght_bot, ME_N));
    check_missing_layer(sharp::pressure_layer_to_height({97500, 82500}, me_pres,
                                                        me_hght_bot, ME_N));

    // layers that don't need the MISSING level, unchanged
    check_layer_bounds(sharp::pressure_layer_to_height({95000, 85000}, me_pres,
                                                       me_hght_top, ME_N, true),
                       500, 1500);
    check_layer_bounds(sharp::pressure_layer_to_height({95000, 85000}, me_pres,
                                                       me_hght_bot, ME_N),
                       500, 1500);

    // an interior MISSING level is bridged, unchanged
    for (const bool agl : {false, true}) {
        CAPTURE(agl);
        check_layer_bounds(
            sharp::pressure_layer_to_height({100000, 80000}, me_pres,
                                            me_hght_mid, ME_N, agl),
            0, 2000);
        check_layer_bounds(sharp::pressure_layer_to_height(
                               {97500, 82500}, me_pres, me_hght_mid, ME_N, agl),
                           246.794525f, 1746.21484f);
    }
}

TEST_CASE("Testing layer conversion change classes") {
    // (b) both endpoints MISSING, no AGL flag: already the sentinel
    constexpr float pres_mm[ME_N] = {100000, 95000, 90000, MISSING, MISSING};
    check_missing_layer(
        sharp::height_layer_to_pressure({1600, 1900}, pres_mm, me_hght, ME_N));
    constexpr float hght_mm[ME_N] = {300, 800, 1300, MISSING, MISSING};
    check_missing_layer(sharp::pressure_layer_to_height({85000, 80000}, me_pres,
                                                        hght_mm, ME_N));

    // (c) toAGL, both endpoints MISSING: was (-10299, -10299)
    check_missing_layer(sharp::pressure_layer_to_height({85000, 80000}, me_pres,
                                                        hght_mm, ME_N, true));
    // (d) toAGL, height[0] MISSING makes one endpoint MISSING: was (0, 11999)
    check_missing_layer(sharp::pressure_layer_to_height(
        {100000, 80000}, me_pres, me_hght_bot, ME_N, true));
    // (e) toAGL, height[0] MISSING, both endpoints valid: was (10499, 11499)
    check_missing_layer(sharp::pressure_layer_to_height(
        {95000, 85000}, me_pres, me_hght_bot, ME_N, true));
    // (f) toAGL, height[0] NaN: was std::range_error
    constexpr float hght_nan[ME_N] = {nanval, 500, 1000, 1500, 2000};
    check_missing_layer(sharp::pressure_layer_to_height({95000, 85000}, me_pres,
                                                        hght_nan, ME_N, true));

    // (g) isAGL, height[0] MISSING: was (99761.9, 99285.6)
    check_missing_layer(sharp::height_layer_to_pressure(
        {500, 1500}, me_pres, me_hght_bot, ME_N, true));
    // was std::range_error
    check_missing_layer(sharp::height_layer_to_pressure(
        {0, 1000}, me_pres, me_hght_bot, ME_N, true));
    // isAGL, height[0] NaN: was already the sentinel
    check_missing_layer(sharp::height_layer_to_pressure({500, 1500}, me_pres,
                                                        hght_nan, ME_N, true));

    // (h) one endpoint MISSING with an AGL flag and a valid origin. The
    // endpoints interpolate to (800, MISSING); was (500, -10299), which threw
    // std::range_error.
    check_missing_layer(sharp::pressure_layer_to_height({95000, 85000}, me_pres,
                                                        hght_mm, ME_N, true));
    // was std::range_error
    constexpr float hght_300[ME_N] = {300, 800, 1300, 1800, 2300};
    check_missing_layer(sharp::height_layer_to_pressure(
        {250, 1750}, me_pres_top, hght_300, ME_N, true));
}

TEST_CASE("Testing layer_mean with MISSING end pressure") {
    for (const float* pres : {me_pres_bot, me_pres_top}) {
        CAPTURE(pres[0]);
        // was std::range_error
        CHECK(sharp::layer_mean(sharp::HeightLayer(0, 2000), me_hght, pres,
                                me_data, ME_N) == MISSING);
        // was std::range_error
        CHECK(sharp::layer_mean(sharp::HeightLayer(250, 1750), me_hght, pres,
                                me_data, ME_N) == MISSING);
        // clipped to the profile; was std::range_error
        CHECK(sharp::layer_mean(sharp::HeightLayer(-500, 3000), me_hght, pres,
                                me_data, ME_N) == MISSING);
    }
    // Pressure is the integration coordinate, so a MISSING pressure[0] makes
    // even this layer MISSING and a MISSING pressure[N-1] doesn't. Both
    // unchanged.
    CHECK(sharp::layer_mean(sharp::HeightLayer(500, 1500), me_hght, me_pres_bot,
                            me_data, ME_N) == MISSING);
    CHECK(sharp::layer_mean(sharp::HeightLayer(500, 1500), me_hght, me_pres_top,
                            me_data, ME_N) == doctest::Approx(294.0f));
}
#endif

