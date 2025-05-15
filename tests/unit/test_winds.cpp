#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/winds.h>

#include <cmath>

#include "doctest.h"

TEST_CASE("Testing wind components (u,v) operations") {
    constexpr sharp::WindComponents s_wind = {0.0, 10.0};
    constexpr sharp::WindComponents w_wind = {10.0, 0.0};
    constexpr sharp::WindComponents n_wind = {0.0, -10.0};
    constexpr sharp::WindComponents e_wind = {-10.0, 0};

    constexpr sharp::WindComponents sw_wind = {10.0, 10.0};
    constexpr sharp::WindComponents se_wind = {-10.0, 10.0};
    constexpr sharp::WindComponents ne_wind = {-10.0, -10.0};
    constexpr sharp::WindComponents nw_wind = {10.0, -10.0};

    CHECK(sharp::vector_angle(s_wind.u, s_wind.v) == 180.0);
    CHECK(sharp::vector_angle(w_wind.u, w_wind.v) == 270.0);
    CHECK(sharp::vector_angle(n_wind.u, n_wind.v) == 0.0);
    CHECK(sharp::vector_angle(e_wind.u, e_wind.v) == 90.0);

    CHECK(sharp::vector_angle(sw_wind.u, sw_wind.v) == 225.0);
    CHECK(sharp::vector_angle(se_wind.u, se_wind.v) == 135.0);
    CHECK(sharp::vector_angle(ne_wind.u, ne_wind.v) == 45.0);
    CHECK(sharp::vector_angle(nw_wind.u, nw_wind.v) == 315.0);

    CHECK(sharp::vector_magnitude(s_wind.u, s_wind.v) == 10.0);
    CHECK(sharp::vector_magnitude(w_wind.u, w_wind.v) == 10.0);
    CHECK(sharp::vector_magnitude(n_wind.u, n_wind.v) == 10.0);
    CHECK(sharp::vector_magnitude(e_wind.u, e_wind.v) == 10.0);

    CHECK(sharp::vector_magnitude(sw_wind.u, sw_wind.v) ==
          doctest::Approx(14.142135));
    CHECK(sharp::vector_magnitude(se_wind.u, se_wind.v) ==
          doctest::Approx(14.142135));
    CHECK(sharp::vector_magnitude(ne_wind.u, ne_wind.v) ==
          doctest::Approx(14.142135));
    CHECK(sharp::vector_magnitude(nw_wind.u, nw_wind.v) ==
          doctest::Approx(14.142135));

    const auto s_vect = sharp::components_to_vector(s_wind.u, s_wind.v);
    const auto w_vect = sharp::components_to_vector(w_wind.u, w_wind.v);
    const auto n_vect = sharp::components_to_vector(n_wind.u, n_wind.v);
    const auto e_vect = sharp::components_to_vector(e_wind.u, e_wind.v);

    const auto sw_vect = sharp::components_to_vector(sw_wind.u, sw_wind.v);
    const auto se_vect = sharp::components_to_vector(se_wind.u, se_wind.v);
    const auto ne_vect = sharp::components_to_vector(ne_wind.u, ne_wind.v);
    const auto nw_vect = sharp::components_to_vector(nw_wind.u, nw_wind.v);

    CHECK(s_vect.speed == 10.0);
    CHECK(w_vect.speed == 10.0);
    CHECK(n_vect.speed == 10.0);
    CHECK(e_vect.speed == 10.0);

    CHECK(s_vect.direction == 180.0);
    CHECK(w_vect.direction == 270);
    CHECK(n_vect.direction == 0);
    CHECK(e_vect.direction == 90);

    CHECK(sw_vect.speed == doctest::Approx(14.142135));
    CHECK(se_vect.speed == doctest::Approx(14.142135));
    CHECK(ne_vect.speed == doctest::Approx(14.142135));
    CHECK(nw_vect.speed == doctest::Approx(14.142135));

    CHECK(sw_vect.direction == 225.0);
    CHECK(se_vect.direction == 135.0);
    CHECK(ne_vect.direction == 45.0);
    CHECK(nw_vect.direction == 315.0);
}

TEST_CASE("Testing vector (speed, direction) operations") {
    constexpr sharp::WindVector s_vect = {10.0, 180.0};
    constexpr sharp::WindVector w_vect = {10.0, 270.0};
    constexpr sharp::WindVector n_vect = {10.0, 0.0};
    constexpr sharp::WindVector e_vect = {10.0, 90.0};

    constexpr sharp::WindVector sw_vect = {10.0, 225.0};
    constexpr sharp::WindVector nw_vect = {10.0, 315.0};
    constexpr sharp::WindVector ne_vect = {10.0, 45.0};
    constexpr sharp::WindVector se_vect = {10.0, 135.0};

    // We have to use aproximations because floating point math with
    // sin and cos is not exact, sadly...
    CHECK(sharp::u_component(s_vect.speed, s_vect.direction) ==
          doctest::Approx(0.0));
    CHECK(sharp::u_component(w_vect.speed, w_vect.direction) ==
          doctest::Approx(10.0));
    CHECK(sharp::u_component(n_vect.speed, n_vect.direction) ==
          doctest::Approx(0.0));
    CHECK(sharp::u_component(e_vect.speed, e_vect.direction) ==
          doctest::Approx(-10.0));

    CHECK(sharp::v_component(s_vect.speed, s_vect.direction) ==
          doctest::Approx(10.0));
    CHECK(sharp::v_component(w_vect.speed, w_vect.direction) ==
          doctest::Approx(0.0));
    CHECK(sharp::v_component(n_vect.speed, n_vect.direction) ==
          doctest::Approx(-10.0));
    CHECK(sharp::v_component(e_vect.speed, e_vect.direction) ==
          doctest::Approx(0.0));

    CHECK(sharp::u_component(sw_vect.speed, sw_vect.direction) ==
          doctest::Approx(7.07107));
    CHECK(sharp::u_component(nw_vect.speed, nw_vect.direction) ==
          doctest::Approx(7.07107));
    CHECK(sharp::u_component(ne_vect.speed, ne_vect.direction) ==
          doctest::Approx(-7.07107));
    CHECK(sharp::u_component(se_vect.speed, se_vect.direction) ==
          doctest::Approx(-7.07107));

    CHECK(sharp::v_component(sw_vect.speed, sw_vect.direction) ==
          doctest::Approx(7.07107));
    CHECK(sharp::v_component(nw_vect.speed, nw_vect.direction) ==
          doctest::Approx(-7.07107));
    CHECK(sharp::v_component(ne_vect.speed, ne_vect.direction) ==
          doctest::Approx(-7.07107));
    CHECK(sharp::v_component(se_vect.speed, se_vect.direction) ==
          doctest::Approx(7.07107));
}

TEST_CASE("Testing mean wind calculations") {
    // create the testing data
    constexpr int NZ = 10;
    float pres[NZ];
    float u_wind[NZ];
    float v_wind[NZ];

    // pressure is always in Pa
    constexpr float snd_bot = 100000.0;
    constexpr float u_start = 10.0;
    constexpr float v_start = 10.0;
    constexpr float dp = 10000.0;

    for (int k = 0; k < NZ; k++) {
        pres[k] = snd_bot - dp * k;
        u_wind[k] = u_start + 5 * k;
        v_wind[k] = v_start + 5 * k;
    }

    sharp::PressureLayer layer1 = {100000.0, 80000.0, -100.0};
    sharp::PressureLayer layer2 = {100000.0, 50000.0, -100.0};
    sharp::PressureLayer layer3 = {100000.0, 10000.0, -100.0};

    const auto mean_layer1 =
        sharp::mean_wind(layer1, pres, u_wind, v_wind, NZ, true);
    const auto mean_layer2 =
        sharp::mean_wind(layer2, pres, u_wind, v_wind, NZ, true);
    const auto mean_layer3 =
        sharp::mean_wind(layer3, pres, u_wind, v_wind, NZ, true);

    const auto mean_layer4 =
        sharp::mean_wind(layer1, pres, u_wind, v_wind, NZ, false);
    const auto mean_layer5 =
        sharp::mean_wind(layer2, pres, u_wind, v_wind, NZ, false);
    const auto mean_layer6 =
        sharp::mean_wind(layer3, pres, u_wind, v_wind, NZ, false);

    CHECK(mean_layer1.u == doctest::Approx(14.7669));
    CHECK(mean_layer1.v == doctest::Approx(14.7669));

    CHECK(mean_layer2.u == doctest::Approx(21.0501));
    CHECK(mean_layer2.v == doctest::Approx(21.0501));

    CHECK(mean_layer3.u == doctest::Approx(26.2741));
    CHECK(mean_layer3.v == doctest::Approx(26.2741));

    CHECK(mean_layer4.u == doctest::Approx(14.9538));
    CHECK(mean_layer4.v == doctest::Approx(14.9538));

    CHECK(mean_layer5.u == doctest::Approx(22.4424));
    CHECK(mean_layer5.v == doctest::Approx(22.4424));

    CHECK(mean_layer6.u == doctest::Approx(32.3939));
    CHECK(mean_layer6.v == doctest::Approx(32.3939));
}
