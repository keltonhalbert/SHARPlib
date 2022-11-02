#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>
#include <limits>

#include "constants.h"
#include "doctest.h"
#include "winds.h"

TEST_CASE("Testing wind components (u,v) operations") {
    sharp::WindComponents s_wind = {0.0, 10.0};
    sharp::WindComponents w_wind = {10.0, 0.0};
    sharp::WindComponents n_wind = {0.0, -10.0};
    sharp::WindComponents e_wind = {-10.0, 0};

    sharp::WindComponents sw_wind = {10.0, 10.0};
    sharp::WindComponents se_wind = {-10.0, 10.0};
    sharp::WindComponents ne_wind = {-10.0, -10.0};
    sharp::WindComponents nw_wind = {10.0, -10.0};


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


    CHECK(sharp::vector_magnitude(sw_wind.u, sw_wind.v) 
                         == doctest::Approx(14.142135));
    CHECK(sharp::vector_magnitude(se_wind.u, se_wind.v) 
                         == doctest::Approx(14.142135));
    CHECK(sharp::vector_magnitude(ne_wind.u, ne_wind.v) 
                         == doctest::Approx(14.142135));
    CHECK(sharp::vector_magnitude(nw_wind.u, nw_wind.v) 
                         == doctest::Approx(14.142135));

    sharp::WindVector s_vect = sharp::components_to_vector(s_wind.u, s_wind.v);
    sharp::WindVector w_vect = sharp::components_to_vector(w_wind.u, w_wind.v);
    sharp::WindVector n_vect = sharp::components_to_vector(n_wind.u, n_wind.v);
    sharp::WindVector e_vect = sharp::components_to_vector(e_wind.u, e_wind.v);

    sharp::WindVector sw_vect = sharp::components_to_vector(sw_wind.u, sw_wind.v);
    sharp::WindVector se_vect = sharp::components_to_vector(se_wind.u, se_wind.v);
    sharp::WindVector ne_vect = sharp::components_to_vector(ne_wind.u, ne_wind.v);
    sharp::WindVector nw_vect = sharp::components_to_vector(nw_wind.u, nw_wind.v);

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
    sharp::WindVector s_vect = {10.0, 180.0};
    sharp::WindVector w_vect = {10.0, 270.0};
    sharp::WindVector n_vect = {10.0, 0.0};
    sharp::WindVector e_vect = {10.0, 90.0};

    sharp::WindVector sw_vect = {10.0, 225.0};
    sharp::WindVector nw_vect = {10.0, 315.0};
    sharp::WindVector ne_vect = {10.0, 45.0};
    sharp::WindVector se_vect = {10.0, 135.0};


    // We have to use aproximations because floating point math with
    // sin and cos is not exact, sadly...
    CHECK(sharp::u_component(s_vect.speed, s_vect.direction) 
                                    == doctest::Approx(0.0));
    CHECK(sharp::u_component(w_vect.speed, w_vect.direction) 
                                    == doctest::Approx(10.0));
    CHECK(sharp::u_component(n_vect.speed, n_vect.direction) 
                                    == doctest::Approx(0.0));
    CHECK(sharp::u_component(e_vect.speed, e_vect.direction) 
                                    == doctest::Approx(-10.0));

    CHECK(sharp::v_component(s_vect.speed, s_vect.direction) 
                                    == doctest::Approx(10.0));
    CHECK(sharp::v_component(w_vect.speed, w_vect.direction) 
                                    == doctest::Approx(0.0));
    CHECK(sharp::v_component(n_vect.speed, n_vect.direction) 
                                    == doctest::Approx(-10.0));
    CHECK(sharp::v_component(e_vect.speed, e_vect.direction) 
                                    == doctest::Approx(0.0));

    CHECK(sharp::u_component(sw_vect.speed, sw_vect.direction)
                                    == doctest::Approx(7.07107));
    CHECK(sharp::u_component(nw_vect.speed, nw_vect.direction)
                                    == doctest::Approx(7.07107));
    CHECK(sharp::u_component(ne_vect.speed, ne_vect.direction)
                                    == doctest::Approx(-7.07107));
    CHECK(sharp::u_component(se_vect.speed, se_vect.direction)
                                    == doctest::Approx(-7.07107));

    CHECK(sharp::v_component(sw_vect.speed, sw_vect.direction)
                                    == doctest::Approx(7.07107));
    CHECK(sharp::v_component(nw_vect.speed, nw_vect.direction)
                                    == doctest::Approx(-7.07107));
    CHECK(sharp::v_component(ne_vect.speed, ne_vect.direction)
                                    == doctest::Approx(-7.07107));
    CHECK(sharp::v_component(se_vect.speed, se_vect.direction)
                                    == doctest::Approx(7.07107));
}
