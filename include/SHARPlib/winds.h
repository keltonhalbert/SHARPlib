/**
 * \file
 * \brief Routines used to compute kinematic attributes and indices<!--
 * --> of vertical sounding profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef __SHARP_WINDS_H__
#define __SHARP_WINDS_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

namespace sharp {

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named floats that represent a wind vector.
 *
 * This is used in the convenience function sharp::components_to_vector
 * for situations where having both the speed and direction are desired
 * from components with a single function call.
 */
struct WindVector {
    /**
     * \brief The magntidue of a wind vector
     */
    float speed;

    /**
     * \brief The direction (degrees from North) of a wind vector.
     */
    float direction;
};

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two names floats that represent<!--
 * --> the components of a wind vector.
 *
 * This is used in the convenience function sharp::vector_to_components
 * for situations where having both the U and V wind components are
 * desired from a wind vector with a single function call.
 *
 */
struct WindComponents {
    /**
     * \brief The zonal (U) component of a wind vector.
     */
    float u;
    /**
     * \brief The meridional (V) component of a wind vector.
     */
    float v;
};

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the zonal (U) wind component from a wind vector.
 *
 * Given the wind speed and direction, compute and return the zonal
 * (U) wind component.
 *
 * \param wind_speed     (distance / time)
 * \param wind_direction (degrees from North)
 * \return u_component   (distance / time)
 */
float u_component(float wind_speed, float wind_direction) noexcept;

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the meridional (V) wind component from a wind vector.
 *
 * Given the wind speed and direction, compute and return the meridional
 * (V) wind component.
 *
 * \param wind_speed     (distance / time)
 * \param wind_direction (degrees from North)
 * \return v_component   (distance / time)
 */
float v_component(float wind_speed, float wind_direction) noexcept;

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the direction (degrees from North) a vector points to.
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the direction (degrees from North) the vector points to.
 *
 * \param u_comp            (distance / time)
 * \param v_comp            (distance / time)
 * \return wind_direction   (degrees from North)
 */
float vector_angle(float u_comp, float v_comp) noexcept;

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the magnitude of a vector given its components.
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the magnitude (distance / time) of the vector.
 *
 * \param u_comp        (distance / time)
 * \param v_comp        (distance / time)
 * \return wind_speed   (distance / time)
 */
float vector_magnitude(float u_comp, float v_comp) noexcept;

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Precisely computes the magnitude of a vector given its components.
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the magnitude (distance / time) of the vector. Instead of using
 * a standard sqrt(u*u + v*v) approach, this uses the c++ Standard Template
 * Library implementation of std::hypot. Using std::hypot is ~20 times slower,
 * but ensures a higher accuracy and checks for overflow/underflow when using
 * extremely large or extremely small values. Chances are, you are fine using
 * the imprecise version as long as values are in the range of
 * 1e-100 < values < 1e100.
 *
 * \param u_comp        (distance / time)
 * \param v_comp        (distance / time)
 * \return wind_speed   (distance / time)
 */
float vector_magnitude_precise(float u_comp, float v_comp) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute wind speed and direction<!--
 * --> from components
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the wind speed (distance / time) and direction (degrees from
 * North) from the components as a struct.
 *
 * \param u_comp                (distance / time)
 * \param v_comp                (distance / time)
 * \return sharp::WindVector    {wind_speed, wind_direction}
 */
WindVector components_to_vector(float u_comp, float v_comp) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute wind speed and direction<!--
 * --> from components
 *
 * Given the components of a vector via sharp::WindComponents, compute
 * and return the wind speed (distance / time) and direction (degrees from
 * North) and return them as sharp::WindVector.
 *
 * \param comp                  sharp::WindComponents {u, v}
 * \return sharp::WindVector    {wind_speed, wind_direction}
 */
WindVector components_to_vector(WindComponents comp) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute the U and V components<!--
 * --> of a wind vector.
 *
 * Given the wind speed and direction of a vector, compute and return the
 * zonal and meridional vector components as a struct.
 *
 * \param wind_speed                (distance / time)
 * \param wind_direction            (degrees from North)
 * \return sharp::WindComponents    {u_comp, v_comp}
 */
WindComponents vector_to_components(float wind_speed,
                                    float wind_direction) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute the U and V components<!--
 * --> of a wind vector.
 *
 * Given the wind speed and direction of a vector via sharp::WindVector,
 * compute and return the zonal and meridional vector components as a struct.
 *
 * \param   vect    sharp::WindVector {speed, direction}
 * \return  sharp::WindComponents    {u_comp, v_comp}
 */
WindComponents vector_to_components(WindVector vect) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the pressure weighted mean wind over the given<!--
 * --> sharp::PressureLayer
 *
 * Computes the pressure weighted mean wind over the given sharp::PressureLayer
 * and arrays of pressure and wind components with a length of num_levs.
 *
 * \param layer     sharp::PressureLayer
 * \param pres      (hPa)
 * \param u_wind    (kts or m/s)
 * \param v_wind    (kts or m/s)
 * \param num_levs  (length of arrays)
 * \return sharp::WindComponents    {mean_u, mean_v}
 */
WindComponents mean_wind(PressureLayer layer, const float* pres,
                         const float* u_wind, const float* v_wind,
                         int num_levs) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the non-pressure weighted mean wind over the given<!--
 * --> sharp::PressureLayer
 *
 * Computes the non-pressure weighted mean wind over the given
 * sharp::PressureLayer and arrays of pressure and wind components
 * with a length of num_levs.
 *
 * \param layer     sharp::PressureLayer
 * \param pres      (hPa)
 * \param u_wind    (kts or m/s)
 * \param v_wind    (kts or m/s)
 * \param num_levs  (length of arrays)
 * \return sharp::WindComponents    {mean_u, mean_v}
 */
WindComponents mean_wind_npw(PressureLayer layer, const float* pres,
                             const float* u_wind, const float* v_wind,
                             int num_levs) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compue the wind shear for a given sharp::PressureLayer
 *
 * Computes the U and V components of the wind shear over a
 * sharp::PressureLayer given the vertical sounding arrays
 * of pressure, u_wind, v_wind, and their length.
 *
 * \param layer     sharp::PressureLayer
 * \param pres      (hPa)
 * \param u_wind    (kts or m/s)
 * \param v_wind    (kts or m/s)
 * \param num_levs  (length of arrays)
 * \return sharp::WindComponents    {shear_u, shear_v}
 */
WindComponents wind_shear(PressureLayer layer, const float* pres,
                          const float* u_wind, const float* v_wind,
                          int num_levs) noexcept;

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compue the wind shear for a given sharp::HeightLayer
 *
 * Computes the U and V components of the wind shear over a
 * sharp::HeightLayer given the vertical sounding arrays
 * of height, u_wind, v_wind, and their length.
 *
 * \param layer     sharp::HeightLayer
 * \param height    (meters)
 * \param u_wind    (kts or m/s)
 * \param v_wind    (kts or m/s)
 * \param num_levs  (length of arrays)
 * \return sharp::WindComponents    {shear_u, shear_v}
 */
WindComponents wind_shear(HeightLayer layer_agl, const float* height,
                          const float* u_wind, const float* v_wind,
                          int num_levs) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Computes the Storm Relative Helicity (SRH) over a given<!--
 * --> sharp::HeightLayer.
 *
 * Computes the Storm Relative Helicity (SRH) over a given sharp::HeightLayer
 * using storm motion vector components stored in sharp::WindComponents. This
 * integration occurs over the given arrays of height, u_wind, and v_wind,
 * with num_levs elements in each. The integration only uses interpolation
 * for the top and bottom of the specified layer.
 *
 * The values in sharp::HeightLayer are expected in units above-ground-level
 * (AGL), and gets converted to above-mean-sea-level (MSL) during the
 * computation.
 *
 * \param layer_agl     sharp::HeightLayer   {bottom_agl, top_agl}
 * \param storm_motion  sharp::WindComonents {storm_u, storm_v}
 * \param height        (meters)
 * \param u_wind        (m/s??)
 * \param v_wind        (m/s??)
 * \param num_levs      (length of arrays)
 *
 */
float helicity(HeightLayer layer_agl, WindComponents storm_motion,
               const float* height, const float* u_wind, const float* v_wind,
               int num_levs) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Computes the Storm Relative Helicity (SRH) over a given<!--
 * --> sharp::PressureLayer.
 *
 * Computes the Storm Relative Helicity (SRH) over a given sharp::PressureLayer
 * using storm motion vector components stored in sharp::WindComponents. This
 * integration occurs over the given arrays of height, u_wind, and v_wind,
 * with num_levs elements in each. The integration only uses interpolation
 * for the top and bottom of the specified layer, and uses raw data levels
 * for everything else.
 *
 * The sharp::PressureLayer implementation of helicity is a wrapper around
 * the sharp::HeightLayer version. The values in sharp::PressureLayer are
 * interpolated to height AGL (m) and then passed to the helicity function
 * that works in height coordinates.
 *
 * \param layer_agl     sharp::PressureLayer   {bottom, top}
 * \param storm_motion  sharp::WindComonents {storm_u, storm_v}
 * \param pressure      (hPa)
 * \param height        (meters)
 * \param u_wind        (m/s??)
 * \param v_wind        (m/s??)
 * \param num_levs      (length of arrays)
 *
 */
float helicity(PressureLayer layer, WindComponents storm_motion,
               const float* pressure, const float* height, const float* u_wind,
               const float* v_wind, int num_levs) noexcept;

}  // end namespace sharp

namespace sharp::exper {}  // end namespace sharp::exper

#endif // __SHARP_WINDS_H__
