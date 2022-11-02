/**
 * \file
 * \brief Routines used to compute kinematic attributes and indices of vertical sounding profiles
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
#pragma once

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
 * \brief A simple structure of two names floats that represent the components of a wind vector. 
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
float u_component(float wind_speed, float wind_direction);


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
float v_component(float wind_speed, float wind_direction);


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
float vector_angle(float u_comp, float v_comp);


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
float vector_magnitude(float u_comp, float v_comp);


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
float vector_magnitude_precise(float u_comp, float v_comp);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute wind speed and direction from components 
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the wind speed (distance / time) and direction (degrees from North)
 * from the components as a struct. 
 *
 * \param u_comp                (distance / time)
 * \param v_comp                (distance / time)
 * \return sharp::WindVector    {wind_speed, wind_direction}
 */
WindVector components_to_vector(float u_comp, float v_comp);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute wind speed and direction from components 
 *
 * Given the components of a vector via sharp::WindComponents, compute
 * and return the wind speed (distance / time) and direction (degrees from North)
 * and return them as sharp::WindVector. Implemented as a wrapper around the
 * sharp::components_to_vector(float, float) function call.  
 *
 * \param u_comp                (distance / time)
 * \param v_comp                (distance / time)
 * \return sharp::WindVector    {wind_speed, wind_direction}
 */
WindVector components_to_vector(const WindComponents& comp);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute the U and V components of a wind vector. 
 *
 * Given the wind speed and direction of a vector, compute and return the 
 * zonal and meridional vector components as a struct. 
 *
 * \param wind_speed                (distance / time)
 * \param wind_direction            (degrees from North)
 * \return sharp::WindComponents    {u_comp, v_comp}
 */
WindComponents vector_to_components(float wind_speed, float wind_direction);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute the U and V components of a wind vector. 
 *
 * Given the wind speed and direction of a vector via sharp::WindVector, 
 * compute and return the zonal and meridional vector components as a struct. 
 * Implemented as a wrapper around sharp::vector_to_components(float, float). 
 *
 * \param wind_speed                (distance / time)
 * \param wind_direction            (degrees from North)
 * \return sharp::WindComponents    {u_comp, v_comp}
 */
WindComponents vector_to_components(const WindVector& vect);


WindComponents vector_to_components(float wind_speed, float wind_direction);
/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the pressure weighted mean wind over the given sharp::PressureLayer 
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
WindComponents mean_wind(const PressureLayer& layer, 
                         const float* pres,   const float* u_wind,
                         const float* v_wind, int num_levs);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the pressure weighted mean wind over the given bottom and top pressure levels 
 *
 * Computes the pressure weighted mean wind over the given bottom and top 
 * pressure levels, and arrays of pressure and wind components with a length 
 * of num_levs. This is a wrapper around the sharp::PressureLayer implementation.
 *
 * \param pressure_bot  (hPa)
 * \param pressure_top  (hPa)
 * \param pres          (hPa)
 * \param u_wind        (kts or m/s)
 * \param v_wind        (kts or m/s)
 * \param num_levs      (length of arrays)
 * \return sharp::WindComponents    {mean_u, mean_v}
 */
WindComponents mean_wind(float pressure_bot,  float pressure_top,
                         const float *pres,   const float* u_wind,
                         const float *v_wind, int num_levs);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the non-pressure weighted mean wind over the given sharp::PressureLayer 
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
WindComponents mean_wind_npw(const PressureLayer& layer, 
                             const float* pres,   const float* u_wind,
                             const float* v_wind, int num_levs);


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the non-pressure weighted mean wind over the given bottom and top pressure levels 
 *
 * Computes the non-pressure weighted mean wind over the given bottom and top 
 * pressure levels, and arrays of pressure and wind components with a length 
 * of num_levs. This is a wrapper around the sharp::PressureLayer implementation.
 *
 * \param pressure_bot  (hPa)
 * \param pressure_top  (hPa)
 * \param pres          (hPa)
 * \param u_wind        (kts or m/s)
 * \param v_wind        (kts or m/s)
 * \param num_levs      (length of arrays)
 * \return sharp::WindComponents    {mean_u, mean_v}
 */
WindComponents mean_wind_npw(float pressure_bot,  float pressure_top,
                             const float *pres,   const float* u_wind,
                             const float *v_wind, int num_levs);


}
