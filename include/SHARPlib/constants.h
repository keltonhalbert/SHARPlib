/**
 * \file
 * \brief Meteorological constants useful for computations with vertical <!--
 * --> atmospheric sounding profiles 
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

#ifndef __SHARP_CONSTANTS_H__
#define __SHARP_CONSTANTS_H__

namespace sharp {

// conversion between celsius and kelvin
static constexpr float ZEROCNK = 273.15f;
// missing array values
static constexpr float MISSING = -9999.0f;  
// Gravitational constant in m/s^2
static constexpr float GRAVITY = 9.80665f;  
// Floating point tolerance for iterative convergence functions
static constexpr float TOL = 1e-10;
// You know... Pi... radians...
static constexpr float PI = 3.14159265;  
// Dry air gas constant, J/(kg*K)
static constexpr float RDGAS = 287.052874;  
// Water vapor gas constant, J/(kg*K)
static constexpr float RVGAS = 461.52f;  
static constexpr float EPSILON = RDGAS / RVGAS;
// specific heat capacity of dry air (J/kg)
static constexpr float CP_DRYAIR = 1005.7f;
// specific heat capacity of water vapor (J/kg) at const pressure
static constexpr float CP_VAPOR = 1870.0f;
// specific heat capacity of water (J/kg)
static constexpr float CP_WATER = 4190.f;
// specific heat capacity of ice (J/kg)
static constexpr float CP_ICE = 2106.0f;
// Latent heat of vaporization of water
// (J/kg) at triple point temperature
static constexpr float LV = 2.501e6f;
// Latent heat of sublimation of water
// ice (J/kg) at triple point temperature
static constexpr float LS = 2.834e6f;

// static constexpr float ROCP = RDGAS / CP_DRYAIR;
static constexpr float ROCP = 0.28571428;
static constexpr float GAMMA_D = GRAVITY / CP_DRYAIR;
// Prandtl number
static constexpr float PRANDTL = 1. / 3.;
//  Von Karman constant (k^2)
static constexpr float VKSQ = 0.18;
// Reference vapor pressure of water vapor at triple point
// temperature (Pa)
static constexpr float VAPPRES_REF = 611.65;

} // end namespace sharp

#endif // __SHARP_CONSTANTS_H__
