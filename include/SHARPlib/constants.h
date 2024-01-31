/**
 * \file
 * \brief Meteorological constants useful for computations with vertical <!--
 * --> atmospheric sounding profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
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
// hPa to Pa conversion
static constexpr float HPA_TO_PA = 100.0;
static constexpr float PA_TO_HPA = 0.01;
// Potential Temperature Reference Pressure (Pa)
static constexpr float THETA_REF_PRESSURE = 100000.0;
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
static constexpr float CP_VAPOR = 1875.0f;
// specific heat capacity of water (J/kg)
static constexpr float CP_LIQUID = 4190.0f;
// specific heat capacity of ice (J/kg)
static constexpr float CP_ICE = 2118.636f;
// Experimentally derived latent heat of vaporization
// of water (J/kg) at triple point temperature
static constexpr float EXP_LV = 2.501e6f;
// Latent heat of sublimation of water
// ice (J/kg) at triple point temperature
static constexpr float EXP_LS = 2.836017e6f;

// Some microphysics stuff I don't quite understand, but I think
// It's meant to help with a more accurate, temperature-dependent
// version of LV and LS to use in the CM1 lifter.
static constexpr float LV1 = EXP_LV + (CP_LIQUID - CP_VAPOR) * ZEROCNK;
static constexpr float LV2 = CP_LIQUID - CP_VAPOR;
static constexpr float LS1 = EXP_LS + (CP_ICE - CP_VAPOR) * ZEROCNK;
static constexpr float LS2 = CP_ICE - CP_VAPOR;

static constexpr float ROCP = RDGAS / CP_DRYAIR;
// static constexpr float ROCP = 0.28571428;
static constexpr float GAMMA_D = GRAVITY / CP_DRYAIR;
// Prandtl number
static constexpr float PRANDTL = 1. / 3.;
//  Von Karman constant (k^2)
static constexpr float VKSQ = 0.18;
// Reference vapor pressure of water vapor at triple point
// temperature (Pa)
static constexpr float VAPPRES_REF = 611.65;

}  // end namespace sharp

#endif  // __SHARP_CONSTANTS_H__
