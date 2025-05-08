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

#ifndef SHARP_CONSTANTS_H
#define SHARP_CONSTANTS_H

namespace sharp {

// Molar gas constant (J / mol / K)
static constexpr float RGAS = 8.314462618f;
// Dry air molecular weight kg / mol
static constexpr float MDRY = 28.96546e-3f;
// Water molecular weight kg / mol
static constexpr float MWATER = 18.015268e-3f;
static constexpr float DRYAIR_SPEC_HEAT_RATIO = 1.4f;  // unitless

// conversion between celsius and kelvin
static constexpr float ZEROCNK = 273.15f;
// hPa to Pa conversion
static constexpr float HPA_TO_PA = 100.0f;
static constexpr float PA_TO_HPA = 0.01f;
// Potential Temperature Reference Pressure (Pa)
static constexpr float THETA_REF_PRESSURE = 100000.0f;
// missing array values
static constexpr float MISSING = -9999.0f;
// Gravitational constant in m/s^2
static constexpr float GRAVITY = 9.80665f;
// Floating point tolerance for iterative convergence functions
static constexpr float TOL = 1e-10f;
// You know... Pi... radians...
static constexpr float PI = 3.14159265f;
// Dry air gas constant, J/(kg*K)
/*static constexpr float RDGAS = 287.052874;*/
static constexpr float RDGAS = RGAS / MDRY;
// Water vapor gas constant, J/(kg*K)
/*static constexpr float RVGAS = 461.52f;*/
static constexpr float RVGAS = RGAS / MWATER;
static constexpr float EPSILON = RDGAS / RVGAS;
// specific heat capacity of dry air (J/kg)
/*static constexpr float CP_DRYAIR = 1005.7f;*/
static constexpr float CP_DRYAIR =
    DRYAIR_SPEC_HEAT_RATIO * RDGAS / (DRYAIR_SPEC_HEAT_RATIO - 1.0f);
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
static constexpr float PRANDTL = 1.f / 3.f;
//  Von Karman constant (k^2)
static constexpr float VKSQ = 0.18f;
// Reference vapor pressure of water vapor at triple point
// temperature (Pa)
static constexpr float VAPPRES_REF = 611.65f;
// Density of liquid water - kg / m^3
static constexpr float RHO_LWAT = 999.97495;

}  // end namespace sharp

#endif  // SHARP_CONSTANTS_H
