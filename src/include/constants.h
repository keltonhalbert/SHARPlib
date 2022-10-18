/**
 * \file
 * \brief Meteorological constants useful for computations with vertical atmospheric sounding profiles 
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

	static const float ZEROCNK = 273.15f; // conversion between celsius and kelvin
	static const float MISSING = -9999.0f; // missing array values
	static const float GRAVITY = 9.80665f; // Gravitational constant in m/s^2
	static const float TOL = 1e-10; // Floating point tolerance for iterative convergence functions
	static const float PI = 3.14159265; // You know... Pi...
	static const float RDGAS = 287.052874; // Dry air gas constant, J/(kg*K)
	static const float RVGAS = 461.52f; // Water vapor gas constant, J/(kg*K)
	static const float CP_DRYAIR = 1005.7f; // specific heat capacity of dry air J/kg
	static const float CP_WATER = 4184.f; // specific heat capacity of water J/kg
	static const float LV = 2.5104e6f; // Latent heat of vaporization of water (J/kg)
	static const float ROCP = RDGAS / CP_DRYAIR;
}
