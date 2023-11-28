/**
 * \file
 * \brief Routines used to compute winter weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2023-11-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef __SHARP_PARAMS_WINTER_H__
#define __SHARP_PARAMS_WINTER_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#include <algorithm>
#include <cmath>

namespace sharp {

/**
 * \author Nathan Dahl - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute positive and negative temperature areas as related to winter weather forecasting
 *
 * The function integrates regions above and below 0 C in the layer between designated pressure lvl
 * and the surface, but only returns a result if a positive area is found overlying a negative area
 *
 * \param   start	   	Top of search layer for positive/negative areas (hPa)
 * \param   pres_arr		Array containing pressure levels of profile (hPa)
 * \param   height_arr		Array containing height levels of profile (m)
 * \param   temperature_arr	Array containing environmental temperature profile (K)
 * \param   N			Length of arrays
 * \param   *pos 		positive area (J/kg)
 * \param   *neg		negative area (J/kg)
 */
void posneg_temperature(float start, const float pres_arr[], const float height_arr[],
              	const float temperature_arr[], const int N, float* pos, float* neg) noexcept;

/**
 * \author Nathan Dahl - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Evaluate initial precip phase as related to winter weathre forecasting
 *
 * The function uses the vertical profile of temperature, relative humidity, and vertical velocity
 * to determine the water phase in the precipitation generation layer (if any)
 *
 * \param   pres_arr            Array containing pressure levels of profile (hPa)
 * \param   height_arr          Array containing height levels of profile (m)
 * \param   temperature_arr     Array containing environmental temperature profile (K)
 * \param   dewpoint_arr     	Array containing environmental dewpoint profile (K)
 * \param   vvel_arr     	Array containing environmental vertical velocity profile (m/s)
 * \param   N                   Length of arrays
 * \param   *plevel             Pressure level of precip origin (Pa)
 * \param   *phase              Integer designating phase (0=liquid,1=freezing rain/mix,3=snow)
 *
 * \return  *init_phase		Precip phase label string
 */
char *init_phase(const float pres_arr[], const float height_arr[], const float temperature_arr[],
                const float dewpoint_arr[], const float vvel_arr[], const int N, float *plevel, short *phase) noexcept;

}  // end namespace sharp

#endif // __SHARP_THERMP_H__
