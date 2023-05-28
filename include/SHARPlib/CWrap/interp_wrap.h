/**
 * \file
 * \brief Header defining the C wrapper interface to interpolation routines
 * \author
 *   Kelton Halbert                     \n
 *   Email: kelton.halbert@noaa.gov     \n
 *   License: Apache 2.0                \n
 * \date    2023-05-19
 *
 * Written for the NWS Storm Prediction Center,
 * based on NSHARP routines written by John Hart 
 * Rich Thompson.
 */

#ifndef __SHARP_INTERP_WRAP_H__
#define __SHARP_INTERP_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

float sharp_interp_height(float height_val, const float* height_arr,
                          const float* data_arr, const int NZ);

float sharp_interp_pressure(float pressure_val, const float* pressure_arr,
                            const float* data_arr, const int NZ);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __SHARP_INTERP_WRAP_H__
