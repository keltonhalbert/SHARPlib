/**
 * \file
 * \brief Header defining the C wrapper interface to interpolation routines
 * \author
 *   Kelton Halbert                     \n
 *   Email: kelton.halbert@noaa.gov     \n
 * \date    2023-05-19
 *
 * Written for the NWS Storm Prediction Center,
 * based on NSHARP routines written by John Hart
 * Rich Thompson.
 */

#ifndef SHARP_INTERP_WRAP_H
#define SHARP_INTERP_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

float sharp_interp_height(float height_val, const float* height_arr,
                          const float* data_arr, const int N);

float sharp_interp_pressure(float pressure_val, const float* pressure_arr,
                            const float* data_arr, const int N);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // SHARP_INTERP_WRAP_H
