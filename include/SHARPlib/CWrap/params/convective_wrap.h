/**
 * \file
 * \brief Header defining the C wrapper interface to derived parameters
 * \author
 *   Kelton Halbert                     \n
 *   Email: kelton.halbert@noaa.gov     \n
 * \date    2023-05-19
 *
 * Written for the NWS Storm Prediction Center,
 * based on NSHARP routines written by John Hart
 * Rich Thompson.
 */

#ifndef SHARP_PARAMS_WRAP_H
#define SHARP_PARAMS_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <SHARPlib/CWrap/layer_wrap.h>
#include <SHARPlib/CWrap/parcel_wrap.h>
#include <SHARPlib/CWrap/profile_wrap.h>
#include <SHARPlib/CWrap/winds_wrap.h>

void sharp_effective_inflow_layer(const float* pressure, const float* height,
                                  const float* temperature,
                                  const float* dewpoint,
                                  const float* virtemp_arr,
                                  const float* buoy_arr, const int N,
                                  sharp_PressureLayer_t* elyr,
                                  float cape_thresh, float cinh_thresh);

void sharp_storm_motion_bunkers_np(const float* pressure, const float* height,
                                   const float* u_wind, const float* v_wind,
                                   const int N,
                                   sharp_HeightLayer_t* mn_wind_lyr_agl,
                                   sharp_HeightLayer_t* wind_shr_lyr_agl,
                                   sharp_WindComponents_t* storm_motion,
                                   int leftMover, int pressureWeighted);

void sharp_storm_motion_bunkers(
    const float* pressure, const float* height, const float* u_wind,
    const float* v_wind, const int N, sharp_PressureLayer_t* eff_infl_lyr,
    sharp_Parcel_t* pcl, sharp_WindComponents_t* storm_motion, int leftMover);

float sharp_entrainment_cape(const float* pressure, const float* height,
                             const float* temperature, const float* mse_arr,
                             const float* u_wind, const float* v_wind,
                             const int N, sharp_Parcel_t* pcl);

float sharp_energy_helicity_index(float cape, float helicity);

float sharp_supercell_composite_parameter(float mu_cape, float eff_srh,
                                          float eff_shear);

float sharp_significant_tornado_parameter(sharp_Parcel_t* pcl,
                                          float lcl_hght_agl, float helicity,
                                          float bulk_shear);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // SHARP_PARAMS_WRAP_H
