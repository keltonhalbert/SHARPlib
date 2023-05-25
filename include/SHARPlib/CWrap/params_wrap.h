/**
 * \file
 * \brief Header defining the C wrapper interface to derived parameters 
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

#ifndef __SHARP_PARAMS_WRAP_H__
#define __SHARP_PARAMS_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <SHARPlib/CWrap/profile_wrap.h>
#include <SHARPlib/CWrap/layer_wrap.h>
#include <SHARPlib/CWrap/winds_wrap.h>
#include <SHARPlib/CWrap/parcel_wrap.h>

void sharp_effective_inflow_layer(sharp_Profile_t* prof,
                                  sharp_PressureLayer_t* elyr,
                                  float cape_thresh, float cinh_thresh);

void sharp_storm_motion_bunkers_np(sharp_Profile_t* prof,
                                   sharp_HeightLayer_t* mn_wind_lyr_agl,
                                   sharp_HeightLayer_t* wind_shr_lyr_agl,
                                   sharp_WindComponents_t* storm_motion,
                                   int leftMover, int pressureWeighted);

void sharp_storm_motion_bunkers(sharp_Profile_t* prof,
                                sharp_WindComponents_t* storm_motion,
                                int leftMover);

float sharp_entrainment_cape(sharp_Profile_t* prof, sharp_Parcel_t* pcl);

float sharp_energy_helicity_index(float cape, float helicity);

float sharp_supercell_composite_parameter(float mu_cape, float eff_srh,
                                          float eff_shear);

float sharp_significant_tornado_parameter(sharp_Profile_t* prof,
                                          sharp_Parcel_t* pcl, float helicity,
                                          float bulk_shear);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __SHARP_PARAMS_WRAP_H__
