/**
 * \file
 * \brief Header defining the C wrapper interface to thermodynamic routines
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

#ifndef __SHARP_THERMO_WRAP_H__
#define __SHARP_THERMO_WRAP_H__


#ifdef __cplusplus
extern "C" {
#endif

#include <SHARPlib/CWrap/layer_wrap.h>

float sharp_wobf(float temperature);

float sharp_vapor_pressure(float pressure, float temperature);

float sharp_vapor_pressure_ice(float pressure, float temperature);

float sharp_lcl_temperature(float temperature, float dewpoint);

float sharp_temperature_at_mixratio(float wv_mixratio, float pressure);

float sharp_theta_level(float potential_temperature, float temperature);

float sharp_theta(float pressure, float temperature, float ref_pressure);

float sharp_mixratio(float pressure, float temperature);

float sharp_mixratio_ice(float pressure, float temperature);

float sharp_specific_humidity(float q);

float sharp_virtual_temperature(float temperature, float qv, float ql,
                                float qi);

float sharp_saturated_lift(float pressure, float theta_sat);

float sharp_wetlift(float pressure, float temperature, float lifted_pressure);

float sharp_moist_adiabat_cm1(float pressure, float temperature,
                              float new_pressure, float* qv_total, float* qv,
                              float* ql, float* qi, const float pres_incr,
                              const float converge, const int ma_type);

void sharp_drylift(float pressure, float temperature, float dewpoint,
                   float* pressure_at_lcl, float* temperature_at_lcl);

float sharp_lifted(float pressure, float temperature, float dewpoint,
                   float lifted_pressure);

float sharp_wetbulb(float pressure, float temperature, float dewpoint);

float sharp_theta_wetbulb(float pressure, float temperature, float dewpoint);

float sharp_thetae(float pressure, float temperature, float dewpoint);

float sharp_HeightLayer_lapse_rate(sharp_HeightLayer_t* layer_agl,
                                   const float* height,
                                   const float* temperature, const int N);

float sharp_PressureLayer_lapse_rate(sharp_PressureLayer_t* layer,
                                     const float* pressure, const float* height,
                                     const float* temperature, const int N);

float sharp_buoyancy(float pcl_temperature, float env_temperature);

float sharp_moist_static_energy(float height_agl, float temperature,
                                float specific_humidity);

float sharp_buoyancy_dilution_potential(float temperature, float mse_bar,
                                        float saturation_mse);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
