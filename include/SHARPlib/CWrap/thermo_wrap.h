#ifndef __SHARP_THERMO_WRAP_H__
#define __SHARP_THERMO_WRAP_H__


#ifdef __cplusplus
extern "C" {
#endif

#include <SHARPlib/CWrap/layer_wrap.h>

float sharp_wobf(float temperature);

float sharp_vapor_pressure(float temperature);

float sharp_lcl_temperature(float temperature, float dewpoint);

float sharp_temperature_at_mixratio(float wv_mixratio, float pressure);

float sharp_theta_level(float potential_temperature, float temperature);

float sharp_theta(float pressure, float temperature, float ref_pressure);

float sharp_mixratio(float pressure, float temperature);

float sharp_virtual_temperature(float pressure, float temperature,
                                float dewpoint);

float sharp_saturated_lift(float pressure, float theta_sat);

float sharp_wetlift(float pressure, float temperature, float lifted_pressure);

void sharp_drylift(float pressure, float temperature, float dewpoint,
                   float* pressure_at_lcl, float* temperature_at_lcl);

float sharp_lifted(float pressure, float temperature, float dewpoint,
                   float lifted_pressure);

float sharp_wetbulb(float pressure, float temperature, float dewpoint);

float sharp_theta_wetbulb(float pressure, float temperature, float dewpoint);

float sharp_thetae(float pressure, float temperature, float dewpoint);

float sharp_HeightLayer_lapse_rate(sharp_HeightLayer_t* layer_agl,
                                   const float* height,
                                   const float* temperature, int num_levs);

float sharp_PressureLayer_lapse_rate(sharp_PressureLayer_t* layer,
                                     const float* pressure, const float* height,
                                     const float* temperature, int num_levs);

float sharp_buoyancy(float pcl_temperature, float env_temperature);

float sharp_moist_static_energy(float height_agl, float temperature,
                                float specific_humidity);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
