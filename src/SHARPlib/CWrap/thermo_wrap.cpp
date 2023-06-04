#include <SHARPlib/CWrap/thermo_wrap.h>
#include <SHARPlib/CWrap/layer_wrap.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/thermo.h>
#include <stdlib.h>

float sharp_wobf(float temperature) {
    return sharp::wobf(temperature);
}

float sharp_vapor_pressure(float pressure, float temperature) {
    return sharp::vapor_pressure(pressure, temperature);
}

float sharp_vapor_pressure_ice(float pressure, float temperature) {
	return sharp::vapor_pressure_ice(pressure, temperature);
}

float sharp_lcl_temperature(float temperature, float dewpoint) {
    return sharp::lcl_temperature(temperature, dewpoint);
}

float sharp_temperature_at_mixratio(float wv_mixratio, float pressure) {
    return sharp::temperature_at_mixratio(wv_mixratio, pressure);
}

float sharp_theta_level(float potential_temperature, float temperature) {
    return sharp::theta_level(potential_temperature, temperature);
}

float sharp_theta(float pressure, float temperature, float ref_pressure) {
    return sharp::theta(pressure, temperature, ref_pressure);
}

float sharp_mixratio(float pressure, float temperature) {
    return sharp::mixratio(pressure, temperature);
}

float sharp_mixratio_ice(float pressure, float temperature) {
	return sharp::mixratio_ice(pressure, temperature);
}

float sharp_specific_humidity(float q) {
	return sharp::specific_humidity(q);
}

float sharp_virtual_temperature(float temperature, float qv, float ql,
                                float qi) {
	return sharp::virtual_temperature(temperature, qv, ql, qi);
}

float sharp_saturated_lift(float pressure, float theta_sat) {
    return sharp::saturated_lift(pressure, theta_sat);
}

float sharp_wetlift(float pressure, float temperature, float lifted_pressure) {
    return sharp::wetlift(pressure, temperature, lifted_pressure);
}

float sharp_moist_adiabat_cm1(float pressure, float temperature,
                              float new_pressure, float* qv_total, float* qv,
                              float* ql, float* qi, const float pres_incr,
                              const float converge, const int ma_type) {
	sharp::adiabat ma = static_cast<sharp::adiabat>(ma_type);
	return sharp::moist_adiabat_cm1(pressure, temperature, new_pressure,
									*qv_total, *qv, *ql, *qi, pres_incr,
									converge, ma);
}

void sharp_drylift(float pressure, float temperature, float dewpoint,
                   float* pressure_at_lcl, float* temperature_at_lcl) {
    if ((pressure_at_lcl == NULL) || (temperature_at_lcl == NULL)) return;
    sharp::drylift(pressure, temperature, dewpoint, *pressure_at_lcl,
                   *temperature_at_lcl);
}

float sharp_lifted(float pressure, float temperature, float dewpoint,
                   float lifted_pressure) {
    return sharp::lifted(pressure, temperature, dewpoint, lifted_pressure);
}

float sharp_wetbulb(float pressure, float temperature, float dewpoint) {
    return sharp::wetbulb(pressure, temperature, dewpoint);
}

float sharp_theta_wetbulb(float pressure, float temperature, float dewpoint) {
    return sharp::theta_wetbulb(pressure, temperature, dewpoint);
}

float sharp_thetae(float pressure, float temperature, float dewpoint) {
    return sharp::thetae(pressure, temperature, dewpoint);
}

float sharp_HeightLayer_lapse_rate(sharp_HeightLayer_t* layer_agl,
                                   const float* height,
                                   const float* temperature, const int N) {
    if (layer_agl == NULL) return sharp::MISSING;
    sharp::HeightLayer* h = static_cast<sharp::HeightLayer*>(layer_agl->obj);
    return sharp::lapse_rate(*h, height, temperature, N);
}

float sharp_PressureLayer_lapse_rate(sharp_PressureLayer_t* layer,
                                     const float* pressure, const float* height,
                                     const float* temperature, const int N) {
    if (layer == NULL) return sharp::MISSING;
    sharp::PressureLayer* p = static_cast<sharp::PressureLayer*>(layer->obj);
    return sharp::lapse_rate(*p, pressure, height, temperature, N);
}

float sharp_buoyancy(float pcl_temperature, float env_temperature) {
    return sharp::buoyancy(pcl_temperature, env_temperature);
}

float sharp_moist_static_energy(float height_agl, float temperature,
                                float specific_humidity) {
    return sharp::moist_static_energy(height_agl, temperature, specific_humidity);
}

float sharp_buoyancy_dilution_potential(float temperature, float mse_bar,
                                        float saturation_mse) {
	return sharp::buoyancy_dilution_potential(temperature, mse_bar,
											  saturation_mse);
}

