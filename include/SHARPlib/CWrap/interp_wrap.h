#ifndef __SHARP_INTERP_WRAP_H__
#define __SHARP_INTERP_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

float sharp_interp_height(float height_val, const float* height_arr,
                          const float* data_arr, int NZ);

float sharp_interp_pressure(float pressure_val, const float* pressure_arr,
                            const float* data_arr, int NZ);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // __SHARP_INTERP_WRAP_H__
