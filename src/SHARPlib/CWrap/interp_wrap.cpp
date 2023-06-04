#include <SHARPlib/CWrap/interp_wrap.h>
#include <SHARPlib/interp.h>
#include <stdlib.h>

float sharp_interp_height(float height_val, const float* height_arr,
                          const float* data_arr, const int N) {
    return sharp::interp_height(height_val, height_arr, data_arr, N);
}

float sharp_interp_pressure(float pressure_val, const float* pressure_arr,
                            const float* data_arr, const int N) {
    return sharp::interp_pressure(pressure_val, pressure_arr, data_arr, N);
}

