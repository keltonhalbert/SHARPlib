/* File: thermo.i */
%module thermo 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/thermo.h>
%}

%include exception.i

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
        return NULL;
    }
}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

%import "../include/SHARPlib/layer.h"

%apply float &OUTPUT { float& pressure_at_lcl };
%apply float &OUTPUT { float& temperature_at_lcl };

/* Set up an argument typemap for our arrays */
%apply (float* IN_ARRAY1, int DIM1) {
    (const float hght_array[], const int N1),
    (const float tmpc_array[], const int N2)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pres_array[], const int N1),
    (const float hght_array[], const int N2),
    (const float tmpc_array[], const int N3)
}

%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (float** out_arr, int* NOUT)
}

%apply (float* IN_ARRAY1, int DIM1) {
	(const float pressure[], const int N1),
	(const float temperature[], const int N2)
}
%apply (float* IN_ARRAY1, int DIM1) {
	(const float temperature[], const int N1),
	(const float dewpoint[], const int N2)
}
%apply (float* IN_ARRAY1, int DIM1) {
	(const float wv_mixratio[], const int N1),
	(const float pressure[], const int N2)
}
%apply (float* IN_ARRAY1, int DIM1) {
	(const float potential_temperature[], const int N1),
	(const float temperature[], const int N2)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float qv[], const int N1)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float temperature[], const int N1),
    (const float qv[], const int N2)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float temperature[], const int N1),
    (const float qv[], const int N2),
    (const float ql[], const int N3),
    (const float qi[], const int N4)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const int N1),
    (const float temperature[], const int N2),
    (const float dewpoint[], const int N3)
}
%apply (float* IN_ARRAY1, int DIM1) {
	(const float pcl_temperature[], const int N1),
	(const float env_temperature[], const int N2)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float height_agl[], const int N1),
    (const float temperature[], const int N2),
    (const float spec_hum[], const int N3)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float temperature[], const int N1),
    (const float mse_bar[], const int N2),
    (const float sat_mse[], const int N3)
}

/**
 *  Wrap the lapse_rate function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (buoyancy_dilution_potential) _buoyancy_dilution_potential;
%rename (temperature_at_mixratio) _temperature_at_mixratio;
%rename (moist_static_energy) _moist_static_energy;
%rename (virtual_temperature) _virtual_temperature;
%rename (vapor_pressure_ice) _vapor_pressure_ice;
%rename (specific_humidity) _specific_humidity;
%rename (lcl_temperature) _lcl_temperature;
%rename (vapor_pressure) _vapor_pressure;
%rename (lapse_rate_max) _lapse_rate_max;
%rename (theta_wetbulb) _theta_wetbulb;
%rename (mixratio_ice) _mixratio_ice;
%rename (theta_level) _theta_level;
%rename (lapse_rate) _lapse_rate;
%rename (buoyancy) _buoyancy;
%rename (mixratio) _mixratio;
%rename (wetbulb) _wetbulb;
%rename (thetae) _thetae;
%rename (theta) _theta;

%ignore buoyancy_dilution_potential;
%ignore temperature_at_mixratio;
%ignore moist_static_energy;
%ignore virtual_temperature;
%ignore vapor_pressure_ice;
%ignore specific_humidity;
%ignore lcl_temperature;
%ignore vapor_pressure;
%ignore lapse_rate_max;
%ignore theta_wetbulb;
%ignore mixratio_ice;
%ignore theta_level;
%ignore lapse_rate;
%ignore buoyancy;
%ignore mixratio;
%ignore wetbulb;
%ignore thetae;
%ignore theta;

%inline %{

float _vapor_pressure(float pressure, float temperature) {
    return sharp::vapor_pressure(pressure, temperature);
}

void _vapor_pressure(const float pressure[], const int N1,
                     const float temperature[], const int N2,
                     float** out_arr, int* NOUT ) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::vapor_pressure(pressure[k], temperature[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

inline float _vapor_pressure_ice(float pressure, float temperature) {
    return sharp::vapor_pressure_ice(pressure, temperature);
}

void _vapor_pressure_ice(const float pressure[], const int N1,
                         const float temperature[], const int N2,
                         float** out_arr, int* NOUT ) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::vapor_pressure_ice(pressure[k], temperature[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _lcl_temperature(float temperature, float dewpoint) {
    return sharp::lcl_temperature(temperature, dewpoint);
}

void _lcl_temperature(const float temperature[], const int N1,
                      const float dewpoint[], const int N2,
                      float** out_arr, int* NOUT ) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::lcl_temperature(temperature[k], dewpoint[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _temperature_at_mixratio(float wv_mixratio, float pressure) {
    return sharp::temperature_at_mixratio(wv_mixratio, pressure);
}

void _temperature_at_mixratio(const float wv_mixratio[], const int N1,
                              const float pressure[], const int N2,
                              float** out_arr, int* NOUT ) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::temperature_at_mixratio(wv_mixratio[k], pressure[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _theta_level(float potential_temperature, float temperature) {
    return sharp::theta_level(potential_temperature, temperature);
}

void _theta_level(const float potential_temperature[], const int N1,
                              const float temperature[], const int N2,
                              float** out_arr, int* NOUT ) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::theta_level(potential_temperature[k], temperature[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _theta(float pressure, float temperature,
             const float ref_pressure = sharp::THETA_REF_PRESSURE) {
    return sharp::theta(pressure, temperature, ref_pressure);
}

void _theta(const float pressure[], const int N1,
            const float temperature[], const int N2,
            float** out_arr, int* NOUT,
            const float ref_pressure=sharp::THETA_REF_PRESSURE) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::theta(pressure[k], temperature[k], ref_pressure);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _mixratio(float pressure, float temperature) {
	return sharp::mixratio(pressure, temperature);
}

void _mixratio(const float pressure[], const int N1, 
               const float temperature[], const int N2, 
               float** out_arr, int* NOUT) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::mixratio(pressure[k], temperature[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}


float _mixratio_ice(float pressure, float temperature) {
	return sharp::mixratio_ice(pressure, temperature);
}

void _mixratio_ice(const float pressure[], const int N1, 
                   const float temperature[], const int N2, 
                   float** out_arr, int* NOUT) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::mixratio_ice(pressure[k], temperature[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _specific_humidity(float qv) {
    return sharp::specific_humidity(qv);
}

void _specific_humidity(const float qv[], const int N1, 
                        float** out_arr, int* NOUT) {

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::specific_humidity(qv[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _virtual_temperature(float temperature, float qv, float ql = 0.0f,
                          float qi = 0.0f) {
    return sharp::virtual_temperature(temperature, qv, ql, qi);
}

void _virtual_temperature(const float temperature[], const int N1, 
                          const float qv[], const int N2, 
                          float** out_arr, int* NOUT) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::virtual_temperature(temperature[k], qv[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

void _virtual_temperature(const float temperature[], const int N1, 
                          const float qv[], const int N2, 
                          const float ql[], const int N3,
                          const float qi[], const int N4,
                          float** out_arr, int* NOUT) {
    if ((N1 != N2) || (N1 != N3) || (N1 != N4)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d, %d)",
            N1, N2, N3, N4
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] =
            sharp::virtual_temperature(temperature[k], qv[k], ql[k], qi[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _wetbulb(float pressure, float temperature, float dewpoint) {
    return sharp::wetbulb(pressure, temperature, dewpoint);
}

void _wetbulb(const float pressure[], const int N1,
              const float temperature[], const int N2,
              const float dewpoint[], const int N3,
              float** out_arr, int* NOUT) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] =
            sharp::wetbulb(pressure[k], temperature[k], dewpoint[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _theta_wetbulb(float pressure, float temperature, float dewpoint) {
    return sharp::theta_wetbulb(pressure, temperature, dewpoint);
}

void _theta_wetbulb(const float pressure[], const int N1,
              const float temperature[], const int N2,
              const float dewpoint[], const int N3,
              float** out_arr, int* NOUT) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] =
            sharp::theta_wetbulb(pressure[k], temperature[k], dewpoint[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _thetae(float pressure, float temperature, float dewpoint) {
    return sharp::thetae(pressure, temperature, dewpoint);
}

void _thetae(const float pressure[], const int N1,
             const float temperature[], const int N2,
             const float dewpoint[], const int N3,
             float** out_arr, int* NOUT) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] =
            sharp::thetae(pressure[k], temperature[k], dewpoint[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _lapse_rate(sharp::HeightLayer layer, 
                  const float hght_array[], const int N1, 
                  const float tmpc_array[], const int N2) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::lapse_rate(layer, hght_array, tmpc_array, N1);
}

float _lapse_rate(sharp::PressureLayer layer, 
                  const float pres_array[], const int N1, 
                  const float hght_array[], const int N2,
                  const float tmpc_array[], const int N3) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::MISSING;
    }

    return sharp::lapse_rate(layer, pres_array, hght_array, tmpc_array, N1);
}

float _lapse_rate_max(sharp::HeightLayer layer_agl, const float depth, 
                      const float hght_array[], const int N1, 
                      const float tmpc_array[], const int N2,
                      sharp::HeightLayer* max_lyr=nullptr) {
    if ((N1 != N2)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::lapse_rate_max(layer_agl, depth, hght_array, tmpc_array, N1, max_lyr);
}

float _lapse_rate_max(sharp::PressureLayer layer, const float depth,
                      const float pres_array[], const int N1,
                      const float hght_array[], const int N2, 
                      const float tmpc_array[], const int N3,
                      sharp::PressureLayer* max_lyr=nullptr) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::MISSING;
    }

    return sharp::lapse_rate_max(layer, depth, pres_array, hght_array, tmpc_array, N1, max_lyr);
}

float _buoyancy(float pcl_temperature, float env_temperature) {
    return sharp::buoyancy(pcl_temperature, env_temperature);
}

void _buoyancy(const float pcl_temperature[], const int N1,
               const float env_temperature[], const int N2,
               float** out_arr, int* NOUT) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] =
            sharp::buoyancy(pcl_temperature[k], env_temperature[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _moist_static_energy(float height_agl, float temperature, float spec_hum) {
    return sharp::moist_static_energy(height_agl, temperature, spec_hum);
}

void _moist_static_energy(const float height_agl[], const int N1,
                          const float temperature[], const int N2,
                          const float spec_hum[], const int N3,
                          float** out_arr, int* NOUT) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::moist_static_energy(height_agl[k], temperature[k],
                                             spec_hum[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _buoyancy_dilution_potential(float temperature, float mse_bar,
                                   float sat_mse) {
    return sharp::buoyancy_dilution_potential(temperature, mse_bar, sat_mse);
}

void _buoyancy_dilution_potential(const float temperature[], const int N1,
                                  const float mse_bar[], const int N2,
                                  const float sat_mse[], const int N3,
                                  float** out_arr, int* NOUT) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
		return; 
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", 
            N1
        );
        return;
    }

    for (int k = 0; k < N1; ++k) {
        temp[k] = sharp::buoyancy_dilution_potential(temperature[k], mse_bar[k],
                                                     sat_mse[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

%}

%include "../include/SHARPlib/thermo.h"
