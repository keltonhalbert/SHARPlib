/* File: thermo.i */
%module thermo 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/thermo.h>
%}

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

/**
 *  Wrap the lapse_rate function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (temperature_at_mixratio) _temperature_at_mixratio;
%rename (virtual_temperature) _virtual_temperature;
%rename (vapor_pressure_ice) _vapor_pressure_ice;
%rename (specific_humidity) _specific_humidity;
%rename (lcl_temperature) _lcl_temperature;
%rename (vapor_pressure) _vapor_pressure;
%rename (mixratio_ice) _mixratio_ice;
%rename (theta_level) _theta_level;
%rename (lapse_rate) _lapse_rate;
%rename (mixratio) _mixratio;
%rename (theta) _theta;

%ignore temperature_at_mixratio;
%ignore virtual_temperature;
%ignore vapor_pressure_ice;
%ignore specific_humidity;
%ignore lcl_temperature;
%ignore vapor_pressure;
%ignore mixratio_ice;
%ignore theta_level;
%ignore lapse_rate;
%ignore mixratio;
%ignore theta;

%exception _mixratio {
	$action
	if (PyErr_Occurred()) SWIG_fail;
}

%exception _lapse_rate {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

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

%}

%include "../include/SHARPlib/thermo.h"
