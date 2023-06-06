/* File: winds.i */
%module winds 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/winds.h>
%}
%import "../include/SHARPlib/layer.h"

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

/* Set up argument typemaps for our arrays */

%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (float** out_arr, int* NOUT)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pres[], const int N1),
    (const float u_wind[], const int N2),
    (const float v_wind[], const int N3)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float height[], const int N1),
    (const float u_wind[], const int N2),
    (const float v_wind[], const int N3)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float wind_speed[], const int N1),
    (const float wind_direction[], const int N2)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float u_comp[], const int N1),
    (const float v_comp[], const int N2)
}

%rename (vector_magnitude) _vector_magnitude;
%rename (vector_angle) _vector_angle;
%rename (u_component) _u_component;
%rename (v_component) _v_component;
%rename (wind_shear) _wind_shear;
%rename (mean_wind) _mean_wind;
%rename (helicity) _helicity;

%ignore vector_magnitude;
%ignore vector_angle;
%ignore u_component;
%ignore v_component;
%ignore wind_shear;
%ignore mean_wind;
%ignore helicity;

%exception _mean_wind {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _mean_wind_npw {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _wind_shear {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _helicity {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

float _u_component(float wind_speed, float wind_direction) {
    return sharp::u_component(wind_speed, wind_direction);
}

void _u_component(const float wind_speed[], const int N1,
                  const float wind_direction[], const int N2,
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
        temp[k] = sharp::u_component(wind_speed[k], wind_direction[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _v_component(float wind_speed, float wind_direction) {
    return sharp::v_component(wind_speed, wind_direction);
}

void _v_component(const float wind_speed[], const int N1,
                  const float wind_direction[], const int N2,
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
        temp[k] = sharp::v_component(wind_speed[k], wind_direction[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _vector_angle(float u_comp, float v_comp) {
    return sharp::vector_angle(u_comp, v_comp);
}

void _vector_angle(const float u_comp[], const int N1,
                   const float v_comp[], const int N2,
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
        temp[k] = sharp::vector_angle(u_comp[k], v_comp[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}

float _vector_magnitude(float u_comp, float v_comp) {
    return sharp::vector_magnitude(u_comp, v_comp);
}

void _vector_magnitude(const float u_comp[], const int N1,
                       const float v_comp[], const int N2,
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
        temp[k] = sharp::vector_magnitude(u_comp[k], v_comp[k]);
    }

    *out_arr = temp;
    *NOUT = N1;
    return;
}


sharp::WindComponents _mean_wind(sharp::PressureLayer layer, 
                 const float pres[], const int N1,
                 const float u_wind[], const int N2,
                 const float v_wind[], const int N3,
                 bool weighted) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::WindComponents();
    }

    return sharp::mean_wind(layer, pres, u_wind, v_wind, N1, weighted);
}

sharp::WindComponents _wind_shear(sharp::PressureLayer layer, 
                 const float pres[], const int N1,
                 const float u_wind[], const int N2,
                 const float v_wind[], const int N3) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::WindComponents();
    }
    
    return sharp::wind_shear(layer, pres, u_wind, v_wind, N1);
}

sharp::WindComponents _wind_shear(sharp::HeightLayer layer, 
                 const float height[], const int N1,
                 const float u_wind[], const int N2,
                 const float v_wind[], const int N3) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::WindComponents();
    }
    
    return sharp::wind_shear(layer, height, u_wind, v_wind, N1);
}

float _helicity(sharp::HeightLayer layer_agl, sharp::WindComponents storm_motion,
               const float height[], const int N1,
               const float u_wind[], const int N2,
               const float v_wind[], const int N3) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return 0.0;
    }
    
    return sharp::helicity(layer_agl, storm_motion, height, u_wind, v_wind, N1); 
}

float _helicity(sharp::PressureLayer layer, 
                sharp::WindComponents storm_motion,
                const float pres[],   const int N1,
                const float u_wind[], const int N2,
                const float v_wind[], const int N3) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return 0.0;
    }
    
    return sharp::helicity(layer, storm_motion, pres, u_wind, v_wind, N1); 
}

%}


%include "../include/SHARPlib/winds.h"
