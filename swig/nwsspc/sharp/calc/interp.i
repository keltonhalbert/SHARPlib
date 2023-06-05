/* File: interp.i */
%module interp 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/interp.h>
%}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

/* Set up an argument typemap for our arrays */
%apply (float* IN_ARRAY1, int DIM1) {
    (const float vert_array[], const int N1),
    (const float data_array[], const int N2)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float vert_array[], const int N2),
    (const float data_array[], const int N3)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float height_vals[], const int N1)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure_vals[], const int N1)
}

%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (float** interp_vals, int* N4)
}

%rename (interp_height) _interp_height;
%ignore interp_height;

%exception _interp_height {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

void _interp_height(const float height_vals[], const int N1,
                    const float vert_array[], const int N2,
                    const float data_array[], const int N3,
                    float** interp_vals, int* N4) {

    if (N2 != N3) {
        PyErr_Format(
            PyExc_ValueError, 
            "Coordinate and data arrays must be same length, got (%d, %d)",
            N2, N3
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", N1
        );
        return;
    }
    *interp_vals = temp; 
    *N4 = N1;

    for (int k = 0; k < N1; ++k) {
        (*interp_vals)[k] =
            sharp::interp_height(height_vals[k], vert_array, data_array, N2);
    }
    return;
}

float _interp_height(float height_val, 
                     const float vert_array[], const int N1, 
                     const float data_array[], const int N2) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::interp_height(height_val, vert_array, data_array, N1);
}

%}


%rename (interp_pressure) _interp_pressure;
%ignore interp_pressure;

%exception _interp_pressure {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

void _interp_pressure(const float pressure_vals[], const int N1,
                      const float vert_array[], const int N2,
                      const float data_array[], const int N3,
                      float** interp_vals, int* N4) {

    if (N2 != N3) {
        PyErr_Format(
            PyExc_ValueError, 
            "Coordinate and data arrays must be same length, got (%d, %d)",
            N2, N3
        );
        return;
    }

    float* temp = (float *)malloc(N1*sizeof(float));
    if (temp == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for output array of size %d.", N1
        );
        return;
    }
    *interp_vals = temp; 
    *N4 = N1;

    for (int k = 0; k < N1; ++k) {
        (*interp_vals)[k] =
            sharp::interp_pressure(pressure_vals[k], vert_array, data_array, N2);
    }
    return;
}

float _interp_pressure(float pres_val, 
                       const float vert_array[], const int N1, 
                       const float data_array[], const int N2) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::interp_pressure(pres_val, vert_array, data_array, N1);
}

%}

%include "../include/SHARPlib/interp.h"
