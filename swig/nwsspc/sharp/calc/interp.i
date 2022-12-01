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
        (float* vert_array, int NZ1),
        (float* data_array, int NZ2)
    }


/**
 *  Wrap the interp_height function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (interp_height) _interp_height;
%ignore interp_height;

%exception _interp_height {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

float _interp_height(float height_val, 
                     float* vert_array, int NZ1, 
                     float* data_array, int NZ2) {
    if (NZ1 != NZ2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            NZ1, NZ2
        );
        return 0.0;
    }

    return sharp::interp_height(height_val, vert_array, data_array, NZ1);
}

%}


/**
 *  Wrap the interp_pressure function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (interp_pressure) _interp_pressure;
%ignore interp_pressure;

%exception _interp_pressure {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

float _interp_pressure(float pres_val, 
                       float* vert_array, int NZ1, 
                       float* data_array, int NZ2) {
    if (NZ1 != NZ2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            NZ1, NZ2
        );
        return 0.0;
    }

    return sharp::interp_pressure(pres_val, vert_array, data_array, NZ1);
}

%}

%include "../include/SHARPlib/interp.h"
