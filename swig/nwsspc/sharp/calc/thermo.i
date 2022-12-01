/* File: thermo.i */
%module thermo 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/thermo.h"
%}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

%import "../include/utils.h"


/* Set up an argument typemap for our arrays */
%apply (float* IN_ARRAY1, int DIM1) {
        (float* hght_array, int NZ1),
        (float* tmpc_array, int NZ2)
    }

%apply (float* IN_ARRAY1, int DIM1) {
        (float* pres_array, int NZ1),
        (float* hght_array, int NZ2),
        (float* tmpc_array, int NZ3)
    }

/**
 *  Wrap the lapse_rate function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (lapse_rate) _lapse_rate;
%ignore lapse_rate;

%exception _lapse_rate {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

float _lapse_rate(sharp::HeightLayer layer, 
                  float* hght_array, int NZ1, 
                  float* tmpc_array, int NZ2) {
    if (NZ1 != NZ2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            NZ1, NZ2
        );
        return 0.0;
    }

    return sharp::lapse_rate(layer, hght_array, tmpc_array, NZ1);
}

float _lapse_rate(sharp::PressureLayer layer, 
                  float* pres_array, int NZ1, 
                  float* hght_array, int NZ2,
                  float* tmpc_array, int NZ3) {
    if ((NZ1 != NZ2) || (NZ1 != NZ3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            NZ1, NZ2, NZ3
        );
        return 0.0;
    }

    return sharp::lapse_rate(layer, pres_array, hght_array, tmpc_array, NZ1);
}

%}


%include "../include/thermo.h"
