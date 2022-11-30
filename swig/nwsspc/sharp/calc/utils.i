/* File: utils.i */
%module utils 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/utils.h"
%}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

/* Set up an argument typemap for our arrays */
%apply (float* IN_ARRAY1, int DIM1) { (const float* pressure, int num_levs) };
%apply (float* IN_ARRAY1, int DIM1) { (const float* height, int num_levs) };


/** 
 * These functions are unecessary,
 * as numpy supports these natively
 */ 
%ignore max_value;
%ignore min_value;
%ignore mean_value;


%include "../include/utils.h"
