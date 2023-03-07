/* File: utils.i */
%module utils 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/utils.h>
%}

%include exception.i

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

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


%include "../include/SHARPlib/utils.h"
