/* File: profile.i */
%module profile 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/SHARPlib.h>

%}


/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

/* Set up an argument typemap for our arrays */
%apply (float *IN_ARRAY1, int DIM1) {
            (float *pres, int NZ1), (float *hght, int NZ2),
            (float *tmpc, int NZ3), (float *dwpc, int NZ4),
            (float *wnd1, int NZ5), (float *wnd2, int NZ6)
        }

%rename (create_profile) _create_profile;
%ignore create_profile;

%exception _create_profile {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

sharp::Profile* _create_profile(float *pres, int NZ1, float *hght, int NZ2,
                                float *tmpc, int NZ3, float *dwpc, int NZ4,
                                float *wnd1, int NZ5, float *wnd2, int NZ6,
                                sharp::Source sounding_type, 
                                bool windComponents) {
    if ( (NZ1 != NZ2) || (NZ1 != NZ3) || (NZ1 != NZ4) ||
         (NZ1 != NZ5) || (NZ1 != NZ6) ) {
        PyErr_Format(PyExc_ValueError, 
                     "Arrays must be same lenght, insead got (%d, %d, %d, %d, %d, %d)",
                     NZ1, NZ2, NZ3, NZ4, NZ5, NZ6);
        return nullptr;
    }

    return sharp::create_profile(pres, hght, tmpc, dwpc, wnd1, wnd2, NZ1,
                                sounding_type, windComponents);
}

%}

%include "../include/SHARPlib/profile.h"
