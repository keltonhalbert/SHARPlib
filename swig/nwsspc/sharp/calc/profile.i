/* File: profile.i */
%module profile 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/profile.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%numpy_typemaps(float, NPY_DOUBLE, int)
%apply (float *IN_ARRAY1, int DIM1) {
            (float *pres, int NZ1), (float *hght, int NZ2),
            (float *tmpc, int NZ3), (float *dwpc, int NZ4),
            (float *wnd1, int NZ5), (float *wnd2, int NZ6)
        }

%import "../include/winds.h"



%rename (create_profile) profile_wrap;
%ignore create_profile;
%exception profile_wrap {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

sharp::Profile* profile_wrap(float *pres, int NZ1, float *hght, int NZ2,
                             float *tmpc, int NZ3, float *dwpc, int NZ4,
                             float *wnd1, int NZ5, float *wnd2, int NZ6,
                             sharp::Source sounding_type, bool windComponents) {
    if ( (NZ1 != NZ2) || (NZ1 != NZ2) || (NZ1 != NZ3) ||
         (NZ1 != NZ4) || (NZ1 != NZ5) || (NZ1 != NZ6) ) {
        PyErr_Format(PyExc_ValueError, "Arrays must be same lenght");
        return nullptr;
    }

    return create_profile(pres, hght, tmpc, dwpc, wnd1, wnd2, NZ1,
                                 sounding_type, windComponents);
}
%}

%include "../include/profile.h"
