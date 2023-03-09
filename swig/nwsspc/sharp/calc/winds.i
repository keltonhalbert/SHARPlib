/* File: winds.i */
%module winds 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/winds.h>
%}
%import "../include/SHARPlib/utils.h"

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

/* Set up an argument typemap for our arrays */
%apply (float* IN_ARRAY1, int DIM1) {
        (const float* pres, int NZ1),
        (const float* u_wind, int NZ2),
        (const float* v_wind, int NZ3)
}

%apply (float* IN_ARRAY1, int DIM1) {
        (const float* height, int NZ1),
        (const float* u_wind, int NZ2),
        (const float* v_wind, int NZ3)
}

%apply (float* IN_ARRAY1, int DIM1) {
        (const float* pres,   int NZ1),
        (const float* height, int NZ2),
        (const float* u_wind, int NZ3),
        (const float* v_wind, int NZ4)
}

/**
 *  Wrap the mean_wind function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (mean_wind) _mean_wind;
%rename (mean_wind_npw) _mean_wind_npw;
%rename (wind_shear) _wind_shear;
%rename (helicity) _helicity;

%ignore mean_wind;
%ignore mean_wind_npw;
%ignore wind_shear;
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

sharp::WindComponents _mean_wind(sharp::PressureLayer layer, 
                 const float* pres, int NZ1,
                 const float* u_wind, int NZ2,
                 const float* v_wind, int NZ3) {
    if ((NZ1 != NZ2) || (NZ1 != NZ3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            NZ1, NZ2, NZ3
        );
        return sharp::WindComponents();
    }

    return sharp::mean_wind(layer, pres, u_wind, v_wind, NZ1);
}

sharp::WindComponents _mean_wind_npw(sharp::PressureLayer layer, 
                 const float* pres, int NZ1,
                 const float* u_wind, int NZ2,
                 const float* v_wind, int NZ3) {
    if ((NZ1 != NZ2) || (NZ1 != NZ3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            NZ1, NZ2, NZ3
        );
        return sharp::WindComponents();
    }

    return sharp::mean_wind_npw(layer, pres, u_wind, v_wind, NZ1);
}

sharp::WindComponents _wind_shear(sharp::PressureLayer layer, 
                 const float* pres, int NZ1,
                 const float* u_wind, int NZ2,
                 const float* v_wind, int NZ3) {
    if ((NZ1 != NZ2) || (NZ1 != NZ3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            NZ1, NZ2, NZ3
        );
        return sharp::WindComponents();
    }
    
    return sharp::wind_shear(layer, pres, u_wind, v_wind, NZ1);
}

float _helicity(sharp::HeightLayer layer_agl, sharp::WindComponents storm_motion,
               const float* height, int NZ1,
               const float* u_wind, int NZ2,
               const float* v_wind, int NZ3) {
    if ((NZ1 != NZ2) || (NZ1 != NZ3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            NZ1, NZ2, NZ3
        );
        return 0.0;
    }
    
    return sharp::helicity(layer_agl, storm_motion, height, u_wind, v_wind, NZ1); 
}

float _helicity(sharp::PressureLayer layer, 
                sharp::WindComponents storm_motion,
                const float* pres,   int NZ1,
                const float* height, int NZ2,
                const float* u_wind, int NZ3,
                const float* v_wind, int NZ4) {
    if ((NZ1 != NZ2) || (NZ1 != NZ3) || (NZ1 != NZ4)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (i%d, %d, %d, %d)",
            NZ1, NZ2, NZ3, NZ4
        );
        return 0.0;
    }
    
    return sharp::helicity(layer, storm_motion, pres, height, u_wind, v_wind, NZ1); 
}

%}


%include "../include/SHARPlib/winds.h"
