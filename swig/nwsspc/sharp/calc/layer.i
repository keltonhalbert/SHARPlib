/* File: layer.i */
%module layer 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/layer.h>
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
%apply (float* IN_ARRAY1, int DIM1) { (const float pressure[], const int N) };
%apply (float* IN_ARRAY1, int DIM1) { (const float height[], const int N) };

%apply (float* IN_ARRAY1, int DIM1) {
        (const float* pressure, const int NZ1),
        (const float* height, const int NZ2)
};

%rename (pressure_layer_to_height) _pres_lyr_to_hght;
%rename (height_layer_to_pressure) _hght_lyr_to_pres;
%ignore pressure_layer_to_height;
%ignore height_layer_to_pressure;

%exception _pres_lyr_to_hght {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _hght_lyr_to_pres {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

sharp::HeightLayer _pres_lyr_to_hght(sharp::PressureLayer layer,
                            const float* pressure, const int NZ1,
                            const float* height, const int NZ2,
                            bool toAGL) {
    if ((NZ1 != NZ2)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            NZ1, NZ2
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    return sharp::pressure_layer_to_height(layer, pressure, height, NZ1, toAGL);
}

sharp::PressureLayer _hght_lyr_to_pres(sharp::HeightLayer layer,
                            const float* pressure, const int NZ1,
                            const float* height, const int NZ2,
                            bool isAGL) {
    if ((NZ1 != NZ2)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            NZ1, NZ2
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    return sharp::height_layer_to_pressure(layer, pressure, height, NZ1, isAGL);
}

%}

/** 
 * These functions are unecessary,
 * as numpy supports these natively
 */ 
%ignore layer_max;
%ignore layer_min;
%ignore layer_mean;


%include "../include/SHARPlib/layer.h"
