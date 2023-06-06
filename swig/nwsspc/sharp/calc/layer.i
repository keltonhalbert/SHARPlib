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
        (const float pressure[], const int N1),
        (const float height[], const int N2)
};

%apply (float* IN_ARRAY1, int DIM1) {
    (const float coord_arr[], const int N1),
    (const float data_arr[], const int N2)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const int N1),
    (const float data_arr[], const int N2)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float height[], const int N1),
    (const float pressure[], const int N2),
    (const float data_arr[], const int N3)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float var_array[], const int N1),
    (const float coord_array[], const int N2)
}

%apply float* OUTPUT { float* lvl_of_min } 
%apply float* OUTPUT { float* lvl_of_max } 

%rename (integrate_layer_trapz) _integrate_layer_trapz;
%rename (pressure_layer_to_height) _pres_lyr_to_hght;
%rename (height_layer_to_pressure) _hght_lyr_to_pres;
%rename (layer_min) _layer_min;
%rename (layer_max) _layer_max;
%rename (layer_mean) _layer_mean;

%ignore pressure_layer_to_height;
%ignore height_layer_to_pressure;
%ignore integrate_layer_trapz;
%ignore layer_min;
%ignore layer_max;
%ignore layer_mean;

%exception _pres_lyr_to_hght {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _hght_lyr_to_pres {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _integ_layer_trapz {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _layer_min {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _layer_max {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%exception _layer_mean {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

sharp::HeightLayer _pres_lyr_to_hght(sharp::PressureLayer layer,
                            const float pressure[], const int N1,
                            const float height[], const int N2,
                            bool toAGL = false) {
    if ((N1 != N2)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            N1, N2
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    return sharp::pressure_layer_to_height(layer, pressure, height, N1, toAGL);
}

sharp::PressureLayer _hght_lyr_to_pres(sharp::HeightLayer layer,
                            const float pressure[], const int N1,
                            const float height[], const int N2,
                            bool isAGL = false) {
    if ((N1 != N2)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            N1, N2
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    return sharp::height_layer_to_pressure(layer, pressure, height, N1, isAGL);
}

float _layer_min(sharp::HeightLayer layer, 
                 const float coord_arr[], const int N1,
                 const float data_arr[], const int N2,
                 float* lvl_of_min) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::layer_min(layer, coord_arr, data_arr, N1, lvl_of_min);
}

float _layer_min(sharp::PressureLayer layer, 
                 const float coord_arr[], const int N1,
                 const float data_arr[], const int N2,
                 float* lvl_of_min) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::layer_min(layer, coord_arr, data_arr, N1, lvl_of_min);
}

float _layer_max(sharp::HeightLayer layer, 
                 const float coord_arr[], const int N1,
                 const float data_arr[], const int N2,
                 float* lvl_of_max) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::layer_max(layer, coord_arr, data_arr, N1, lvl_of_max);
}

float _layer_max(sharp::PressureLayer layer, 
                 const float coord_arr[], const int N1,
                 const float data_arr[], const int N2,
                 float* lvl_of_max) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d,)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::layer_max(layer, coord_arr, data_arr, N1, lvl_of_max);
}

float _layer_mean(sharp::HeightLayer layer, 
                  const float height[], const int N1,
                  const float pressure[], const int N2,
                  const float data_arr[], const int N3,
                  const bool isAGL = false) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::MISSING;
    }

    return sharp::layer_mean(layer, height, pressure, data_arr, N1, isAGL);
}

float _layer_mean(sharp::PressureLayer layer, 
                  const float pressure[], const int N1,
                  const float data_arr[], const int N2) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::layer_mean(layer, pressure, data_arr, N1);
}


float _integrate_layer_trapz(sharp::HeightLayer layer, 
                            const float var_array[], const int N1,
                            const float coord_array[], const int N2,
                            const int integ_sign = 0,
                            const bool weighted = false) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::integrate_layer_trapz(layer, var_array, coord_array, N1,
                                        integ_sign, weighted);
}

%} /*end inline*/

%include "../include/SHARPlib/layer.h"
