/* File: params.i */
%module params 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
    #include <SHARPlib/params/convective.h>
%}

%include exception.i
/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const int N1),
    (const float height[], const int N2),
    (const float temperature[], const int N3),
    (const float dewpoint[], const int N4),
    (const float virtemp_arr[], const int N5)
}

%rename (effective_inflow_layer) _effective_inflow_layer;
%ignore effective_inflow_layer;

%exception _effective_inflow_layer {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

sharp::PressureLayer _effective_inflow_layer(
                        const float pressure[], const int N1,
                        const float height[], const int N2,
                        const float temperature[], const int N3,
                        const float dewpoint[], const int N4,
                        const float virtemp_arr[], const int N5,
                        sharp::Parcel* mupcl = NULL,
                        const float cape_thresh = 100.0,
                        const float cinh_thresh = -250.0) {
	if ( (N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5) ) {
        PyErr_Format(PyExc_ValueError, 
            "Arrays must be same lenght, insead got (%d, %d, %d, %d, %d)",
            N1, N2, N3, N4, N5
        );
        return {sharp::MISSING, sharp::MISSING};
	}

    /*Allocate temporary array for holding Buoyancy*/
    float* buoy = (float *)malloc(N1*sizeof(float));
    if (buoy == NULL) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for temporary array of size %d.", 
            N1
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    sharp::PressureLayer eil = sharp::effective_inflow_layer(
        pressure, height, temperature, dewpoint, virtemp_arr, buoy, N1,
        cape_thresh, cinh_thresh, mupcl);
    free(buoy);
    return eil;
}

%}

%import "../include/SHARPlib/layer.h"
%import "../include/SHARPlib/winds.h"
%import "../include/SHARPlib/parcel.h"
%include "../include/SHARPlib/params/convective.h"
