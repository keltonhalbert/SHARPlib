/* File: parcel.i */
%module parcel 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
%}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const int N1),
    (const float height[], const int N2),
    (const float pot_temperature[], const int N3),
    (const float wv_mixratio[], const int N4)
};

%import "../include/SHARPlib/layer.h"
%import "../include/SHARPlib/thermo.h"
%include "../include/SHARPlib/parcel.h"

/* We need to instantiate the templated functions for the types we want to use.*/ 
/* Python supports overloading the same function name. */

%extend sharp::Parcel {

    static sharp::Parcel mixed_layer_parcel(
        const float pressure[], const int N1,
        const float height[], const int N2, 
        const float pot_temperature[], const int N3, 
        const float wv_mixratio[], const int N4,
        sharp::PressureLayer& mix_layer
    ) {
        if (N1 != N2) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        return sharp::Parcel::mixed_layer_parcel(pressure, height, pot_temperature, wv_mixratio, N1, mix_layer);
    }
}
