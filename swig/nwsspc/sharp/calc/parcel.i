/* File: parcel.i */
%module parcel 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
%}

%include exception.i

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
        return NULL;
    }
}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}
%import "../include/SHARPlib/parcel.h"

%apply (float* IN_ARRAY1, int DIM1) {
	(const float pressure[], const int N1),
	(const float height[], const int N2),
  (const float pot_temperature[], const int N3),
  (const float wv_mixratio[], const int N4)
}

%ignore sharp::Parcel::mixed_layer_parcel;

%rename("%s") sharp::Parcel::mixed_layer_parcel;

%extend sharp::Parcel {
    static sharp::Parcel sharp::Parcel::mixed_layer_parcel(
        const float pressure[], const int N1,
        const float height[], const int N2,
        const float pot_temperature[], const int N3,
        const float wv_mixratio[], const int N4,
        sharp::PressureLayer mix_lyr
    ) {
        if ((N1 != N2) || (N1 != N3) || (N1 != N4)) {
            PyErr_Format(
                PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d, %d)",
                N1, N2, N3, N4
            );
            return sharp::Parcel(); 
        }

        return sharp::Parcel::mixed_layer_parcel(
            pressure, height, pot_temperature, 
            wv_mixratio, N1, mix_lyr
        );
    }
}

%import "../include/SHARPlib/thermo.h"
%import "../include/SHARPlib/layer.h"

%include "../include/SHARPlib/parcel.h"
