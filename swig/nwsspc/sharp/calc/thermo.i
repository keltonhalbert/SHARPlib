/* File: thermo.i */
%module thermo 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/thermo.h>
%}

/* Import Numpy Array Functionality */
%include "numpy.i"
%init %{
import_array();
%}

%import "../include/SHARPlib/layer.h"

%apply float &OUTPUT { float& pressure_at_lcl };
%apply float &OUTPUT { float& temperature_at_lcl };

/* Set up an argument typemap for our arrays */
%apply (float* IN_ARRAY1, int DIM1) {
        (const float hght_array[], const int N1),
        (const float tmpc_array[], const int N2)
    }

%apply (float* IN_ARRAY1, int DIM1) {
        (const float pres_array[], const int N1),
        (const float hght_array[], const int N2),
        (const float tmpc_array[], const int N3)
    }

%apply (float* IN_ARRAY1, int DIM1) {
		(const float pressure[], const int N1),
		(const float temperature[], const int N2)
}
%apply (float* INPLACE_ARRAY1, int DIM1) {
		(float mixr_arr[], int N3)
}

/**
 *  Wrap the lapse_rate function in such a way
 *  that we can use our numpy %apply map, and make
 *  sure arrays are of the same size
 */
%rename (lapse_rate) _lapse_rate;
%rename (mixratio) _mixratio;

%ignore lapse_rate;
%ignore mixratio;

%exception _mixratio {
	$action
	if (PyErr_Occurred()) SWIG_fail;
}

%exception _lapse_rate {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

float _lapse_rate(sharp::HeightLayer layer, 
                  const float hght_array[], const int N1, 
                  const float tmpc_array[], const int N2) {
    if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
        return sharp::MISSING;
    }

    return sharp::lapse_rate(layer, hght_array, tmpc_array, N1);
}

float _lapse_rate(sharp::PressureLayer layer, 
                  const float pres_array[], const int N1, 
                  const float hght_array[], const int N2,
                  const float tmpc_array[], const int N3) {
    if ((N1 != N2) || (N1 != N3)) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d, %d)",
            N1, N2, N3
        );
        return sharp::MISSING;
    }

    return sharp::lapse_rate(layer, pres_array, hght_array, tmpc_array, N1);
}

float _mixratio(float pressure, float temperature) {
	return sharp::mixratio(pressure, temperature);
}

void _mixratio(const float pressure[], const int N1, const float temperature[],
               const int N2, float mixr_arr[], const int N3) {
        if (N1 != N2) {
        PyErr_Format(
            PyExc_ValueError, "Arrays must be same length, got (%d, %d)",
            N1, N2
        );
		return; 
    }
	
	for (int i = 0; i < N1; ++i) {
		mixr_arr[i] = sharp::mixratio(pressure[i], temperature[i]);
	};

	return ;
}

%}

%include "../include/SHARPlib/thermo.h"
