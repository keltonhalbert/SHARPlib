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

%import "../include/SHARPlib/layer.h"

%include exception.i

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

/* Set up argument typemap */
%apply(float IN_ARRAY1, int DIM1){
	(const float pressure[], const int N1),
	(const float temperature[], const int N2),
	(const float dewpoint[], const int N3),
	(const float wv_mixratio[], const int N4),
	(const float theta_arr[], const int N5),
	(const float thetae[], const int N6)
}

%apply (float IN_ARRAY1, int DIM1) {
	(const float pres_arr[], const int N1),
	(const float hght_arr[], const int N2),
	(const float buoy_arr[], const int N3)
}


%rename(define_parcel) _define_parcel;
%rename(find_lfc_el) _find_lfc_el;
%rename(cape_cinh) _cape_cinh;

%ignore define_parcel;
%ignore find_lfc_el;
%ignore cape_cinh;

%inline %{

void _define_parcel(const float pressure[], const int N1, 
					const float temperature[], const int N2,
                    const float dewpoint[], const int N3,
					const float wv_mixratio[], const int N4,
					const float theta_arr[], const int N5,
					const float thetae[], const int N6,
					sharp::Parcel& pcl, sharp::LPL source) {
	if ( (N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5) || (N1 != N6) ) {
        PyErr_Format(PyExc_ValueError, 
                     "Arrays must be same lenght, insead got (%d, %d, %d, %d, %d, %d)",
                     N1, N2, N3, N4, N5, N6);
	}
	sharp::define_parcel(pressure, temperature, dewpoint, wv_mixratio,
						 theta_arr, thetae, N1, pcl, source);
}

void find_lfc_el(sharp::Parcel* pcl, const float pres_arr[], const int N1,
                 const float hght_arr[], const int N2, const float buoy_arr[],
                 const int N3){
	if ( (N1 != N2) || (N1 != N3) ) {
        PyErr_Format(PyExc_ValueError, 
                     "Arrays must be same lenght, insead got (%d, %d, %d)",
                     N1, N2, N3);
	}
	sharp::find_lfc_el(pcl, pres_arr, hght_arr, buoy_arr, N1);
}

void cape_cinh(const float pres_arr[], const int N1, 
               const float hght_arr[], const int N2,
			   const float buoy_arr[], const int N3,
			   sharp::Parcel* pcl) {
	if ( (N1 != N2) || (N1 != N3) ) {
		PyErr_Format(PyExc_ValueError, 
					 "Arrays must be same lenght, insead got (%d, %d, %d)",
					 N1, N2, N3);
	}
	sharp::cape_cinh(pres_arr, hght_arr, buoy_arr, N1, pcl);
}

%}
















%include "../include/SHARPlib/parcel.h"
