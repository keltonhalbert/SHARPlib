/* File: params.i */
%module params 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
    #include <SHARPlib/params/convective.h>
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

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1),
    (const float height[], const std::ptrdiff_t N2),
    (const float temperature[], const std::ptrdiff_t N3),
    (const float dewpoint[], const std::ptrdiff_t N4),
    (const float virtemp_arr[], const std::ptrdiff_t N5)
}

%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1),
    (const float height[], const std::ptrdiff_t N2),
    (const float u_wind[], const std::ptrdiff_t N3),
    (const float v_wind[], const std::ptrdiff_t N4)
}
%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1),
    (const float height[], const std::ptrdiff_t N2),
    (const float temperature[], const std::ptrdiff_t N3),
    (const float mse_arr[], const std::ptrdiff_t N4),
    (const float u_wind[], const std::ptrdiff_t N5),
    (const float v_wind[], const std::ptrdiff_t N6)
}

%rename (effective_inflow_layer) _effective_inflow_layer;
%rename (storm_motion_bunkers) _storm_motion_bunkers;
%rename (entrainment_cape) _entrainment_cape;
%ignore effective_inflow_layer;
%ignore storm_motion_bunkers;
%ignore entrainment_cape;

%inline %{

sharp::PressureLayer _effective_inflow_layer(
                        sharp::lifter_wobus& lifter,
                        const float pressure[], const std::ptrdiff_t N1,
                        const float height[], const std::ptrdiff_t N2,
                        const float temperature[], const std::ptrdiff_t N3,
                        const float dewpoint[], const std::ptrdiff_t N4,
                        const float virtemp_arr[], const std::ptrdiff_t N5,
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
    float* pcl_vtmp_arr = (float *)malloc(N1*sizeof(float));
    if ((buoy == NULL) || (pcl_vtmp_arr == NULL)) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for temporary array of size %d.", 
            N1
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    sharp::PressureLayer eil = sharp::effective_inflow_layer(lifter, 
        pressure, height, temperature, dewpoint, virtemp_arr, 
        pcl_vtmp_arr, buoy, N1, cape_thresh, cinh_thresh, mupcl);
    free(buoy);
    free(pcl_vtmp_arr);
    return eil;
}

sharp::PressureLayer _effective_inflow_layer(
                        sharp::lifter_cm1& lifter,
                        const float pressure[], const std::ptrdiff_t N1,
                        const float height[], const std::ptrdiff_t N2,
                        const float temperature[], const std::ptrdiff_t N3,
                        const float dewpoint[], const std::ptrdiff_t N4,
                        const float virtemp_arr[], const std::ptrdiff_t N5,
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
    float* pcl_vtmp_arr = (float *)malloc(N1*sizeof(float));
    if ((buoy == NULL) || (pcl_vtmp_arr == NULL)) {
        PyErr_Format(
            PyExc_MemoryError, 
            "Could not allocate memory for temporary array of size %d.", 
            N1
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    sharp::PressureLayer eil = sharp::effective_inflow_layer(lifter, 
        pressure, height, temperature, dewpoint, virtemp_arr, 
        pcl_vtmp_arr, buoy, N1, cape_thresh, cinh_thresh, mupcl);
    free(buoy);
    return eil;
}

sharp::WindComponents _storm_motion_bunkers(
                            const float pressure[], const std::ptrdiff_t N1,
                            const float height[], const std::ptrdiff_t N2,
                            const float u_wind[], const std::ptrdiff_t N3,
                            const float v_wind[], const std::ptrdiff_t N4,
                            sharp::HeightLayer mean_wind_layer_agl,
                            sharp::HeightLayer wind_shear_layer_agl,
                            const bool leftMover = false,
                            const bool pressureWeighted = false) {
    if ((N1 != N2) || (N1 != N3) || (N1 != N4)) {
        PyErr_Format(PyExc_ValueError, 
            "Arrays must be same lenght, insead got (%d, %d, %d, %d)",
            N1, N2, N3, N4
        );
        return {sharp::MISSING, sharp::MISSING};
    }

    return sharp::storm_motion_bunkers(pressure, height, u_wind, v_wind, N1,
                                       mean_wind_layer_agl, wind_shear_layer_agl,
                                       leftMover, pressureWeighted);
}

sharp::WindComponents _storm_motion_bunkers(
                            const float pressure[], const std::ptrdiff_t N1,
                            const float height[], const std::ptrdiff_t N2,
                            const float u_wind[], const std::ptrdiff_t N3,
                            const float v_wind[], const std::ptrdiff_t N4,
                            sharp::PressureLayer eff_lyr, 
                            sharp::Parcel* mupcl, 
                            const bool leftMover=false) {
    if ((N1 != N2) || (N1 != N3) || (N1 != N4)) {
        PyErr_Format(PyExc_ValueError, 
            "Arrays must be same lenght, insead got (%d, %d, %d, %d)",
            N1, N2, N3, N4
        );
        return {sharp::MISSING, sharp::MISSING};
    }
    return sharp::storm_motion_bunkers(pressure, height, u_wind, v_wind, N1,
                                       eff_lyr, mupcl, leftMover);
}

float _entrainment_cape(const float pressure[], const std::ptrdiff_t N1,
                        const float height[], const std::ptrdiff_t N2,
                        const float temperature[], const std::ptrdiff_t N3,
                        const float mse_arr[], const std::ptrdiff_t N4,
                        const float u_wind[], const std::ptrdiff_t N5,
                        const float v_wind[], const std::ptrdiff_t N6,
                        sharp::Parcel* pcl) {
	if ( (N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5) || (N1 != N6) ) {
        PyErr_Format(PyExc_ValueError, 
            "Arrays must be same lenght, insead got (%d, %d, %d, %d, %d, %d)",
            N1, N2, N3, N4, N5, N6
        );
        return sharp::MISSING; 
	}
        return sharp::entrainment_cape(pressure, height, temperature, mse_arr,
                                       u_wind, v_wind, N1, pcl);
}

%} /* end inline */

%import "../../include/SHARPlib/layer.h"
%import "../../include/SHARPlib/winds.h"
%import "../../include/SHARPlib/parcel.h"
%include "../../include/SHARPlib/params/convective.h"
