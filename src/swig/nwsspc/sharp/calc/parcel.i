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

// output array typemap
%apply (float** ARGOUTVIEWM_ARRAY1, int* DIM1) {
    (float** out_arr, int* NOUT)
};

// output value typemaps
%apply float& OUTPUT { float& lfc_pres };
%apply float& OUTPUT { float& el_pres };
%apply float& OUTPUT { float& cape };
%apply float& OUTPUT { float& cinh };

// mixed_layer_parcel typemap
%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1),
    (const float height[], const std::ptrdiff_t N2),
    (const float pot_temperature[], const std::ptrdiff_t N3),
    (const float wv_mixratio[], const std::ptrdiff_t N4)
};

// most_unstable_parcel typemap
%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1),
    (const float height[], const std::ptrdiff_t N2),
    (const float temperature[], const std::ptrdiff_t N3),
    (const float virtemp[], const std::ptrdiff_t N4),
    (const float dewpoint[], const std::ptrdiff_t N5)
}

// lift_parcel typemape
%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1)
};

// find_lfc_el and cape_cinh typemap
%apply (float* IN_ARRAY1, int DIM1) {
    (const float pressure[], const std::ptrdiff_t N1),
    (const float height[], const std::ptrdiff_t N2),
    (const float buoyancy[], const std::ptrdiff_t N3)
}

// don't wrap the current C++ version -- 
// we need to restructure the function for numpy arrays
%ignore cape_cinh;
%ignore find_lfc_el;

// this allows for us to call our wrapped code below
%rename("%s") find_lfc_el;
%rename("%s") cape_cinh;

%import "../../include/SHARPlib/layer.h"
%import "../../include/SHARPlib/thermo.h"
%include "../../include/SHARPlib/parcel.h"

%extend sharp::Parcel {

    static sharp::Parcel mixed_layer_parcel(
        sharp::PressureLayer& mix_layer,
        const float pressure[], const std::ptrdiff_t N1,
        const float height[], const std::ptrdiff_t N2, 
        const float pot_temperature[], const std::ptrdiff_t N3, 
        const float wv_mixratio[], const std::ptrdiff_t N4
    ) {
        if (N1 != N2) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        return sharp::Parcel::mixed_layer_parcel(
            mix_layer,
            pressure, 
            height, 
            pot_temperature, 
            wv_mixratio, 
            N1
        );
    }

    static sharp::Parcel mixed_layer_parcel(
        sharp::HeightLayer& mix_layer,
        const float pressure[], const std::ptrdiff_t N1,
        const float height[], const std::ptrdiff_t N2, 
        const float pot_temperature[], const std::ptrdiff_t N3, 
        const float wv_mixratio[], const std::ptrdiff_t N4
    ) {
        if (N1 != N2) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        return sharp::Parcel::mixed_layer_parcel(
            mix_layer,
            pressure, 
            height, 
            pot_temperature, 
            wv_mixratio, 
            N1
        );
    }

    static sharp::Parcel most_unstable_parcel(
        sharp::lifter_wobus& lifter,
        sharp::PressureLayer& search_layer,
        const float pressure[], const std::ptrdiff_t N1, 
        const float height[], const std::ptrdiff_t N2,
        const float temperature[], const std::ptrdiff_t N3, 
        const float virtemp[], const std::ptrdiff_t N4, 
        const float dewpoint[], const std::ptrdiff_t N5) {
        if ((N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5)) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        float* buoy_arr = (float *)malloc(N1*sizeof(float));
        float* pcl_vtmpk_arr = (float *)malloc(N1*sizeof(float));
        if ((buoy_arr == NULL) || (pcl_vtmpk_arr == NULL)) {
            PyErr_Format(
                PyExc_MemoryError,
                "Could not allocate memory for output array of size %d.",
                N1
            );
            return sharp::Parcel();
        }

        const sharp::Parcel mu_pcl = 
        sharp::Parcel::most_unstable_parcel(
            search_layer,
            lifter,
            pressure, 
            height, 
            temperature, 
            virtemp, 
            dewpoint, 
            pcl_vtmpk_arr,
            buoy_arr, 
            N1
        );
        delete[] buoy_arr;
        delete[] pcl_vtmpk_arr;

        return mu_pcl;
    }

    static sharp::Parcel most_unstable_parcel(
        sharp::lifter_wobus& lifter,
        sharp::HeightLayer& search_layer,
        const float pressure[], const std::ptrdiff_t N1, 
        const float height[], const std::ptrdiff_t N2,
        const float temperature[], const std::ptrdiff_t N3, 
        const float virtemp[], const std::ptrdiff_t N4, 
        const float dewpoint[], const std::ptrdiff_t N5) {
        if ((N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5)) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        float* buoy_arr = (float *)malloc(N1*sizeof(float));
        float* pcl_vtmpk_arr = (float *)malloc(N1*sizeof(float));
        if ((buoy_arr == NULL) || (pcl_vtmpk_arr == NULL)) {
            PyErr_Format(
                PyExc_MemoryError,
                "Could not allocate memory for output array of size %d.",
                N1
            );
            return sharp::Parcel();
        }

        const sharp::Parcel mu_pcl = 
        sharp::Parcel::most_unstable_parcel(
            search_layer,
            lifter,
            pressure, 
            height, 
            temperature, 
            virtemp, 
            dewpoint, 
            pcl_vtmpk_arr,
            buoy_arr, 
            N1
        );
        delete[] buoy_arr;
        delete[] pcl_vtmpk_arr;

        return mu_pcl;
    }

    static sharp::Parcel most_unstable_parcel(
        sharp::lifter_cm1& lifter,
        sharp::PressureLayer& search_layer,
        const float pressure[], const std::ptrdiff_t N1, 
        const float height[], const std::ptrdiff_t N2,
        const float temperature[], const std::ptrdiff_t N3, 
        const float virtemp[], const std::ptrdiff_t N4, 
        const float dewpoint[], const std::ptrdiff_t N5) {
        if ((N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5)) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        float* buoy_arr = (float *)malloc(N1*sizeof(float));
        float* pcl_vtmpk_arr = (float *)malloc(N1*sizeof(float));
        if ((buoy_arr == NULL) || (pcl_vtmpk_arr == NULL)) {
            PyErr_Format(
                PyExc_MemoryError,
                "Could not allocate memory for output array of size %d.",
                N1
            );
            return sharp::Parcel();
        }

        const sharp::Parcel mu_pcl = 
        sharp::Parcel::most_unstable_parcel(
            search_layer,
            lifter,
            pressure, 
            height, 
            temperature, 
            virtemp, 
            dewpoint, 
            pcl_vtmpk_arr,
            buoy_arr, 
            N1
        );
        delete[] buoy_arr;
        delete[] pcl_vtmpk_arr;

        return mu_pcl;
    }

    static sharp::Parcel most_unstable_parcel(
        sharp::lifter_cm1& lifter,
        sharp::HeightLayer& search_layer,
        const float pressure[], const std::ptrdiff_t N1, 
        const float height[], const std::ptrdiff_t N2,
        const float temperature[], const std::ptrdiff_t N3, 
        const float virtemp[], const std::ptrdiff_t N4, 
        const float dewpoint[], const std::ptrdiff_t N5) {
        if ((N1 != N2) || (N1 != N3) || (N1 != N4) || (N1 != N5)) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d)",
                  N1, N2
            );
            return sharp::Parcel();
        }

        float* buoy_arr = (float *)malloc(N1*sizeof(float));
        float* pcl_vtmpk_arr = (float *)malloc(N1*sizeof(float));
        if ((buoy_arr == NULL) || (pcl_vtmpk_arr == NULL)) {
            PyErr_Format(
                PyExc_MemoryError,
                "Could not allocate memory for output array of size %d.",
                N1
            );
            return sharp::Parcel();
        }

        const sharp::Parcel mu_pcl = 
        sharp::Parcel::most_unstable_parcel(
            search_layer,
            lifter,
            pressure, 
            height, 
            temperature, 
            virtemp, 
            dewpoint, 
            pcl_vtmpk_arr,
            buoy_arr, 
            N1
        );
        delete[] buoy_arr;
        delete[] pcl_vtmpk_arr;

        return mu_pcl;
    }

    void lift_parcel(sharp::lifter_wobus& liftpcl, 
                    const float pressure[], const std::ptrdiff_t N1,
                    float** out_arr, int* NOUT) {

        float* temp = (float *)malloc(N1*sizeof(float));
        if (temp == NULL) {
            PyErr_Format(
                PyExc_MemoryError,
                "Could not allocate memory for output array of size %d.",
                N1
            );
            return;
        }

        *out_arr = temp;
        *NOUT = N1;

        $self->lift_parcel(liftpcl, pressure, *out_arr, N1);
    }

    void lift_parcel(sharp::lifter_cm1& liftpcl, 
                    const float pressure[], const std::ptrdiff_t N1,
                    float** out_arr, int* NOUT) {

        float* temp = (float *)malloc(N1*sizeof(float));
        if (temp == NULL) {
            PyErr_Format(
                PyExc_MemoryError,
                "Could not allocate memory for output array of size %d.",
                N1
            );
            return;
        }

        *out_arr = temp;
        *NOUT = N1;

        $self->lift_parcel(liftpcl, pressure, *out_arr, N1);
    }

    void find_lfc_el(const float pressure[], const std::ptrdiff_t N1, 
                     const float height[], const std::ptrdiff_t N2, 
                     const float buoyancy[], const std::ptrdiff_t N3,
                     float& lfc_pres, float& el_pres) {
        if ((N1 != N2) || (N1 != N3)) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d, %d)",
                  N1, N2, N3
            );
            return;
        }

        $self->find_lfc_el(pressure, height, buoyancy, N1);
        lfc_pres = $self->lfc_pressure;
        el_pres = $self->eql_pressure;
    }

    void cape_cinh(const float pressure[], const std::ptrdiff_t N1,
                   const float height[], const std::ptrdiff_t N2,
                   const float buoyancy[], const std::ptrdiff_t N3,
                   float& cape, float& cinh) {
        if ((N1 != N2) || (N1 != N3)) {
            PyErr_Format(
                PyExc_ValueError, 
                "Arrays must be same lenght, insead got (%d, %d, %d)",
                  N1, N2, N3
            );
            return;
        }

        $self->cape_cinh(pressure, height, buoyancy, N1);
        cape = $self->cape;
        cinh = $self->cinh;
    }

}
