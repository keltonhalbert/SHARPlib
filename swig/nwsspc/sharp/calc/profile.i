/* File: profile.i */
%module profile 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/constants.h"
    #include "../include/utils.h"
    #include "../include/interp.h"
    #include "../include/thermo.h"
    #include "../include/winds.h"
    #include "../include/profile.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%numpy_typemaps(float, NPY_DOUBLE, int)
%apply (float *IN_ARRAY1, int DIM1) {
            (float *pres, int NZ1), (float *hght, int NZ2),
            (float *tmpc, int NZ3), (float *dwpc, int NZ4),
            (float *wnd1, int NZ5), (float *wnd2, int NZ6)
        }

%import "../include/constants.h"
%import "../include/interp.h"
%import "../include/utils.h"
%import "../include/thermo.h"
%import "../include/winds.h"



%rename (create_profile) profile_wrap;
%ignore create_profile;
%exception profile_wrap {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \breif Creates a Profile from arrays
 *
 * \param pres              Array of pressure (hPa)
 * \param hght              Array of height (meters)
 * \param tmpc              Array of temperature (degC)
 * \param dwpc              Array of dewpoiont (degC)
 * \param wspd_or_u         Array of wind speed or u-component (kts or m/s)
 * \param wdir_or_v         Array of wind direction or v-component (deg or m/s)
 * \param NZ                The number of vertical levels
 * \param sounding_type     The sharp::Source of the profile
 * \param windComponents    A boolean if winds are vectors or components
 */
sharp::Profile* profile_wrap(float *pres, int NZ1, float *hght, int NZ2,
                             float *tmpc, int NZ3, float *dwpc, int NZ4,
                             float *wnd1, int NZ5, float *wnd2, int NZ6,
                             sharp::Source sounding_type, bool windComponents) {
    if ( (NZ1 != NZ2) || (NZ1 != NZ2) || (NZ1 != NZ3) ||
         (NZ1 != NZ4) || (NZ1 != NZ5) || (NZ1 != NZ6) ) {
        PyErr_Format(PyExc_ValueError, "Arrays must be same lenght");
        return nullptr;
    }

    return create_profile(pres, hght, tmpc, dwpc, wnd1, wnd2, NZ1,
                                 sounding_type, windComponents);
}
%}

%include "../include/profile.h"
