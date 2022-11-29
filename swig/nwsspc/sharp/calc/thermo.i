/* File: thermo.i */
%module thermo 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/constants.h"
    #include "../include/profile.h"
    #include "../include/utils.h"
    #include "../include/interp.h"
    #include "../include/thermo.h"
%}
/*
%include "./numpy.i"
%init %{
import_array();
%}
*/
%import "../include/constants.h"
%import "../include/profile.h"
%import "../include/interp.h"
%import "../include/utils.h"
%include "../include/thermo.h"
