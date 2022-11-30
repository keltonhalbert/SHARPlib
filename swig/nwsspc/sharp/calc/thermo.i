/* File: thermo.i */
%module thermo 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/thermo.h"
%}
%import "../include/utils.h"
%include "../include/thermo.h"
