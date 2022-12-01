/* File: constants.i */
%module constants
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/constants.h>
%}
%include "../include/SHARPlib/constants.h"
