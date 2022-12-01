/* File: params.i */
%module params 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/params.h>
%}
%include "../include/SHARPlib/params.h"
