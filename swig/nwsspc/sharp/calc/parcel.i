/* File: parcel.i */
%module parcel 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
%}
%import "../include/SHARPlib/profile.h"
%include "../include/SHARPlib/parcel.h"
