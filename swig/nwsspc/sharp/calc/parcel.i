/* File: parcel.i */
%module parcel 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/parcel.h"
%}
%import "../include/profile.h"
%include "../include/parcel.h"
