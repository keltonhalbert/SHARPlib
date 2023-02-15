/* File: params.i */
%module params 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
    #include <SHARPlib/params.h>
%}
%import "../include/SHARPlib/utils.h"
%import "../include/SHARPlib/winds.h"
%import "../include/SHARPlib/profile.h"
%import "../include/SHARPlib/parcel.h"
%include "../include/SHARPlib/params.h"
