/* File: params.i */
%module params 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
    #include <SHARPlib/params.h>
%}

%include exception.i

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

%import "../include/SHARPlib/layer.h"
%import "../include/SHARPlib/winds.h"
%import "../include/SHARPlib/profile.h"
%import "../include/SHARPlib/parcel.h"
%include "../include/SHARPlib/params.h"
