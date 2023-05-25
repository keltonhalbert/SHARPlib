/* File: parcel.i */
%module parcel 
%{
    #define SWIG_FILE_WITH_INIT
    #include <SHARPlib/parcel.h>
%}
%import "../include/SHARPlib/layer.h"

%include exception.i

%exception {
    try {
        $action
    }
    catch (const std::runtime_error& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

%import "../include/SHARPlib/profile.h"
%include "../include/SHARPlib/parcel.h"
