/* File: profile.i */
%module profile 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/constants.h"
    #include "../include/utils.h"
    #include "../include/interp.h"
    #include "../include/profile.h"
%}
%import "../include/constants.h"
%import "../include/interp.h"
%import "../include/utils.h"
%include "../include/profile.h"
