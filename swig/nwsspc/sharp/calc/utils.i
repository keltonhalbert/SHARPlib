/* File: utils.i */
%module utils 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/constants.h"
    #include "../include/interp.h"
    #include "../include/utils.h"
%}
%import "../include/constants.h"
%import "../include/interp.h"
%include "../include/utils.h"
