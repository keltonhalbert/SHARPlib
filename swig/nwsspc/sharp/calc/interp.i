/* File: interp.i */
%module interp 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/constants.h"
    #include "../include/profile.h"
    #include "../include/utils.h"
    #include "../include/interp.h"
%}
%import "../include/constants.h"
%import "../include/profile.h"
%import "../include/utils.h"
%include "../include/interp.h"
