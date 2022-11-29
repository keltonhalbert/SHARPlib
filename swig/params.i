/* File: params.i */
%module params 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../include/constants.h"
    #include "../include/profile.h"
    #include "../include/utils.h"
    #include "../include/interp.h"
    #include "../include/winds.h"
    #include "../include/params.h"
%}
%import "../include/constants.h"
%import "../include/profile.h"
%import "../include/interp.h"
%import "../include/utils.h"
%import "../include/winds.h"
%include "../include/params.h"
