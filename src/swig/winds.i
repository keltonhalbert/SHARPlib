/* File: winds.i */
%module winds 
%{
    #define SWIG_FILE_WITH_INIT
    #include "../../include/constants.h"
    #include "../../include/profile.h"
    #include "../../include/utils.h"
    #include "../../include/interp.h"
    #include "../../include/winds.h"
%}
%import "../../include/constants.h"
%import "../../include/profile.h"
%import "../../include/interp.h"
%import "../../include/utils.h"
%include "../../include/winds.h"
