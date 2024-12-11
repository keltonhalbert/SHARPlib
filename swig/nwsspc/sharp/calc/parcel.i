/* File: parcel.i */
%module parcel 
%{
    #include <SHARPlib/parcel.h>
%}

%import "../include/SHARPlib/layer.h"
%import "../include/SHARPlib/thermo.h"
%include "../include/SHARPlib/parcel.h"

/* We need to instantiate the templated functions for the types we want to use.*/ 
/* Python supports overloading the same function name. */
%template(most_unstable_parcel) sharp::Parcel::most_unstable_parcel<sharp::PressureLayer, sharp::lifter_wobus>;
%template(most_unstable_parcel) sharp::Parcel::most_unstable_parcel<sharp::HeightLayer, sharp::lifter_wobus>;
%template(most_unstable_parcel) sharp::Parcel::most_unstable_parcel<sharp::PressureLayer, sharp::lifter_cm1>;
%template(most_unstable_parcel) sharp::Parcel::most_unstable_parcel<sharp::HeightLayer, sharp::lifter_cm1>;
%template(mixed_layer_parcel) sharp::Parcel::mixed_layer_parcel<sharp::PressureLayer>; 
%template(mixed_layer_parcel) sharp::Parcel::mixed_layer_parcel<sharp::HeightLayer>; 
%template(lift_parcel) sharp::Parcel::lift_parcel<sharp::lifter_wobus>;
%template(lift_parcel) sharp::Parcel::lift_parcel<sharp::lifter_cm1>;
