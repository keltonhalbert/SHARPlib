
/**
 * \file
 * \brief Routines used for parcel lifting and integration 
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-11-09
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */

namespace sharp {

// Trying to mentally sketch out how I want parcel lifting
// and CAPE/CINH integration to work. I want it to be modular
// such that different methods of computing moist adiabats can
// be supported. I want to separate the parcel lifting from the
// actual numerical integration. I want to be able to store and
// reuse the temperature/pressure traces of parcel lifting. Need 
// to support a "fast CAPE/CINH" for effective inflow calculations.
// 
//
// 1) Create the parcel struct
// 2) Define its starting attributes (MU/ML/SFC/etc)
// 3) Lift the parcel/compute temperature trace
// 4) Set LCL/LFC/EL values
// 5) Integrate CAPE/CINH

} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


