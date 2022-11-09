
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
#include "constants.h"
#include "interp.h"
#include "parcel.h"
#include "utils.h"

namespace sharp {


Parcel::Parcel() {
    pcl_pres = MISSING;
    pcl_tmpc = MISSING;
    pcl_dwpc = MISSING;
    
    lcl_pressure = MISSING;
    lfc_pressure = MISSING;
    eql_pressure = MISSING;
    mpl_pressure = MISSING;

    cape = 0.0;
    cinh = 0.0;
}


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

void define_parcel(Profile* prof, Parcel* pcl, LPL source) {
    pcl->source = source;

    switch(source) {
        case LPL::SFC:
            // set the parcel attributes to the surface
            pcl->pcl_pres = prof->pres[0];
            pcl->pcl_tmpc = prof->tmpc[0];
            pcl->pcl_dwpc = prof->dwpc[0]; 
            break;
        case LPL::FCST:
            // call forecast parcel routine
            break;
        case LPL::MU:
            // call the most unstable parcel routine
            break;
        case LPL::ML:
            // call the mixed layer parcel routine
            break;
        case LPL::USR:
            // Parcel attributes are user defined, do
            // nothing then?
            break;
        case LPL::EIL: 
            // Call the mean effective inflow layer parcel
            break;
    };
}






} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


