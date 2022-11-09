
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
#ifndef __SHARP_PARCEL
#define __SHARP_PARCEL

#include "profile.h"
#include "utils.h"

namespace sharp {


/**
 * \brief Enum that defines the lifted parcel level (LPL) of origin.
 */
enum class LPL : int {
    SFC  = 1, // surface
    FCST = 2, // forecast surface
    MU   = 3, // most unstable
    ML   = 4, // 100mb mixed layer
    USR  = 5, // user-defined
    EIL  = 6, // Mean effective inflow layer
};


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Defines a parcel to be lifted
 *
 * Contains information about a Parcel's starting level and
 * thermodynamic attributes, as well as paramaters computed
 * using the parcel.  
 */
struct Parcel {

    /**
     * \brief Parcel starting pressure (hPa)
     */
    float pres;

    /**
     * \brief Parcel starting temperature (degC)
     */
    float tmpc;

    /**
     * \brief Parcel starting dewpoint (degC)
     */
    float dwpc;

    /**
     * \brief Pressure at the Lifted Condensation Level (hPa)
     */
    float lcl_pressure;

    /**
     * \brief Pressure at the Level of Free Convection (hPa)
     */
    float lfc_pressure;

    /**
     * \brief Pressure at the parcel Equilibrium Level (hPa)
     */
    float eql_pressure;

    /**
     * \brief Pressure at the Maximum Parcel Level (hPa)
     */
    float mpl_pressure;

    /**
     * \brief Parcel Convective Available Potential Energy (J/kg)
     */
    float cape;

    /**
     * \brief Parcel Convective Inhibition (J/kg)
     */
    float cinh;

    /**
     * \brief The type of parcel this is
     */
    LPL source; 
    
    /**
     * \brief Parcel empty constructor that sets all values to sharp::MISSING
     */
    Parcel();
};


/**
 *
 */
void define_parcel(Profile* prof, Parcel* pcl, LPL source);


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
