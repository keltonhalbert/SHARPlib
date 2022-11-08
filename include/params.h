/**
 * \file
 * \brief Routines used to computed derived sounding parameters from vertical atmospheric profiles.  
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */
#ifndef __SHARP_PARAMS
#define __SHARP_PARAMS

#include "utils.h"

namespace sharp {

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
    float pcl_pres;

    /**
     * \brief Parcel starting temperature (degC)
     */
    float pcl_tmpc;

    /**
     * \brief Parcel starting dewpoint (degC)
     */
    float pcl_dwpc;

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
};


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
