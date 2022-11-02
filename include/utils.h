/**
 * \file
 * \brief Utilities and structures useful for grouping objects, conducting QC, and things that otherwise don't have a home.  
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-11-02
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */
#pragma once

namespace sharp {

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named floats that represent a height layer. 
 *
 */
struct HeightLayer {

    /**
     * \brief The bottom of the height layer (meters)
     */
    float zbot;

    /**
     * \brief The top of the height layer (meters)
     */
    float ztop;

    /**
     * \brief The height interval with which to iterate over the layer (meters) 
     */
    float dz = 100.0;
};



/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named floats that represent a pressure layer. 
 *
 */
struct PressureLayer {

    /**
     * \brief The bottom of the pressure layer (hPa)
     */
    float pbot;

    /**
     * \brief The top of the pressure layer (hPa)
     */
    float ptop;

    /**
     * \brief The pressure interval with which to iterate over the layer (hPa)
     */
    float dp = 1;
};


}


namespace sharp::exper {


}


