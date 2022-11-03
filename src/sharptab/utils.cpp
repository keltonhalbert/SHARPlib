
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
#include <stdexcept>

#include "utils.h"

namespace sharp {


HeightLayer::HeightLayer(float bot, float top, float delta) {

    if (bot > top) {
        throw std::range_error("The top of the layer must be > the bottom of the layer.");
    }
    zbot = bot;
    ztop = top;
    dz = delta;
}


HeightLayer::HeightLayer(float bot, float top) {
    if (bot > top) {
        throw std::range_error("The top of the layer must be > the bottom of the layer.");
    }
    zbot = bot;
    ztop = top;
    dz = 100.0;
}


PressureLayer::PressureLayer(float bot, float top, float delta) {
    if (bot < top) {
        throw std::range_error("The top of the pressurelayer must be < the bottom of the pressure layer.");
    }
    pbot = bot;
    ptop = top;
    dp = delta;
}


PressureLayer::PressureLayer(float bot, float top) {
    if (bot < top) {
        throw std::range_error("The top of the pressurelayer must be < the bottom of the pressure layer.");
    }
    pbot = bot;
    ptop = top;
    dp = 1;
}







} // end namespace sharp


namespace sharp::exper {

} // end namespace sharp::exper
