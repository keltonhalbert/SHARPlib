// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
#pragma once

namespace sharp {

inline float wobf(float temperature);
inline float vappres(float temperature);
inline float lcltemp(float temperature, float dewpoint);
inline float temperature_at_mixratio(float mixratio, float pressure);
inline float thalvl(float potential_temperature, float temperature);
inline float theta(float pressure, float temperature, float ref_pressure);


}
