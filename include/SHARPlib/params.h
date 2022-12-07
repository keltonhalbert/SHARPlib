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


#include <SHARPlib/constants.h>
#include <SHARPlib/utils.h>

namespace sharp {


float energy_helicity_index(float cape, float helicity);


float supercell_composite_parameter(float mu_cape, float eff_srh, 
                                                   float eff_shear);


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
