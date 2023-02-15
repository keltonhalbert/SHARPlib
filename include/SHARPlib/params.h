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
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>
#include <SHARPlib/profile.h>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Estimates supercell storm motion using the Bunkers et al. 2000 method.
 *
 * Estimates supercell storm motion using the Bunkers et al. 2000 method
 * described in the following paper: 
 *      https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2
 *
 * This does not use any of the updated methods described by Bunkers et al. 2014,
 * which uses Effective Inflow Layer metricks to get better estimates of storm
 * motion, especially when considering elevated convection. That method can
 * be found in params.cpp/.h.
 *
 * TO-DO: Need both left and right motions handled, currently only handles 
 * the right movers.
 *
 * \param prof  sharp::Profile of atmospheric data
 * \return sharp::WindComponents    {u_storm, v_storm}
 */
WindComponents storm_motion_bunkers_np(Profile *prof);



float energy_helicity_index(float cape, float helicity);


float supercell_composite_parameter(float mu_cape, float eff_srh, 
                                                   float eff_shear);


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
