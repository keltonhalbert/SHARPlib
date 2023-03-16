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
#include <SHARPlib/parcel.h>

namespace sharp {


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Computes the Effective Inflow Layer for a given atmospheric sounding.
 *
 * Computes the Effective Inflow Layer, or the layer of the atmosphere believed
 * to be the primary source of inflow for supercell thunderstorms. The Effective
 * Inflow Layer, and its use in computing shear and storm relative helicity, is 
 * described by Thompson et al. 2007:
 *     https://www.spc.noaa.gov/publications/thompson/effective.pdf
 *
 * Standard/default values for cape_thresh and cinh_thresh have been 
 * experimentally determined to be cape_thresh = 100 J/kg and 
 * cinh_thresh = -250 J/kg.
 * 
 * \param Profile       sharp::Profile of atmospheric data
 * \param cape_thresh   The CAPE threshold that defines the EIL. Default is 100 J/kg.
 * \param cinh_thresh   The CINH threshold that defines the EIL. Default is -250 J/kg.
 */
PressureLayer effective_inflow_layer(Profile *prof, float cape_thresh, 
                                     float cinh_thresh);

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
 * motion, especially when considering elevated convection.
 *
 * \param pressure  Vertical array of Air Pressure (hPa)
 * \param height    Vertical array of heights (meters)
 * \param u_wind    Vertical array of U wind components (m/s)
 * \param v_wind    Vertical array of V wind components (m/s)
 * \param num_levs  The number of vertical levels in the arrays 
 * \param leftMover Boolean flag to return left or right mover
 * \return sharp::WindComponents    {u_storm, v_storm}
 */
WindComponents storm_motion_bunkers_np(const float *pressure,
        const float *height, const float *u_wind, const float *v_wind,
        int num_levs, bool leftMover);


/**
 *  \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 *  \brief Estimates supercell storm motion by using the Bunkers et al. 2014 method.
 *
 *  Estimates supercell storm motion using the Bunkers et al. 2014 method
 *  described in the following paper:
 *      https://doi.org/10.1175/WAF-D-21-0153.1
 *
 *  This method is parcel based, using a mean-wind vector defined as the 
 *  pressure-weighted mean wind between the Effective Inflow Layer surface 
 *  (see sharp::effective_inflow_layer) and 65% of the depth between that 
 *  surface and the most-unstable parcel's Equilibrium Level. This method
 *  produces the same storm motion estimate for surface based supercells,
 *  and captures the motion of elevated supercells much better than 
 *  the Bunkers 2000 method. 
 *
 *  \param Profile      sharp::Profile of atmospheric data
 *  \param leftMover    Boolean flag for left or right deviant supercell
 */
WindComponents storm_motion_bunkers(Profile* prof, bool leftMover);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Computes Entrainment CAPE using a previously lifted parcel. 
 *
 *
 * Computes Entrainment CAPE, or ECAPE, as described by Peters et al. 2023,
 * "An analytic formula for entraining CAPE in mid-latitude storm environments".
 *
 * \param prof      A sharp::Profile of atmospheric data
 * \param pcl       A sharp::Parcel with its sharp::LPL/attributes defined.
 */
float entrainment_cape(Profile* prof, Parcel *pcl);


float energy_helicity_index(float cape, float helicity);


float supercell_composite_parameter(float mu_cape, float eff_srh, 
                                                   float eff_shear);


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
