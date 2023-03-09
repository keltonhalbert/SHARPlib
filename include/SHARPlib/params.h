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
 * motion, especially when considering elevated convection. That method can
 * be found in params.cpp/.h.
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


WindComponents storm_motion_bunkers(const float* pressure, const float* height,
                                    const float* u_wind, const float* v_wind, 
                                    int num_levs, bool leftMover);

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
