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


#include <SHARPlib/constants.h>
#include <SHARPlib/utils.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>
#include <SHARPlib/profile.h>

namespace sharp {


WindComponents storm_motion_bunkers_np(Profile *prof) {
    // derived storm deviation from the mean wind
    // given in m/s 
    float deviation = 7.5;

    // get the pressure values of the AGL heights required
    float pressure_sfc = prof->pres[0];
    float pressure_500m = interp_height(prof->hght[0] + 500.0, 
                                        prof->hght, prof->pres, prof->NZ);
    float pressure_5500m = interp_height(prof->hght[0] + 5500.0, 
                                         prof->hght, prof->pres, prof->NZ);
    float pressure_6000m = interp_height(prof->hght[0] + 6000.0, 
                                         prof->hght, prof->pres, prof->NZ);

    // set up the layers
    PressureLayer layer_lo = {pressure_sfc, pressure_500m};
    PressureLayer layer_hi = {pressure_5500m, pressure_6000m};
    PressureLayer layer_tot = {pressure_sfc, pressure_6000m};


    // get the mean wind of these two layers
    WindComponents mean_wind_0_500m = mean_wind(layer_lo, prof->pres, 
                                                prof->uwin, prof->vwin, 
                                                prof->NZ); 

    WindComponents mean_wind_5500_6000m = mean_wind(layer_hi, prof->pres, 
                                                    prof->uwin, prof->vwin, 
                                                    prof->NZ); 

    WindComponents mean_wind_0_6000m = mean_wind(layer_tot, prof->pres, 
                                                 prof->uwin, prof->vwin, 
                                                 prof->NZ); 

    float shear_u = mean_wind_5500_6000m.u - mean_wind_0_500m.u;
    float shear_v = mean_wind_5500_6000m.v - mean_wind_0_500m.v;
    float mag = vector_magnitude(shear_u, shear_v);

    float right_storm_u = mean_wind_0_6000m.u + ( (deviation / mag) * shear_v); 
    float right_storm_v = mean_wind_0_6000m.v + ( (deviation / mag) * shear_u); 

    return {right_storm_u, right_storm_v};
}


float energy_helicity_index(float cape, float helicity) {
#ifndef NO_QC
    if ((cape == MISSING) || (helicity == MISSING)) {
        return MISSING;
    }
#endif
    return (cape * helicity) / 160000.0;
}


float supercell_composite_parameter(float mu_cape, float eff_srh,
                                                   float eff_shear) {
#ifndef NO_QC
    if ((mu_cape == MISSING) || (eff_srh == MISSING) || (eff_shear == MISSING)) {
        return MISSING;
    }
#endif

    if (eff_shear > 20.0) {
        eff_shear = 20.0;
    }
    else if (eff_shear < 10.0) {
        eff_shear = 0.0;
    }

    float mu_cape_term = mu_cape / 1000.0;
    float eff_srh_term = eff_srh / 50.0;
    float eff_shear_term = eff_shear / 20.0;

    return mu_cape_term * eff_srh_term * eff_shear_term;
}


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


