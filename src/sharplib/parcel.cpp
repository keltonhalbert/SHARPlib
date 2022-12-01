
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
#include "constants.h"
#include "interp.h"
#include "thermo.h"
#include "utils.h"
#include "profile.h"
#include "parcel.h"

namespace sharp {


Parcel::Parcel() {
    // Set LPL values to MISSING
    pres = MISSING;
    tmpc = MISSING;
    dwpc = MISSING;
    
    // Set derived values to MISSING
    lcl_pressure = MISSING;
    lfc_pressure = MISSING;
    eql_pressure = MISSING;
    mpl_pressure = MISSING;

    cape = 0.0;
    cinh = 0.0;
}


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the parcel initial values to the surface of the profile
 *
 * \param prof
 * \param pcl
 *
 */
void _sfc(Profile* prof, Parcel* pcl) {
            pcl->pres = prof->pres[0];
            pcl->tmpc = prof->tmpc[0];
            pcl->dwpc = prof->dwpc[0]; 
}


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the parcel initial values to the most unstable parcel level 
 *
 * \param prof
 * \param pcl
 *
 */
void _mu(Profile* prof, Parcel* pcl) {
    // Search for the most unstable parcel in the bottom
    // 300 hPa of the profile
    PressureLayer mu_layer(prof->pres[0], prof->pres[0]-300.0);

    // max_value returns the max, and will set the pressure
    // of the max via a pointer to a float. 
    max_value(mu_layer, prof->pres, prof->theta_e, prof->NZ, &(pcl->pres));
    pcl->tmpc = interp_pressure(pcl->pres, prof->pres, prof->tmpc, prof->NZ);
    pcl->dwpc = interp_pressure(pcl->pres, prof->pres, prof->dwpc, prof->NZ);
}


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the parcel initial values to the bottom 100mb mixed layer 
 *
 * \param prof
 * \param pcl
 *
 */
void _ml(Profile* prof, Parcel* pcl) {
    PressureLayer mix_layer(prof->pres[0], prof->pres[0]-100.0);

    // get the mean attributes of the lowest 100 hPa
    float mean_mixr = mean_value(mix_layer, prof->pres, 
                                 prof->mixr, prof->NZ);
    float mean_thta = mean_value(mix_layer, prof->pres, 
                                 prof->theta, prof->NZ);

    // set the parcel attributes
    pcl->pres = prof->pres[0];
    pcl->tmpc = theta(1000.0, mean_thta, prof->pres[0]);
    pcl->dwpc = temperature_at_mixratio(mean_mixr, prof->pres[0]); 
}


void define_parcel(Profile* prof, Parcel* pcl, LPL source) {
    pcl->source = source;

    if (source == LPL::SFC) {
        _sfc(prof, pcl);
        return;
    }
    else if (source == LPL::FCST) {
        // TO-DO: Write the forecast_surface routine
        return;
    }
    else if (source == LPL::MU) {
        _mu(prof, pcl);
        return;
    }
    else if (source == LPL::ML) {
        _ml(prof, pcl);
        return;
    }
    else if (source == LPL::EIL) {
        // TO-DO: Write the EIL routine
        return;
    }
    else if (source == LPL::USR) {
        // do nothing - its already been set!
        return;
    }
    else {
        // TO-DO: probably should raise an error or something
        return;
    }
}


float cinh_below_lcl(Profile* prof, Parcel* pcl, float pres_lcl, float tmpc_lcl) {
    // get the virtual temoerature of the LCL
    float vtmp_lcl = virtual_temperature(pres_lcl, tmpc_lcl, tmpc_lcl);

    // parcel thetav, or virtual potential temperature,
    // is constant from the LPL to the LCL.
	float pcl_thetav = theta(pres_lcl, vtmp_lcl, 1000.0);

	// Accumulate CINH in the mixing layer below the LCL.
	// This will be done in 10 mb increments and will use the
	// virtual temperature correction. 
    float cinh = 0.0;
	for (float pbot = pcl->pres; pbot > pres_lcl; pbot -= 10.0) {
		float ptop = pbot - 10.0;
		// don't accidentally go above the LCL
		if (ptop < pres_lcl) ptop = pres_lcl;

        // _pcb - parcel bottom layer
        // _pct - parcel top layer
        // _enb - environment bottom layer
        // _ent - environment top layer

        // "lift" the thetav of the parcel to a new pressure level
		float vtmp_pcb = theta(1000.0, pcl_thetav, pbot);
		float vtmp_pct = theta(1000.0, pcl_thetav, ptop);

        // get the virtual temperature of the environment
		float vtmp_enb = interp_pressure(pbot, prof->pres, prof->vtmp, prof->NZ);
		float vtmp_ent = interp_pressure(ptop, prof->pres, prof->vtmp, prof->NZ);

		float buoy_bot = buoyancy(vtmp_pcb, vtmp_enb);
		float buoy_top = buoyancy(vtmp_pct, vtmp_ent);

		float hbot = interp_pressure(pbot, prof->pres, prof->hght, prof->NZ);
		float htop = interp_pressure(ptop, prof->pres, prof->hght, prof->NZ);

		float dz = htop - hbot;

        // integrate using trapezoid method 
		float lyre = ( (buoy_bot + buoy_top) / 2.0 ) * dz;
        if (lyre < 0.0) cinh += lyre;
	}
    return cinh;
}

void parcel_wobf(Profile* prof, Parcel* pcl) {
    lifter_wobus lifter;
    lift_parcel<lifter_wobus>(lifter, prof, pcl);
    return;
}

} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


