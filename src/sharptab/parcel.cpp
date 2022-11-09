
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
#include "parcel.h"
#include "utils.h"
#include "thermo.h"

namespace sharp {


Parcel::Parcel() {
    pres = MISSING;
    tmpc = MISSING;
    dwpc = MISSING;
    
    lcl_pressure = MISSING;
    lfc_pressure = MISSING;
    eql_pressure = MISSING;
    mpl_pressure = MISSING;

    cape = 0.0;
    cinh = 0.0;
}


// Trying to mentally sketch out how I want parcel lifting
// and CAPE/CINH integration to work. I want it to be modular
// such that different methods of computing moist adiabats can
// be supported. I want to separate the parcel lifting from the
// actual numerical integration. I want to be able to store and
// reuse the temperature/pressure traces of parcel lifting. Need 
// to support a "fast CAPE/CINH" for effective inflow calculations.
// 
//
// 1) Create the parcel struct
// 2) Define its starting attributes (MU/ML/SFC/etc)
// 3) Lift the parcel/compute temperature trace
// 4) Set LCL/LFC/EL values
// 5) Integrate CAPE/CINH

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
    }
    else if (source == LPL::FCST) {

    }
    else if (source == LPL::MU) {
        _mu(prof, pcl);
    }
    else if (source == LPL::ML) {
        _ml(prof, pcl);
    }
    else if (source == LPL::EIL) {
    }
    else if (source == LPL::USR) {
        // do nothing - its already
        // been set!
    }
    else {
        // TO-DO: probably should raise an error or something
    }
}


void lift_parcel_wobus(const float* pressure, const float* height, 
                       int num_levs, Parcel* pcl) {

    float vtmp_pcl = virtual_temperature(pcl->pres, pcl->tmpc, pcl->dwpc);

    // int pcl_idx = 0;
    // pcl->temperature_trace[pcl_idx] = vtmp_pcl;
    // pcl->pressure_trace[pcl_idx] = pcl->pres;
    // pcl->height_trace[pcl_idx] = interp_pressure(pcl->pres, pressure, height, num_levs);
    // pcl_idx += 1;

    // Lift the parcel from the LPL to the LCL
    float pres_at_lcl; 
    float tmpc_at_lcl;

    // does not need to use virtual temperature, because this isn't
    // a quantity that has to do with buoyancy...
    drylift(pcl->pres, pcl->tmpc, pcl->dwpc, pres_at_lcl, tmpc_at_lcl);

    // vtmp_pcl = virtual_temperature(pres_at_lcl, tmpc_at_lcl, tmpc_at_lcl);
    // pcl->temperature_trace[pcl_idx] = vtmp_pcl;
    // pcl->pressure_trace[pcl_idx] = pres_at_lcl;
    // pcl->height_trace[pcl_idx] = interp_pressure(pcl->pres, pressure, height, num_levs);
    // pcl_idx += 1;

    // define the parcel saturated lift layer to be
    // from the LCL to the top of the profile available
    PressureLayer sat_layer(pres_at_lcl, pressure[num_levs-1]);

    // excludes the indices that would correspond to the exact top and
    // bottom of this layer - default is that bottom and top is interpolated
    LayerIndex sat_index = get_layer_index(sat_layer, pressure, num_levs);

    // iterate from the LCL to the top of the profile
    float pbot = pres_at_lcl;
    float tbot = tmpc_at_lcl;
    float ptop, ttop;
    for (int k = sat_index.kbot; k <= sat_index.ktop; k++) {
       ptop = pressure[k]; 
       ttop = wetlift(pbot, tbot, ptop);
       vtmp_pcl = virtual_temperature(ptop, ttop, ttop);
       // pcl->temperature_trace[pcl_idx] = vtmp_pcl;
       // pcl->pressure_trace[pcl_idx] = ptop;
       // pcl->height_trace[pcl_idx] = height[k];
       // pcl_idx += 1; 

       // set the top of the current layer to the
       // bottom of the next layer
       pbot = ptop;
       tbot = ttop;
    }

    // lift final level
    ttop = wetlift(pbot, tbot, sat_layer.ptop);
    vtmp_pcl = virtual_temperature(ptop, ttop, ttop);
   // pcl->temperature_trace[pcl_idx] = vtmp_pcl;
   // pcl->pressure_trace[pcl_idx] = sat_layer.ptop;
   // pcl->height_trace[pcl_idx] = height[sat_index.ktop+1]; 
   // pcl_idx += 1; 
}



} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


