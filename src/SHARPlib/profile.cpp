/**
 * \file
 * \brief Data structures for containing data from vertical atmospheric sounding profiles 
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
#include <SHARPlib/profile.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>
#include <SHARPlib/utils.h>

namespace sharp {

/*
 * Construct a new Profile who's arrays have a length
 * of num_levs, and an integer describing what kind
 * of sounding it is (oberved, model, etc). 
 */
Profile::Profile(int num_levs, Source sounding_type) {
    pres = new float[num_levs];
    hght = new float[num_levs];
    tmpc = new float[num_levs];
    dwpc = new float[num_levs];
    wspd = new float[num_levs];
    wdir = new float[num_levs];

    mixr = new float[num_levs];
    relh = new float[num_levs];
    vtmp = new float[num_levs];
    uwin = new float[num_levs];
    vwin = new float[num_levs];
    vvel = new float[num_levs];
    theta = new float[num_levs];
    theta_e = new float[num_levs];

    NZ = num_levs;
    snd_type = sounding_type;
}

/*
 * Handle the memory management of deallocating
 * arrays when we destroy a Profile. 
 */
Profile::~Profile() {
    delete[] pres;
    delete[] hght;
    delete[] tmpc;
    delete[] dwpc;
    delete[] wspd;
    delete[] wdir;

    delete[] mixr;
    delete[] relh;
    delete[] vtmp;
    delete[] uwin;
    delete[] vwin;
    delete[] vvel;
    delete[] theta;
    delete[] theta_e;
}

Profile* create_profile(float *pres, float *hght, 
                        float *tmpc, float *dwpc,
                        float *wspd_or_u, float *wdir_or_v,
                        int NZ, Source sounding_type, bool windComponents) {

    Profile* prof = new Profile(NZ, sounding_type);
    
    for (int k = 0; k < NZ; k++) {
        prof->pres[k] = pres[k];
        prof->hght[k] = hght[k];
        prof->tmpc[k] = tmpc[k];
        prof->dwpc[k] = dwpc[k];

        float vtmp = virtual_temperature(pres[k], tmpc[k], dwpc[k]);
        float mixr = mixratio(pres[k], dwpc[k]);
        float thta = theta(pres[k], tmpc[k], 1000.0);
        float thte = thetae(pres[k], tmpc[k], dwpc[k]);

        prof->vtmp[k] = vtmp;
        prof->mixr[k] = mixr;
        prof->theta[k] = thta;
        prof->theta_e[k] = thte;


        if (windComponents) {
            prof->uwin[k] = wspd_or_u[k];
            prof->vwin[k] = wdir_or_v[k];

            WindVector vec = components_to_vector(prof->uwin[k], prof->vwin[k]);
            prof->wspd[k] = vec.speed;
            prof->wdir[k] = vec.direction;
        }
        else {
            prof->wspd[k] = wspd_or_u[k];
            prof->wdir[k] = wdir_or_v[k];

            WindComponents cmp = vector_to_components(prof->wspd[k], 
                                                      prof->wdir[k]);
            prof->uwin[k] = cmp.u;
            prof->vwin[k] = cmp.v;
        }
    }

    return prof;
}


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


