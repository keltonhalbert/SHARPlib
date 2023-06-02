/**
 * \file
 * \brief Data structures for containing data from vertical atmospheric<!--
 * --> sounding profiles
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

namespace sharp {

/*
 * Construct a new Profile who's arrays have a length
 * of num_levs, and an integer describing what kind
 * of sounding it is (oberved, model, etc).
 */
Profile::Profile(const int num_levs, Source sounding_type) : NZ(num_levs) {
    pres = new float[num_levs];
    hght = new float[num_levs];
    tmpk = new float[num_levs];
    dwpk = new float[num_levs];
    wspd = new float[num_levs];
    wdir = new float[num_levs];

    mixr = new float[num_levs];
    vtmp = new float[num_levs];
    uwin = new float[num_levs];
    vwin = new float[num_levs];
    vvel = new float[num_levs];
    theta = new float[num_levs];
    theta_e = new float[num_levs];
    moist_static_energy = new float[num_levs];
    buoyancy = new float[num_levs];

    //NZ = num_levs;
    snd_type = sounding_type;
}

/*
 * Handle the memory management of deallocating
 * arrays when we destroy a Profile.
 */
Profile::~Profile() {
    delete[] pres;
    delete[] hght;
    delete[] tmpk;
    delete[] dwpk;
    delete[] wspd;
    delete[] wdir;

    delete[] mixr;
    delete[] vtmp;
    delete[] uwin;
    delete[] vwin;
    delete[] vvel;
    delete[] theta;
    delete[] theta_e;
    delete[] moist_static_energy;
    delete[] buoyancy;
}

}  // end namespace sharp

