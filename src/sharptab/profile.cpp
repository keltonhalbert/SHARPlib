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
#include "profile.h"

namespace sharp {

/*
 * Construct a new Profile who's arrays have a length
 * of num_levs, and an integer describing what kind
 * of sounding it is (oberved, model, etc). 
 */
Profile::Profile(int num_levs, Source sounding_type) {

    m_pres = new float[num_levs];
    m_hght = new float[num_levs];
    m_tmpc = new float[num_levs];
    m_dwpc = new float[num_levs];
    m_mixr = new float[num_levs];
    m_relh = new float[num_levs];
    m_vtmp = new float[num_levs];
    m_wspd = new float[num_levs];
    m_wdir = new float[num_levs];
    m_uwin = new float[num_levs];
    m_vwin = new float[num_levs];
    m_vvel = new float[num_levs];

    m_theta = new float[num_levs];
    m_theta_e = new float[num_levs];

    m_nlevs = num_levs;
    m_snd_type = sounding_type;

}



/*
 * Handle the memory management of deallocating
 * arrays when we destroy a Profile. 
 */
Profile::~Profile() {

    delete[] m_pres;
    delete[] m_hght;
    delete[] m_tmpc;
    delete[] m_dwpc;
    delete[] m_mixr;
    delete[] m_relh;
    delete[] m_vtmp;
    delete[] m_wspd;
    delete[] m_wdir;
    delete[] m_uwin;
    delete[] m_vwin;
    delete[] m_vvel;
    delete[] m_theta;
    delete[] m_theta_e;

}

} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


