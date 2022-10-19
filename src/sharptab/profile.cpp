// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
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
		m_vtmp = new float[num_levs];
		m_uwin = new float[num_levs];
		m_vwin = new float[num_levs];
		m_omeg = new float[num_levs];

		m_nlevs = num_levs;
		m_snd_type = sounding_type;

	}



	/*
	 * Construct a new Profile who's arrays have already been allocated,
	 * passing the pointers of memory to be assigned to those arrays. 
	 * Includes the length of those arrays via num_levs, and the integer
	 * sounding type flag that describes what it is (observed, model, etc). 
	 */
	Profile::Profile(float* pres, float* hght, float* tmpc, float* dwpc, \
			float* mixr, float* vtmp, float* uwin, float* vwin,
			float* omeg, int num_levs, Source sounding_type) {

		m_pres = pres;
		m_hght = hght;
		m_tmpc = tmpc;
		m_dwpc = dwpc;
		m_mixr = mixr;
		m_vtmp = vtmp;
		m_uwin = uwin;
		m_vwin = vwin;
		m_omeg = omeg;

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
		delete[] m_vtmp;
		delete[] m_uwin;
		delete[] m_vwin;
		delete[] m_omeg;

	}

}
