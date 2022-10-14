// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
#pragma once

namespace sharp {
    
    struct Profile {

        float* pres;
        float* hght;
        float* tmpc;
        float* dwpc;
        float* mixr;
        float* vtmp;
        float* uwin;
        float* vwin;

        int nlevs;
        int type; // used to determine what kind of profile (e.g. observed, model)

		// Constructor that allocates memory for
		// the profile arrays
		Profile(int num_levels, int sounding_type) {
			this->pres = new float[num_levels];
			this->hght = new float[num_levels];
			this->tmpc = new float[num_levels];
			this->dwpc = new float[num_levels];
			this->mixr = new float[num_levels];
			this->vtmp = new float[num_levels];
			this->uwin = new float[num_levels];
			this->vwin = new float[num_levels];

			this->nlevs = num_levels;
			this->type = sounding_type;
			
		}


		// Constructor that takes previously allocated
		// arrays and assigns them to the struct
		Profile(float* pres, float* hght, float* tmpc, float* dwpc, \
				float* mixr, float* vtmp, float* uwin, float* vwin,
				int nlevs, int type) {

			this->pres = pres;
			this->hght = hght;
			this->tmpc = tmpc;
			this->dwpc = dwpc;
			this->mixr = mixr;
			this->vtmp = vtmp;
			this->uwin = uwin;
			this->vwin = vwin;
			this->nlevs = nlevs;
			this->type = type
		}

		// destructor that frees memory 
		// when we destroy the Profile
		~Profile() {

			delete[] this->pres;
			delete[] this->hght;
			delete[] this->tmpc;
			delete[] this->dwpc;
			delete[] this->mixr;
			delete[] this->vtmp;
			delete[] this->uwin;
			delete[] this->vwin;

		}
    }

}
