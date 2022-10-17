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
    
	/*
	 * The Profile struct is used to store the arrays
	 * of vertical data associated with either a weather
	 * balloon or forecast model sounding profile. In 
	 * C++, structs and classes are completely identical. 
	 * We use a struct in this case, however, since really
	 * it is just a data container, and we don't need
	 * things like inheritence, polymorphism, etc that are
	 * typically associated with classes. Additionally,
	 * in C++ structs have all members and functions 
	 * set to 'public' as default, which is desireable!
	 */
    struct Profile {

		/*
		 * The m_ convention is a C++ convention
		 * used to differenciate between variables
		 * that are the members of a class/struct, 
		 * and everything else. This is helpful for
		 * keeping track of what is a part of the 
		 * class, and what might be a local variable 
		 * or an argument passed to a function. 
		 */
        float* m_pres;
        float* m_hght;
        float* m_tmpc;
        float* m_dwpc;
        float* m_mixr;
        float* m_vtmp;
        float* m_uwin;
        float* m_vwin;

        int m_nlevs;
        int m_snd_type; // used to determine what kind of profile (e.g. observed, model)

		// Empty constructor
		Profile();

		// Constructor that allocates memory for
		// the profile arrays
		Profile(int num_levels, int sounding_type);

		// Constructor that takes previously allocated
		// arrays and assigns them to the struct
		Profile(float* pres, float* hght, float* tmpc, float* dwpc, \
				float* mixr, float* vtmp, float* uwin, float* vwin,
				int nlevs, int type);

		// Destructor that deallocates arrays
		~Profile();

    }

}
