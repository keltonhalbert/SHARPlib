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
#pragma once

namespace sharp {
    
/**
 * \brief Stores arrays of vertical atmospheric sounding profile data. 
 * \author
 *   Kelton Halbert                 \n
 *   Email: kelton.halbert@noaa.gov \n
 *   License: Apache 2.0            \n
 * \date 2022-10-18
 *
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
 * \n\n
 * The m_ convention is a C++ convention
 * used to differenciate between variables
 * that are the members of a class/struct, 
 * and everything else. This is helpful for
 * keeping track of what is a part of the 
 * class, and what might be a local variable 
 * or an argument passed to a function. 
 */
struct Profile {

    /**
     * \brief Enum used to differentiate sounding profile types
     *
     * Defines what type of soundings we can use
     * and helps differentiate those types.
     */
    enum Source : unsigned int {
        Observed = 0,
        PFC = 1,
        ACARS = 2,
    };

    /**
     * \brief Vertical array of pressure in millibars (descending)
     */
    float* m_pres;
    /**
     * \brief Vertical array of height in meters (ascending)
     */
    float* m_hght;
    /**
     * \brief Vertical array of temperature in degrees Celsius 
     */
    float* m_tmpc;
    /**
     * \brief Vertical array of dewpoint in degrees Celsius 
     */
    float* m_dwpc;
    /**
     * \brief Vertical array of water vapor mixing ratio in degrees Celsius 
     */
    float* m_mixr;
    /**
     * \brief Vertical array of virtual temperature in degrees Celsius 
     */
    float* m_vtmp;
    /**
     * \brief Vertical array of the U wind component in knots 
     */
    float* m_uwin;
    /**
     * \brief Vertical array of the V wind component in knots 
     */
    float* m_vwin;
    /**
     * \brief Vertical array of omega (vertical velocity on pressure surfaces) in millibars/hour 
     */
    float* m_omeg;



    /**
     * \brief The number of vertical levels present / the length of the arrays. 
     */
    int m_nlevs;

    /**
     * \brief The type of soungind profile (e.g. Observed, PFC) 
     */
    Source m_snd_type; 



    // Constructor that allocates memory for
    // the profile arrays
    Profile(int num_levels, Source sounding_type);

    // Constructor that takes previously allocated
    // arrays and assigns them to the struct
    Profile(float* pres, float* hght, float* tmpc, float* dwpc, \
            float* mixr, float* vtmp, float* uwin, float* vwin,
            float* omeg, int nlevs, Source sounding_type);

    // Destructor that deallocates arrays
    ~Profile();

};


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


