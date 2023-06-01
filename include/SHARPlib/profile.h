/**
 * \file
 * \brief Data structures for containing data from vertical<!--
 * --> atmospheric sounding profiles
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
#ifndef __SHARP_PROFILE_H__
#define __SHARP_PROFILE_H__

namespace sharp {

/**
 * \brief Enum used to differentiate sounding profile types
 *
 * Defines what type of soundings we can use
 * and helps differentiate those types.
 */
enum class Source : int {
    Observed = 0,
    PFC = 1,
    ACARS = 2,
	END,
};

/**
 * \brief Stores arrays of vertical atmospheric sounding profile data.
 *
 * The Profile struct is used to store the arrays
 * of vertical data associated with either a weather
 * balloon or forecast model sounding profile. In
 * C++, structs and classes are completely identical.
 * We use a struct in this case, however, since really
 * it is just a data container. Additionally,
 * in C++ structs have all members and functions
 * set to 'public' as default, which is desireable!
 */
struct Profile {
    /**
     * \brief Vertical array of pressure in Pa (descending)
     */
    float* pres;

    /**
     * \brief Vertical array of height in meters (ascending)
     */
    float* hght;

    /**
     * \brief Vertical array of temperature in degrees Kelvin
     */
    float* tmpk;

    /**
     * \brief Vertical array of dewpoint in degrees Kelvin
     */
    float* dwpk;

    /**
     * \brief Vertical array of water vapor mixing ratio in kg/kg
     */
    float* mixr;

    /**
     * \brief Vertical array of virtual temperature in degrees Kelvin
     */
    float* vtmp;

    /**
     * \brief Vertical array of wind speed in m/s
     */
    float* wspd;

    /**
     * \brief Vertical array of wind direction in degrees
     */
    float* wdir;

    /**
     * \brief Vertical array of the U wind component in m/s
     */
    float* uwin;

    /**
     * \brief Vertical array of the V wind component in m/s
     */
    float* vwin;

    /**
     * \brief Vertical array of vertical velocity in m/s
     */
    float* vvel;

    /**
     * \brief Vertical array of potential temperature in degrees Kelvin
     */
    float* theta;

    /**
     * \brief Vertical array of equiv. potential temperature in degrees Kelvin
     */
    float* theta_e;

    /**
     * \brief Vertical array of Moist Static Energy
     */
    float* moist_static_energy;

    /*
     * \brief Buoyancy for a lifted parcel (m/s^2)
     */
    float* buoyancy;

    /**
     * \brief The number of vertical levels present / the length of the arrays.
     */
    const int NZ;

    /**
     * \brief The type of soungind profile (e.g. Observed, PFC)
     */
    Source snd_type;

    /**
     * \brief Constructor that allocates arrays for our sounding
     */
    Profile(const int num_levels, Source sounding_type);

    Profile() = delete;

    // Destructor that deallocates arrays
    ~Profile();
};

}  // end namespace sharp

#endif // __SHARP_PROFILE_H__
