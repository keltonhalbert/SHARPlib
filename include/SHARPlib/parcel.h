
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
#ifndef __SHARP_PARCEL_H__
#define __SHARP_PARCEL_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/profile.h>
#include <SHARPlib/thermo.h>


namespace sharp {

////////////    FUNCTORS    ///////////
//
/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A functor that calls the Wobus Wetlift funtion
 *
 * This functor is used to wrap the Wobus Wetlift function for parcel
 * lifting routines. These functions are wrapped by functors - classes
 * with ther operator() overloaded - so that functions can be
 * passed to templates in a way that the compiler can still
 * optimize, rather than using function pointers or lambdas.
 *
 * Specifically, this functor is designed to be passed as a template
 * argument to sharp::lift_parcel, so that the method of computing
 * moist adiabats can be changed without changing the overall parcel
 * lifting algorithm. The reason this is awesome is that the compiler
 * can still optimize and inline this code, while the user can configure
 * the parcel lifting algorithm to their specifications.
 *
 */
struct lifter_wobus {
    /**
     * \brief Overloads operator() to call sharp::wetlift.
     * \param pres      Initial percel pressure (hPa)
     * \param tmpk      Initial parcel temperature (degC)
     * \param new_pres  Final level of parcel after lift (hPa)
     */
    [[nodiscard]] float operator()(float pres, float tmpk,
                                   float new_pres) const noexcept {
        return wetlift(pres, tmpk, new_pres);
    }
};
//
////////////  END FUNCTORS   ///////////

/**
 * \brief Enum that defines the lifted parcel level (LPL) of origin.
 *
 * The SFC parcel is a surface-based parcel, where the parcel initial attributes
 * are set to the surface pressure, temperature, and dewpoint.
 *
 * The FCST parcel is a forecast-surface-based parcel, in which the afternoon
 * surface temperature and dewpoint are estimated and set as the parcel starting
 * values.
 *
 * The MU parcel is the most unstable parcel, in which the parcel attributes are
 * set to the pressure, temperature, and dewpoint of the maximum Theta-E level
 * within the bottom 300 hPa of the profile.
 *
 * The ML parcel is the mixed-layer parcel, in which the mean theta and water
 * vapor mixing ratio within the lowest 100 hPa are used to estimate a boundary
 * layer mean, and lifted from the surface.
 *
 * The EIL parcel is the mean Effective Inflow Layer parcel, in which a parcel
 * is lifted from the center of the Effective Inflow Layer
 *
 * The USR parcel means that the parcel initial lifting attributes have already
 * been set by the programmer or user, and there is no need for them to be
 * set or modified.
 */
enum class LPL : int {
    /**
     * \brief Surface Based Parcel
     */
    SFC = 1,

    /**
     * \brief Forecast Surface Parcel
     */
    FCST = 2,

    /**
     * \brief Most Unstable Parcel
     */
    MU = 3,  // most unstable

    /**
     * \brief 100mb Mixed Layer Parcel
     */
    ML = 4,  // 100mb mixed layer

    /**
     * \brief User-defined Parcel
     */
    USR = 5,  // user-defined

    /**
     * \brief Mean Effective Inflow Layer Parcel
     */
    EIL = 6,  // Mean effective inflow layer
	END,
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Data that defines a Parcel, its attributes, and derived quantities.
 *
 * Contains information about a Parcel's starting level and
 * thermodynamic attributes, as well as paramaters computed
 * using the parcel.
 */
struct Parcel {
    /**
     * \brief Parcel starting pressure (hPa)
     */
    float pres = MISSING;

    /**
     * \brief Parcel starting temperature (degC)
     */
    float tmpk = MISSING;

    /**
     * \brief Parcel starting dewpoint (degC)
     */
    float dwpk = MISSING;

    /**
     * \brief Pressure at the Lifted Condensation Level (hPa)
     */
    float lcl_pressure = MISSING;

    /**
     * \brief Pressure at the Level of Free Convection (hPa)
     */
    float lfc_pressure = MISSING;

    /**
     * \brief Pressure at the parcel Equilibrium Level (hPa)
     */
    float eql_pressure = MISSING;

    /**
     * \brief Parcel Convective Available Potential Energy (J/kg) between the
     * LFC and EL
     */
    float cape = 0.0;

    /**
     * \brief Parcel Convective Inhibition (J/kg) between the LFC and EL
     */
    float cinh = 0.0;

    /**
     * \brief The type of parcel this is
     */
    LPL source;

};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Sets the Lifted Parcel Level (LPL) attributes for a parcel.
 *
 * Before computing CAPE and CINH, the parcel's level of origin, or the
 * Lifted Parcel Level (LPL), must be set. sharp::LPL defines common
 * lifting levels, and passing the appropriate enum member
 * will set the parcel attributes to that type of parcel or level.
 *
 * If you wish to set a custom LPL, you can do so and then set the
 * source to sharp::LPL::USR.
 *
 * \param prof      sharp::Profile
 * \param pcl       sharp::Parcel
 * \param source    sharp::LPL
 */
void define_parcel(Profile* prof, Parcel* pcl, LPL source) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Lifts a parcel to compute buoyancy
 *
 * Lifts a parcel dry adiabatically from its sharp::LPL to its
 * LCL dry adiabatically, and then moist adiabatically from the
 * LCL to the top of the profile. This fills the buoyancy array
 * within a sharp::Profile, and that buoyancy array can be used
 * to find the LFC, EL, and integrate CAPE over various layers.
 *
 */
template <typename Lft>
void lift_parcel(Lft liftpcl, Profile* prof, Parcel* pcl) noexcept {
    // Lift the parcel from the LPL to the LCL
    float pres_lcl;
    float tmpk_lcl;
    drylift(pcl->pres, pcl->tmpk, pcl->dwpk, pres_lcl, tmpk_lcl);
    pcl->lcl_pressure = pres_lcl;

    const float thetav_lcl =
        theta(pres_lcl, virtual_temperature(pcl->pres, tmpk_lcl, tmpk_lcl),
              THETA_REF_PRESSURE);

    // Define the dry and saturated lift layers
    PressureLayer dry_lyr = {pcl->pres, pcl->lcl_pressure};
    PressureLayer sat_lyr = {pcl->lcl_pressure, prof->pres[prof->NZ - 1]};
    // The LayerIndex excludes the top and bottom for interpolation reasons
    const LayerIndex dry_idx = get_layer_index(dry_lyr, prof->pres, prof->NZ);
    const LayerIndex sat_idx = get_layer_index(sat_lyr, prof->pres, prof->NZ);

	// zero out any residual buoyancy from
	// other parcels that may have been lifted
	for (int k = 0; k < dry_idx.kbot; ++k) {
		prof->buoyancy[k] = 0.0;
	}

    // Fill the array with dry parcel buoyancy.
    // Virtual potential temperature (Theta-V)
    // is conserved for a parcels dry ascent to the LCL
    for (int k = dry_idx.kbot; k < dry_idx.ktop+1; ++k) {
        const float pcl_vtmp =
            theta(THETA_REF_PRESSURE, thetav_lcl, prof->pres[k]);
        prof->buoyancy[k] = buoyancy(pcl_vtmp, prof->vtmp[k]);
    }

    // fill the array with the moist parcel buoyancy
    for (int k = sat_idx.kbot; k < prof->NZ; ++k) {
        // compute above-lcl buoyancy here
        const float pcl_pres = prof->pres[k];
        const float pcl_tmpk = liftpcl(pres_lcl, tmpk_lcl, pcl_pres);
        // parcel is saturated, so temperature and dewpoint are same
        const float pcl_vtmp =
            virtual_temperature(pcl_pres, pcl_tmpk, pcl_tmpk);
        const float env_vtmp = prof->vtmp[k];
        const float buoy = buoyancy(pcl_vtmp, env_vtmp);
		prof->buoyancy[k] = buoy;
    }
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Find the LFC and EL that bounds the layer with the maximum CAPE
 *
 * Searches the buoyancy array for the LFC and EL layer combination that
 * results in the most CAPE in the given profile. The buoyancy array is
 * typically computed by calling sharp::lift_parcel. Once the LFC and EL
 * are found, the values are set in pcl->lfc_pres and pcl->eql_pres.
 *
 * \param pcl   a sharp::Parcel with its sharp::LPL/attributes defined
 * \param pres_arr  The pressure coordinate array
 * \param hght_arr  The height coordinate array
 * \param buoy_arr  The profile buoyancy array
 * \param NZ        The length of the arrays
 */
void find_lfc_el(Parcel* pcl, const float pres_arr[], const float hght_arr[],
                 const float buoy_arr[], const int N) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute CAPE and CINH for a previously lifted parcel.
 *
 * Assuming that sharp::lift_parcel has been called, cape_cinh
 * will integrate the area between the LFC and EL to compute CAPE,
 * and integrate the area between the LPL and LCL to compute CINH.
 *
 * The results are set in pcl->cape and pcl->cinh. 
 *
 * \param prof  A sharp::Profile of sounding data
 * \param pcl   A sharp::Parcel corresponding to the profile buoyancy array. 
 */
void cape_cinh(Profile* prof, Parcel *pcl) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Lifts a parcel using the Wobus method of computing moist adiabats
 *
 *
 * Lifts a parcel from its sharp::LPL to the Equilibrium Level (EL) using the
 * Wobus method for computing moist adiabats. See sharp::wobf and sharp::wetlift
 * for implementation details. Attributes such as CAPE, CINH, and specific
 * parcel levels are stored within the sharp::Parcel.
 *
 * \param prof      A sharp::Profile of atmospheric data
 * \param pcl       A sharp::Parcel with its sharp::LPL/attributes defined.
 */
void parcel_wobf(Profile* prof, Parcel* pcl) noexcept;


}  // end namespace sharp

namespace sharp::exper {}  // end namespace sharp::exper

#endif // __SHARP_PARCEL_H__
