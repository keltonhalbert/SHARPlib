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
 * with their operator() overloaded - so that functions can be
 * passed to templates in a way that the compiler can still
 * optimize, rather than using function pointers or lambdas.
 *
 * Specifically, this functor is designed to be passed as a template
 * argument to sharp::lift_parcel, so that the method of computing
 * moist adiabats can be changed without changing the overall parcel
 * lifting code. The reason this is awesome is that the compiler
 * can still optimize and inline this code, while the user can configure
 * the parcel lifting algorithm to their specifications.
 *
 */
struct lifter_wobus {
    /**
     * \brief Overloads operator() to call sharp::wetlift.
     * \param pres      Parcel pressure (Pa)
     * \param tmpk      Parcel temperature (degK)
     * \param new_pres  Final level of parcel after lift (Pa)
     *
     * \return          The virtual temperature of the lifted parcel 
     */
    [[nodiscard]] inline float operator()(float pres, float tmpk,
                                   float new_pres) const noexcept {
        float pcl_tmpk = wetlift(pres, tmpk, new_pres);
        return virtual_temperature(pcl_tmpk, mixratio(new_pres, pcl_tmpk));
    }
};

struct lifter_cm1 {
    /**
     * \brief The type of moist adiabat to use, as defined by sharp::adiabat
     */
    adiabat ma_type = adiabat::pseudo_liq;

    /**
     * \brief The pressure increment (Pa) to use for the iterative solver
     */
    float pressure_incr = 500.0f; 

    /**
     * \brief The iterative convergence criteria
     */
    float converge = 0.0002f;

    /**
     * \brief Used to keep track of mixing ratio for conserved/adiabatic lifting
     */
    float qv_total = MISSING;

    /**
     * \brief Water vapor mixing ratio variable updated during parcel lifts
     */
    float qv = 0.0;

    /**
     * \brief Liquid water mixing ratio variable updated during parcel lifts
     */
    float ql = 0.0;

    /**
     * \brief Ice water mixing ratio variable updated during parcel lifts
     */
    float qi = 0.0;

    /**
     * \brief Overloads operator() to call sharp::moist_adiabat_cm1
     *
     * \param pres      Parcel pressure (Pa)
     * \param tmpk      Parcel temperature (degK)
     * \param new_pres  Final level of parcel after lift (Pa)
     *
     * \return          The virtual temperature of the lifted parcel 
     */
    [[nodiscard]] inline float operator()(float pres, float tmpk,
                                   float new_pres) noexcept {
        float pcl_tmpk = moist_adiabat_cm1(
            pres, tmpk, new_pres, this->qv_total, this->qv, this->ql, this->qi,
            this->pressure_incr, this->converge, this->ma_type);
        return virtual_temperature(pcl_tmpk, this->qv, this->ql, this->qi);
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
 * within the bottom 400 hPa of the profile.
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
 * \brief Sets the Lifted Parcel Level (LPL) attributes for a sharp::Parcel.
 *
 * Before computing CAPE and CINH, the parcel's level of origin, or the
 * Lifted Parcel Level (LPL), must be set. The enum sharp::LPL defines
 * common lifting levels, and passing the appropriate enum member
 * will set the sharp::Parcel attributes to that type of parcel or level.
 *
 * If you wish to set a custom LPL, you can do so and then set the
 * source to sharp::LPL::USR.
 *
 * \param pressure      Array of pressure (Pa)
 * \param temperature   Array of temperature (degK)
 * \param dewpoint      Array of dewpoint temperature (degK)
 * \param wv_mixratio   Array of water vapor mixing ratio (kg/kg)
 * \param theta_arr     Array of potential temperature (degK)
 * \param thetae_arr    Array of eqiv. potential temperature (degK)
 * \param N             The length of the arrays
 * \param pcl           The sharp::Parcel to set the attributes to
 * \param LPL           The type of sharp::Parcel to define 
 */
void define_parcel(const float pressure[], const float temperature[],
                   const float dewpoint[], const float wv_mixratio[],
                   const float theta_arr[], const float thetae_arr[],
                   const int N, Parcel& pcl, LPL source) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Lifts a sharp::Parcel to compute buoyancy
 *
 * Lifts a sharp::Parcel dry adiabatically from its sharp::LPL to its
 * LCL dry adiabatically, and then moist adiabatically from the
 * LCL to the top of the profile. The moist adiabat used is determined
 * bu the type of lifting functor passed to the function (i.e. 
 * sharp::lifter_wobus or sharp::lifter_cm1). 
 *
 * \param liftpcl
 * \param pressure_arr              Array of env pressure (Pa)
 * \param virtual_temperature_arr   Array of env virtual temperature (degK)
 * \param buoyancy_arr              The array to fill with Buoyancy (m/s^2)
 * \param N                         The length of the arrays
 * \param pcl                       The sharp::Parcel to lift
 */
template <typename Lft>
void lift_parcel(Lft liftpcl, const float pressure_arr[],
                 const float virtual_temperature_arr[], float buoyancy_arr[],
                 const int N, Parcel* pcl) noexcept {
    // Lift the parcel from the LPL to the LCL
    float pres_lcl;
    float tmpk_lcl;
    drylift(pcl->pres, pcl->tmpk, pcl->dwpk, pres_lcl, tmpk_lcl);
    pcl->lcl_pressure = pres_lcl;

    const float qv_lcl = mixratio(pres_lcl, tmpk_lcl);
    const float thetav_lcl = theta(
        pres_lcl, virtual_temperature(tmpk_lcl, qv_lcl), THETA_REF_PRESSURE);

    // Define the dry and saturated lift layers
    PressureLayer dry_lyr = {pcl->pres, pcl->lcl_pressure};
    PressureLayer sat_lyr = {pcl->lcl_pressure, pressure_arr[N - 1]};
    // The LayerIndex excludes the top and bottom for interpolation reasons
    const LayerIndex dry_idx = get_layer_index(dry_lyr, pressure_arr, N);
    const LayerIndex sat_idx = get_layer_index(sat_lyr, pressure_arr, N);

	// zero out any residual buoyancy from
	// other parcels that may have been lifted
	for (int k = 0; k < dry_idx.kbot; ++k) {
		buoyancy_arr[k] = 0.0;
	}

    // Fill the array with dry parcel buoyancy.
    // Virtual potential temperature (Theta-V)
    // is conserved for a parcels dry ascent to the LCL
    for (int k = dry_idx.kbot; k < dry_idx.ktop+1; ++k) {
        const float pcl_vtmp =
            theta(THETA_REF_PRESSURE, thetav_lcl, pressure_arr[k]);
        buoyancy_arr[k] = buoyancy(pcl_vtmp, virtual_temperature_arr[k]);
    }

    // fill the array with the moist parcel buoyancy
    for (int k = sat_idx.kbot; k < N; ++k) {
        // compute above-lcl buoyancy here
        const float pcl_pres = pressure_arr[k];
        const float pcl_vtmpk = liftpcl(pres_lcl, tmpk_lcl, pcl_pres);
        const float env_vtmpk = virtual_temperature_arr[k];
        buoyancy_arr[k] = buoyancy(pcl_vtmpk, env_vtmpk);

    }
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Find the LFC and EL that bounds the layer with the maximum CAPE
 *
 * Searches the buoyancy array for the LFC and EL combination that
 * results in the most CAPE in the given profile. The buoyancy array is
 * typically computed by calling sharp::lift_parcel. Once the LFC and EL
 * are found, the values are set in sharp::Parcel.lfc_pres and 
 * sharp::Parcel.eql_pres.
 *
 * \param pcl       a sharp::Parcel with its sharp::LPL/attributes defined
 * \param pres_arr  The pressure coordinate array (Pa)
 * \param hght_arr  The height coordinate array (meters)
 * \param buoy_arr  The profile buoyancy array (m/s^2)
 * \param N        The length of the arrays
 */
void find_lfc_el(Parcel* pcl, const float pres_arr[], const float hght_arr[],
                 const float buoy_arr[], const int N) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute CAPE and CINH for a previously lifted sharp::Parcel.
 *
 * Assuming that sharp::lift_parcel has been called, cape_cinh
 * will integrate the area between the LFC and EL to compute CAPE,
 * and integrate the area between the LPL and LCL to compute CINH.
 *
 * The results are set in pcl->cape and pcl->cinh. 
 *
 * \param pres_arr  Array of pressure (Pa)
 * \param hght_arr  Array of height (meters)
 * \param buoy_arr  Array of buoyancy (m/s^2)
 * \param N         Length of arrays
 * \param pcl       A sharp::Parcel corresponding to the buoyancy array. 
 */
void cape_cinh(const float pres_arr[], const float hght_arr[],
               const float buoy_arr[], const int N, Parcel* pcl) noexcept;

}  // end namespace sharp


#endif // __SHARP_PARCEL_H__
