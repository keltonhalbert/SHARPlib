/**
 * \file
 * \brief Routines used for parcel lifting and integration
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
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
#include <iostream>

#include <cstddef>

namespace sharp {

////////////    FUNCTORS    ///////////
//

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A functor that calls the Wobus Wetlift funtion
 *
 * This functor is used to wrap the Wobus Wetlift function for parcel
 * lifting routines. Functors - classes with their operator()
 * overloaded - are used so that functions can be
 * passed to templates in a way that the compiler can still
 * optimize, rather than using function pointers or lambdas.
 *
 * Specifically, this functor is designed to be passed as a template
 * argument to sharp::Parcel::lift_parcel, so that the method of computing
 * moist adiabats can be changed without changing the overall parcel
 * lifting code. The reason this is awesome is that the compiler
 * can still optimize and inline this code, while the user can configure
 * the parcel lifting algorithm to their specifications.
 */
struct lifter_wobus {
    static constexpr bool lift_from_lcl = true;

    /**
     * \brief Overloads operator() to call sharp::wetlift.
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(const float pres, const float tmpk,
                                          const float new_pres) const {
        return wetlift(pres, tmpk, new_pres);
    }

    /**
     * \brief Computes the virtual temperature of the parcel
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (degK)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature(const float pres, const float tmpk) const {
        return virtual_temperature(tmpk, mixratio(pres, tmpk));
    }
};

struct lifter_cm1 {
    static constexpr bool lift_from_lcl = false;

    /**
     * \brief The type of moist adiabat to use, as defined by sharp::adiabat
     */
    adiabat ma_type = adiabat::pseudo_liq;

    /**
     * \brief The pressure increment (Pa) to use for the iterative solver
     */
    float pressure_incr = 500.0f;

    /**
     * \brief The iterative convergence criteria (K)
     */
    float converge = 0.001f;

    /**
     * \brief Used to keep track of mixing ratio for conserved/adiabatic lifting
     */
    float rv_total = MISSING;

    /**
     * \brief Water vapor mixing ratio variable updated during parcel lifts
     */
    float rv = 0.0;

    /**
     * \brief Liquid water mixing ratio variable updated during parcel lifts
     */
    float rl = 0.0;

    /**
     * \brief Ice water mixing ratio variable updated during parcel lifts
     */
    float ri = 0.0;

    /**
     * \brief Overloads operator() to call sharp::moist_adiabat_cm1
     *
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The virtual temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(const float pres, const float tmpk,
                                          const float new_pres) {
        return moist_adiabat_cm1(
            pres, tmpk, new_pres, this->rv_total, this->rv, this->rl, this->ri,
            this->pressure_incr, this->converge, this->ma_type);
    }

    /**
     * \brief Computes the virtual temperature of the parcel
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (degK)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature([[maybe_unused]] const float pres, const float tmpk) const {
        return virtual_temperature(tmpk, this->rv, this->rl, this->ri);
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
     * \brief Mixed Layer Parcel
     */
    ML = 4,  // Mixed layer

    /**
     * \brief User-defined Parcel
     */
    USR = 5,  // user-defined
    END,
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Data that defines a Parcel, its attributes, and derived quantities.
 *
 * Contains information about a Parcel's starting level and
 * thermodynamic attributes, as well as paramaters computed
 * using the parcel.
 */
struct Parcel {
    /**
     * \brief Parcel starting pressure (Pa)
     */
    float pres = MISSING;

    /**
     * \brief Parcel starting temperature (K)
     */
    float tmpk = MISSING;

    /**
     * \brief Parcel starting dewpoint (K)
     */
    float dwpk = MISSING;

    /**
     * \brief Pressure at the Lifted Condensation Level (Pa)
     */
    float lcl_pressure = MISSING;

    /**
     * \brief Pressure at the Level of Free Convection (Pa)
     */
    float lfc_pressure = MISSING;

    /**
     * \brief Pressure at the parcel Equilibrium Level (Pa)
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

    Parcel();
    Parcel(float pressure, float temperature, float dewpoint, LPL lpl);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Lifts a sharp::Parcel to compute buoyancy
     *
     * Lifts a sharp::Parcel dry adiabatically from its sharp::LPL to its
     * LCL dry adiabatically, and then moist adiabatically from the
     * LCL to the top of the profile. The moist adiabat used is determined
     * bu the type of lifting functor passed to the function (i.e.
     * sharp::lifter_wobus or sharp::lifter_cm1).
     *
     * \param   liftpcl         Parcel lifting function/functor
     * \param   pressure_arr    Array of env pressure (Pa)
     * \param   virtemp_arr     Array of env virtual temperature (K)
     * \param   buoyancy_arr    The array to fill with Buoyancy (m/s^2)
     * \param   N               The length of the arrays
     * \param   pcl_vtmpk_arr   Optional array to fill virtual temp (K)
     */
    template <typename Lft>
    void lift_parcel(Lft& liftpcl, const float pressure_arr[],
                     const float virtemp_arr[], float buoyancy_arr[],
                     const std::ptrdiff_t N, float pcl_vtmpk_arr[] = nullptr) {
        // Lift the parcel from the LPL to the LCL
        float pres_lcl;
        float tmpk_lcl;
        drylift(this->pres, this->tmpk, this->dwpk, pres_lcl, tmpk_lcl);
        // If we are lifting elevated parcel (i.e. EIL), we need to make
        // sure our LCL isnt above the top of our data.
        if (pres_lcl < pressure_arr[N - 1]) return;

        this->lcl_pressure = pres_lcl;

        const float qv_lcl = mixratio(pres_lcl, tmpk_lcl);
        const float thetav_lcl =
            theta(pres_lcl, virtual_temperature(tmpk_lcl, qv_lcl),
                  THETA_REF_PRESSURE);

        // Define the dry and saturated lift layers
        PressureLayer dry_lyr = {this->pres, this->lcl_pressure};
        PressureLayer sat_lyr = {this->lcl_pressure, pressure_arr[N - 1]};

        // The LayerIndex excludes the top and bottom for interpolation reasons
        const LayerIndex dry_idx = get_layer_index(dry_lyr, pressure_arr, N);
        const LayerIndex sat_idx = get_layer_index(sat_lyr, pressure_arr, N);

        // zero out any residual buoyancy from
        // other parcels that may have been lifted
        for (std::ptrdiff_t k = 0; k < dry_idx.kbot; ++k) {
            buoyancy_arr[k] = 0.0f;
            if (pcl_vtmpk_arr) pcl_vtmpk_arr[k] = 0.0f;
        }

        // Fill the array with dry parcel buoyancy.
        // Virtual potential temperature (Theta-V)
        // is conserved for a parcels dry ascent to the LCL
        for (std::ptrdiff_t k = dry_idx.kbot; k < dry_idx.ktop + 1; ++k) {
            const float pcl_vtmp =
                theta(THETA_REF_PRESSURE, thetav_lcl, pressure_arr[k]);
            buoyancy_arr[k] = buoyancy(pcl_vtmp, virtemp_arr[k]);
            if (pcl_vtmpk_arr) pcl_vtmpk_arr[k] = pcl_vtmp;
        }

        float pres_bot = pres_lcl;
        float tmpk_bot = tmpk_lcl;

        // fill the array with the moist parcel buoyancy
        for (std::ptrdiff_t k = sat_idx.kbot; k < N; ++k) {
            // compute above-lcl buoyancy here
            const float pcl_pres = pressure_arr[k];
            const float pcl_tmpk = liftpcl(pres_bot, tmpk_bot, pcl_pres);

            if constexpr (!Lft::lift_from_lcl) {
                pres_bot = pcl_pres;
                tmpk_bot = pcl_tmpk;
            }

            const float pcl_vtmpk = liftpcl.parcel_virtual_temperature(pcl_pres, pcl_tmpk);
            const float env_vtmpk = virtemp_arr[k];
            buoyancy_arr[k] = buoyancy(pcl_vtmpk, env_vtmpk);
            if (pcl_vtmpk_arr) pcl_vtmpk_arr[k] = pcl_vtmpk;
        }
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Find the LFC and EL that bounds the layer with the maximum CAPE
     *
     * Searches the buoyancy array for the LFC and EL combination that
     * results in the most CAPE in the given profile. The buoyancy array is
     * typically computed by calling sharp::lift_parcel. Once the LFC and EL
     * are found, the values are set in sharp::Parcel.lfc_pres and
     * sharp::Parcel.eql_pres.
     *
     * \param   pres_arr    The pressure coordinate array (Pa)
     * \param   hght_arr    The height coordinate array (meters)
     * \param   buoy_arr    The profile buoyancy array (m/s^2)
     * \param   N           The length of the arrays
     */
    void find_lfc_el(const float pres_arr[], const float hght_arr[],
                     const float buoy_arr[], const std::ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Compute CAPE and CINH for a previously lifted sharp::Parcel.
     *
     * Assuming that sharp::lift_parcel has been called, cape_cinh
     * will integrate the area between the LFC and EL to compute CAPE,
     * and integrate the area between the LPL and LCL to compute CINH.
     *
     * The results are set in pcl->cape and pcl->cinh.
     *
     * \param   pres_arr    Array of pressure (Pa)
     * \param   hght_arr    Array of height (meters)
     * \param   buoy_arr    Array of buoyancy (m/s^2)
     * \param   N           Length of arrays
     */
    void cape_cinh(const float pres_arr[], const float hght_arr[],
                   const float buoy_arr[], const std::ptrdiff_t N);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Construct and return a surface-based Parcel
     *
     * Given input values of surface pressure, temperature,
     * and dewpoint temperature, construct and return a Surface
     * Based Parcel.
     *
     * \param    pressure        Surface pressure (Pa)
     * \param    temperature     Surface temperature (K)
     * \param    dewpoint        Surface dewpoint (K)
     * \return   sharp::Parcel with surface values
     */
    static inline Parcel surface_parcel(const float pressure,
                                        const float temperature,
                                        const float dewpoint) noexcept {
        return Parcel(pressure, temperature, dewpoint, LPL::SFC);
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Construct and return a mixed-layer Parcel
     *
     * Given input arrays of pressure, height, potential
     * temperature, and water vapor mixing ratio, as well
     * as a defined sharp::PressureLayer or
     * sharp::HeightLayer, compute and return a mixed-layer
     * Parcel.
     *
     * \param   pressure            Array of pressure (Pa)
     * \param   height              Array of height (meters)
     * \param   pot_temperature     Array of potential temperature (K)
     * \param   wv_mixratio         Array of water vapor mixing ratio (unitless)
     * \param   N                   Length of arrays
     * \param   Lyr                 sharp::PressureLayer or sharp::HeightLayer
     * \return sharp::Parcel with mixed-layer values
     */
    template <typename Lyr>
    static Parcel mixed_layer_parcel(Lyr& mix_layer, const float pressure[],
                                     const float height[],
                                     const float pot_temperature[],
                                     const float wv_mixratio[],
                                     const std::ptrdiff_t N) noexcept {
        float mean_mixr, mean_thta, pcl_pres;
        if constexpr (Lyr::coord == LayerCoordinate::pressure) {
            mean_mixr = layer_mean(mix_layer, pressure, wv_mixratio, N);
            mean_thta = layer_mean(mix_layer, pressure, pot_temperature, N);
            pcl_pres = mix_layer.bottom;

        } else {
            mean_mixr = layer_mean(mix_layer, height, pressure, wv_mixratio, N);
            mean_thta =
                layer_mean(mix_layer, height, pressure, pot_temperature, N);
            pcl_pres =
                sharp::interp_height(mix_layer.bottom, height, pressure, N);
        }
        const float tmpk = theta(THETA_REF_PRESSURE, mean_thta, pcl_pres);
        const float dwpk = temperature_at_mixratio(mean_mixr, pcl_pres);

        return Parcel(pcl_pres, tmpk, dwpk, LPL::ML);
    }

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Construct and return a most-unstable Parcel
     *
     * Given input arrays of pressure, height, temperature,
     * virtual temperature, and dewpoint, a scratch/working
     * array to store values of buoyancy, and a defined
     * sharp::PressureLayer or sharp::HeightLayer to search
     * over, find and return the most-unstable parcel.
     *
     * \param   pressure            Array of pressure (Pa)
     * \param   height              Array of height (meters)
     * \param   temperature         Array of temperature (K)
     * \param   virtemp             Array of virtual temperature (K)
     * \param   dewpoint            Array of dewpoint temperature (K)
     * \param   buoyancy            Writeable array for buoyancy calcs (m/s^2)
     * \param   N                   Length of arrays
     * \param   search_layer        sharp::PressureLayer or sharp::HeightLay
     * \param   lifter              The parcel moist adiabatic ascent function
     * \return The most unstable sharp::Parcel within the search layer
     */
    template <typename Lyr, typename Lft>
    static Parcel most_unstable_parcel(Lyr& search_layer, Lft& lifter,
                                       const float pressure[],
                                       const float height[],
                                       const float temperature[],
                                       const float virtemp[],
                                       const float dewpoint[], float buoyancy[],
                                       const std::ptrdiff_t N) noexcept {
        LayerIndex lyr_idx;
        if constexpr (Lyr::coord == LayerCoordinate::pressure) {
            lyr_idx = get_layer_index(search_layer, pressure, N);

        } else {
            lyr_idx = get_layer_index(search_layer, height, N);
        }

        Parcel max_parcel;
        for (std::ptrdiff_t idx = lyr_idx.kbot; idx <= lyr_idx.ktop; ++idx) {
            const float pres = pressure[idx];
            const float tmpk = temperature[idx];
            const float dwpk = dewpoint[idx];
            Parcel pcl(pres, tmpk, dwpk, LPL::MU);

            pcl.lift_parcel(lifter, pressure, virtemp, buoyancy, N);
            pcl.cape_cinh(pressure, height, buoyancy, N);
            if (pcl.cape > max_parcel.cape) max_parcel = pcl;
        }

        return max_parcel;
    }
};

}  // end namespace sharp

#endif  // __SHARP_PARCEL_H__
