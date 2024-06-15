/**
 * \file
 * \brief Routines used for parcel lifting and integration
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-11-09
 *
 * \contributor
 *   Amelia Urquhart                 \n
 *   Email: amelia.r.h.urquhart-1@ou.edu\n
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
#include <SHARPlib/winds.h>

namespace sharp {

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
 * \param   pressure        Array of pressure (Pa)
 * \param   temperature     Array of temperature (degK)
 * \param   dewpoint        Array of dewpoint temperature (degK)
 * \param   wv_mixratio     Array of water vapor mixing ratio (kg/kg)
 * \param   theta_arr       Array of potential temperature (degK)
 * \param   thetae_arr      Array of eqiv. potential temperature (degK)
 * \param   N               The length of the arrays
 * \param   pcl             The sharp::Parcel to set the attributes to
 * \param   source          The type of sharp::Parcel to define
 */
void define_parcel(const float pressure[], const float temperature[],
                   const float dewpoint[], const float wv_mixratio[],
                   const float theta_arr[], const float thetae_arr[],
                   const int N, Parcel& pcl, LPL source);

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
 * \param   liftpcl
 * \param   pressure_arr                Array of env pressure (Pa)
 * \param   virtual_temperature_arr     Array of env virtual temperature (degK)
 * \param   buoyancy_arr                The array to fill with Buoyancy (m/s^2)
 * \param   N                           The length of the arrays
 * \param   pcl                         The sharp::Parcel to lift
 */
template <typename Lft>
void lift_parcel(Lft liftpcl, const float pressure_arr[],
                 const float virtual_temperature_arr[], float buoyancy_arr[],
                 const int N, Parcel* pcl) {
    // Lift the parcel from the LPL to the LCL
    float pres_lcl;
    float tmpk_lcl;
    drylift(pcl->pres, pcl->tmpk, pcl->dwpk, pres_lcl, tmpk_lcl);
    // If we are lifting elevated parcel (i.e. EIL), we need to make
    // sure out LCL isnt above the top of our data.
    if (pres_lcl < pressure_arr[N - 1]) return;

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
    for (int k = dry_idx.kbot; k < dry_idx.ktop + 1; ++k) {
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
 */
struct lifter_wobus {
    /**
     * \brief Overloads operator() to call sharp::wetlift.
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (degK)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The virtual temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(float pres, float tmpk,
                                          float new_pres) const {
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
     * \param   tmpk        Parcel temperature (degK)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The virtual temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(float pres, float tmpk,
                                          float new_pres) {
        float pcl_tmpk = moist_adiabat_cm1(
            pres, tmpk, new_pres, this->rv_total, this->rv, this->rl, this->ri,
            this->pressure_incr, this->converge, this->ma_type);
        return virtual_temperature(pcl_tmpk, this->rv, this->rl, this->ri);
    }
};

/**
 * \author Amelia Urquhart - OU-SoM
 *
 * \brief Lifts a sharp::Parcel according to the lifting formulas from Peters et al. 2022
 */
struct lifter_peters_et_al {
    /**
     * \brief The type of moist adiabat to use, as defined by sharp::ascent_type
     */
    ascent_type ma_type = ascent_type::adiab_entr;

    /**
     * \brief The entrainment rate used by the lifter
     */
    float entr_rate = MISSING;

    /**
     * \brief The pressure increment (Pa) to use for the iterative solver
     */
    float pressure_incr = 500.0f;

    /**
     * \brief The iterative convergence criteria
     */
    float converge = 0.001f;

    /**
     * \brief The warm limit of the mixed phase range of temperatures (Kelvin)
     */
    float mixed_phase_warm = 273.15;

    /**
     * \brief The cold limit of the mixed phase range of temperatures (Kelvin)
     */
    float mixed_phase_cold = 253.15; 

    /**
     * \brief Internal temperature variable updated during parcel lifts [WIP]
     * 
     * Keeps track of the last temperature lifted to for efficiency's sake
     */
    float temperature = MISSING;

    /**
     * \brief Internal pressure variable updated during parcel lifts [WIP]
     * 
     * Keeps track of the last temperature lifted to for efficiency's sake
     */
    float pressure = MISSING;

    /**
     * \brief Water vapor mass fraction variable updated during parcel lifts [WIP]
     */
    float qv = MISSING;

    /**
     * \brief Total vapor (water vapor + cloud condensate) mass fraction variable updated during parcel lifts [WIP]
     */
    float qt = MISSING;

    /**
     * \brief Environmental profile used by the lifter
     */
    Profile* profile;

    /**
     * \brief Overloads operator() to call sharp::moist_adiabat_peters_et_al
     *
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (degK)
     * \param   new_pres    Final level of parcel after lift (Pa)
     *
     * \return  The density temperature of the lifted parcel
     */
    [[nodiscard]] inline float operator()(float pres, float tmpk,
                                          float new_pres) {

        bool can_do_efficient_lift = true;

        // Checks whether lifting can be continued from the last time operator 
        // was called. Massively increases efficiency if profile is sorted in
        // order of increasing height.
        if(new_pres > pressure) {
            can_do_efficient_lift = false;
        }

        float pcl_tmpk;

        // TODO!!! Write actual moist adiabat solver code.
        // THIS WILL BE PAINFUL.
        if(can_do_efficient_lift) {
            pcl_tmpk = moist_adiabat_peters_et_al(pressure, temperature, new_pres, 
                this->qv, this->qt, this->profile, this->pressure_incr, 
                this->entr_rate, this->ma_type);
        } else {
            pcl_tmpk = moist_adiabat_peters_et_al(pres, tmpk, new_pres, 
                this->qv, this->qt, this->profile, this->pressure_incr, 
                this->entr_rate, this->ma_type);
        }

        pressure = new_pres;
        temperature = pcl_tmpk;

        return density_temperature(pcl_tmpk, qv, qt);
    }

    /**
     * If using entrainment, this must be run to give the lifter an 
     * environmental profile. If not using entrainment, this can be ignored.
     */
    void set_profile(Profile* prof) {
        profile = prof;
    }

    /**
     * If using entrainment, this must be run to determine the entrainment rate
     */
    void determine_entrainment_rate(Profile* prof, LPL lpl) {
    // Written at 4 AM, will definitely need to be double-checked
        Parcel pcl;

        define_parcel(prof->pres, prof->tmpk, prof->dwpk, prof->mixr,
                             prof->theta, prof->theta_e, prof->NZ, pcl,
                             lpl);

        lifter_peters_et_al lifter;

        if(ma_type == ascent_type::adiab_entr || ma_type == ascent_type::adiab_nonentr) {
            lifter.ma_type = ascent_type::adiab_nonentr;
        }

        if(ma_type == ascent_type::pseudo_entr || ma_type == ascent_type::pseudo_nonentr) {
            lifter.ma_type = ascent_type::pseudo_nonentr;
        }

        lift_parcel(lifter, prof->pres, prof->vtmp, prof->buoyancy,
                           prof->NZ, &pcl);

        entr_rate = entrainment_rate(prof->pres, prof->hght, prof->tmpk, prof->moist_static_energy, prof->uwin, prof->vwin, prof->NZ, &pcl);
    }

    // Has to be copied here outside of convective.h to avoid circular include
    float entrainment_rate(const float pressure[], const float height[],
                        const float temperature[], const float mse_arr[],
                        const float u_wind[], const float v_wind[], const int N,
                        Parcel *pcl) {
        // if cape is zero, we get a divide by zero issue later.
        // there can "technically" be LFC/EL without CAPE because of very,
        // very shallow buoyancy near zero when searching for LFC/EL.
        // To-Do: maybe refactor LFC/EL search?
        if ((pcl->lfc_pressure == MISSING) || (pcl->lfc_pressure == MISSING) ||
            (pcl->cape == 0)) {
            return 0.0;
        }

        float *mse_diff = new float[N];
        float *mse_bar = new float[N];

        // compute MSE_bar
        mse_bar[0] = mse_arr[0];
        const float psfc = pressure[0];
        for (int k = 1; k < N; ++k) {
            PressureLayer mn_lyr = {psfc, pressure[k]};
            mse_bar[k] = layer_mean(mn_lyr, pressure, mse_arr, N);
        }

        // compute MSE_star
        const float hsfc = height[0];
        for (int k = 0; k < N; ++k) {
            const float tmpk = temperature[k];
            const float rsat = mixratio(pressure[k], tmpk);
            const float qsat = specific_humidity(rsat);
            const float height_agl = height[k] - hsfc;
            const float mse_star = moist_static_energy(height_agl, tmpk, qsat);
            mse_diff[k] = buoyancy_dilution_potential(tmpk, mse_bar[k], mse_star);
        }

        // compute NCAPE
        PressureLayer plyr = {pcl->lfc_pressure, pcl->eql_pressure};
        HeightLayer hlyr = pressure_layer_to_height(plyr, pressure, height, N);
        const float NCAPE = integrate_layer_trapz(hlyr, mse_diff, height, N);

        // Compute the Bunkers non-parcel based storm motion
        WindComponents strm_mtn =
            storm_motion_bunkers(pressure, height, u_wind, v_wind, N, {0.0, 6000.0},
                                {0.0, 6000.0}, false, false);

        // get the mean 0-1km storm relative wind
        HeightLayer layer = {hsfc + 0.0f, hsfc + 1000.0f};
        LayerIndex layer_idx = get_layer_index(layer, height, N);

        // loop from the surface to the last level before 1km AGL.
        float V_sr_mean = 0.0;
        int count = 0;
        for (int k = 0; k < layer_idx.ktop + 1; ++k) {
    #ifndef NO_QC
            if ((u_wind[k] == MISSING) || (v_wind[k] == MISSING)) {
                continue;
            }
    #endif
            V_sr_mean +=
                vector_magnitude(u_wind[k] - strm_mtn.u, v_wind[k] - strm_mtn.v);
            count += 1;
        }
        const float u_1km = interp_height(layer.top, height, u_wind, N);
        const float v_1km = interp_height(layer.top, height, v_wind, N);

        V_sr_mean += vector_magnitude(u_1km - strm_mtn.u, v_1km - strm_mtn.v);
        V_sr_mean = V_sr_mean / (count + 1);

        // now for all of the ECAPE nonsense
        constexpr float L = 120.0;
        const float H =
            interp_pressure(pcl->eql_pressure, pressure, height, N) - hsfc;
        constexpr float sigma = 1.6;
        constexpr float alpha = 0.8;
        const float pitchfork = (VKSQ * (alpha * alpha) * (PI * PI) * L) /
                                (PRANDTL * (sigma * sigma) * H);
        const float V_sr_tilde = V_sr_mean / std::sqrt(2.0f * pcl->cape);
        const float V_sr_tilde_sq = V_sr_tilde * V_sr_tilde;
        const float N_tilde = NCAPE / pcl->cape;

        const float term1 = pitchfork / V_sr_tilde_sq;
        const float term2 = 1.0f + pitchfork + term1 * N_tilde;
        const float term3 = 4.0f * term1 * (1.0f - pitchfork * N_tilde);
        const float sqrt_term = term2 * term2 + term3;

        // in the case of a negative solution,
        // set ECAPE to 0;
        if (sqrt_term < 0) {
            return 0;
        }

        const float E_tilde = V_sr_tilde_sq + (-1.0f - pitchfork - (term1)*N_tilde +
                                            std::sqrt(sqrt_term)) /
                                                (2.0f * term1);
        delete[] mse_diff;
        delete[] mse_bar;

        float E_tilde_no_vsr = E_tilde - V_sr_tilde_sq;

        float entrainment_rate = ((2 * (1 - E_tilde_no_vsr)) / (E_tilde_no_vsr + N_tilde)) / H;

        return entrainment_rate;
    }

    // Has to be copied here outside of convective.h to avoid circular include
    WindComponents storm_motion_bunkers(const float pressure[],
                                        const float height[], const float u_wind[],
                                        const float v_wind[], const int N,
                                        HeightLayer mean_wind_layer_agl,
                                        HeightLayer wind_shear_layer_agl,
                                        const bool leftMover = false,
                                        const bool pressureWeighted = false) {
        constexpr float deviation = 7.5;  // deviation from mean wind in m/s

        PressureLayer mw_lyr = height_layer_to_pressure(mean_wind_layer_agl,
                                                        pressure, height, N, true);

        WindComponents layer_mean_wind = {MISSING, MISSING};
        layer_mean_wind =
            mean_wind(mw_lyr, pressure, u_wind, v_wind, N, pressureWeighted);

        // The shear is computed by finding the 500m deep
        // mean winds at the top and bottom of the wind_shear_layer
        // and then differencing the two. Means are not pressure weighted
        // for the non-parcel based method.
        HeightLayer h_layer_lo = {wind_shear_layer_agl.bottom,
                                wind_shear_layer_agl.bottom + 500.0f};
        HeightLayer h_layer_hi = {wind_shear_layer_agl.top - 500.0f,
                                wind_shear_layer_agl.top};

        PressureLayer p_layer_lo =
            height_layer_to_pressure(h_layer_lo, pressure, height, N, true);
        PressureLayer p_layer_hi =
            height_layer_to_pressure(h_layer_hi, pressure, height, N, true);

        WindComponents winds_lo =
            mean_wind(p_layer_lo, pressure, u_wind, v_wind, N, false);
        WindComponents winds_hi =
            mean_wind(p_layer_hi, pressure, u_wind, v_wind, N, false);

        const float shear_u = winds_hi.u - winds_lo.u;
        const float shear_v = winds_hi.v - winds_lo.v;
        const float mag = vector_magnitude(shear_u, shear_v);

        float storm_u = MISSING;
        float storm_v = MISSING;

        if (leftMover) {
            storm_u = layer_mean_wind.u - ((deviation / mag) * shear_v);
            storm_v = layer_mean_wind.v + ((deviation / mag) * shear_u);
        } else {
            storm_u = layer_mean_wind.u + ((deviation / mag) * shear_v);
            storm_v = layer_mean_wind.v - ((deviation / mag) * shear_u);
        }
        return {storm_u, storm_v};
    }
};

//
////////////  END FUNCTORS   ///////////

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
 * \param   pcl         a sharp::Parcel with its sharp::LPL/attributes defined
 * \param   pres_arr    The pressure coordinate array (Pa)
 * \param   hght_arr    The height coordinate array (meters)
 * \param   buoy_arr    The profile buoyancy array (m/s^2)
 * \param   N           The length of the arrays
 */
void find_lfc_el(Parcel* pcl, const float pres_arr[], const float hght_arr[],
                 const float buoy_arr[], const int N);

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
 * \param   pres_arr    Array of pressure (Pa)
 * \param   hght_arr    Array of height (meters)
 * \param   buoy_arr    Array of buoyancy (m/s^2)
 * \param   N           Length of arrays
 * \param   pcl         A sharp::Parcel corresponding to the buoyancy array.
 */
void cape_cinh(const float pres_arr[], const float hght_arr[],
               const float buoy_arr[], const int N, Parcel* pcl);

}  // end namespace sharp

#endif  // __SHARP_PARCEL_H__
