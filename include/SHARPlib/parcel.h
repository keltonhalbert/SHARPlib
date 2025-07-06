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
#ifndef SHARP_PARCEL_H
#define SHARP_PARCEL_H

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

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
     * \brief Perform the setup step for the parcel lifter.
     *
     * Some parcel lifters require setup in order to handle
     * adiabatic ascent and tracking of vapor, liquid, and
     * ice mixing ratios. The Wobus lifter does not require
     * this, however, so this function does nothing.
     */
    inline void setup([[maybe_unused]] const float lcl_pres,
                      [[maybe_unused]] const float lcl_tmpk) const {}

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
     * \param   tmpk        Parcel temperature (K)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature(
        const float pres, const float tmpk) const {
        return virtual_temperature(tmpk, mixratio(pres, tmpk));
    }
};

struct lifter_cm1 {
   private:
    /**
     * \brief Used to keep track of mixing ratio for conserved/adiabatic lifting
     */
    float rv_total = MISSING;

    /**
     * \brief Water vapor mixing ratio variable updated during parcel lifts
     */
    float rv = MISSING;

    /**
     * \brief Liquid water mixing ratio variable updated during parcel lifts
     */
    float rl = MISSING;

    /**
     * \brief Ice water mixing ratio variable updated during parcel lifts
     */
    float ri = MISSING;

   public:
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
     * \brief perform the necessary setup for parcel ascent.
     *
     * This function sets the total water mixing ratio for
     * adiabatic parcel ascent, and zeroes out the vapor,
     * liquid, and ice mixing ratios from previous parcel
     * ascents.
     */
    inline void setup(const float lcl_pres, const float lcl_tmpk) {
        this->rv_total = mixratio(lcl_pres, lcl_tmpk);
        this->rv = 0.0f;
        this->rl = 0.0f;
        this->ri = 0.0f;
    }

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
        return moist_adiabat_cm1(pres, tmpk, new_pres, this->rv_total, this->rv,
                                 this->rl, this->ri, this->pressure_incr,
                                 this->converge, this->ma_type);
    }

    /**
     * \brief Computes the virtual temperature of the parcel
     * \param   pres        Parcel pressure (Pa)
     * \param   tmpk        Parcel temperature (K)
     *
     * \return  The virtual temperature of the parcel
     */
    [[nodiscard]] inline float parcel_virtual_temperature(
        [[maybe_unused]] const float pres, const float tmpk) const {
        return virtual_temperature(tmpk, this->rv, this->rl, this->ri);
    }
};

/**
 * \author Amelia Urquhart - OU-SoM
 *
 * \brief Lifts a sharp::Parcel according to the lifting formulas from 
 * Peters et al. 2022
 */
struct lifter_peters_et_al {
    /**
     * \brief The type of moist ascent to use.
     */
    ascent_type ma_type = ascent_type::adiab_entr;

    /**
     * \brief The entrainment rate used by the lifter.
     * 
     * Can be either automatically determined or manually set.
     */
    float entr_rate = MISSING;

    /**
     * \brief The pressure increment (Pa) to use for the iterative solver.
     */
    float pressure_incr = 500.0f;

    /**
     * \brief The warm limit of the mixed phase range of temperatures.
     * 
     * The recommended value is 273.15 K, but this is customizable.
     */
    float mixed_phase_warm = 273.15;

    /**
     * \brief The cold limit of the mixed phase range of temperatures.
     * 
     * The recommended value is 253.15 K, but this is customizable
     */
    float mixed_phase_cold = 253.15; 

    /**
     * \brief Internal temperature updated during parcel lifts.
     * 
     * Keeps track of the last temperature lifted to.
     */
    float temperature = MISSING;

    /**
     * \brief Internal pressure updated during parcel lifts.
     * 
     * Keeps track of the last temperature lifted to.
     */
    float pressure = MISSING;

    /**
     * \brief Water vapor mass fraction updated during parcel lifts.
     */
    float qv = MISSING;

    /**
     * \brief Total vapor (water vapor + cloud condensate) mass fraction 
     * variable updated during parcel lifts.
     */
    float qt = MISSING;

    /**
     * \brief Environmental profile used by the lifter.
     */
    // Input variables
    float* prof_pres;
    float* prof_hght;
    float* prof_tmpk;
    float* prof_dwpk;
    float* prof_uwin;
    float* prof_vwin;
    int prof_NZ;

    // Derived variables (NEED TO WRITE CODE THAT DEFINES THESE)
    float* prof_mixr;
    float* prof_vtmp;
    float* prof_theta;
    float* prof_theta_e;
    float* prof_moist_static_energy;

    /**
     * \brief Overloads operator() to call moist_adiabat_peters_et_al()
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

        // Checks whether lifting can be continued from the last time 
        // operator was called. Massively increases efficiency if profile 
        // is sorted in order of increasing height.
        if(new_pres > pressure) {
            can_do_efficient_lift = false;
        }

        float pcl_tmpk;

        if(can_do_efficient_lift) {
            pcl_tmpk = moist_adiabat_peters_et_al(pressure, temperature, 
                new_pres,  this->qv, this->qt, prof_pres, prof_tmpk, prof_mixr, 
                prof_NZ, this->pressure_incr, this->entr_rate, this->ma_type, 
                this->mixed_phase_warm, this->mixed_phase_cold);
        } else {
            pcl_tmpk = moist_adiabat_peters_et_al(pres, tmpk, new_pres, 
                this->qv, this->qt, prof_pres, prof_tmpk, prof_mixr, 
                prof_NZ, this->pressure_incr, this->entr_rate, this->ma_type, 
                this->mixed_phase_warm, this->mixed_phase_cold);
        }

        pressure = new_pres;
        temperature = pcl_tmpk;

        return density_temperature(pcl_tmpk, qv, qt);
    }

    /**
     * This must be run to give the lifter an environmental profile. The
     * profile is used in the formula given by Peters whether entrainment is
     * used or not, so this must be given either way. Otherwise
     * a segmentation fault will occur.
     */
    void set_profile(float* pres, float* hght, float* tmpk, float* dwpk, 
        float* uwin, float* vwin, int NZ) {
        prof_pres = pres;
        prof_hght = hght;
        prof_tmpk = tmpk;
        prof_dwpk = dwpk;
        prof_uwin = uwin;
        prof_vwin = vwin;
        prof_NZ = NZ;

        prof_mixr = new float[NZ];
        prof_vtmp = new float[NZ];
        prof_theta = new float[NZ];
        prof_theta_e = new float[NZ];
        prof_moist_static_energy = new float[NZ];
        
        // VERY IMPORTANT!!!! must compute needed derived state variables
        for(int idx = 0; idx < NZ; idx++) {
            prof_mixr[idx] = sharp::mixratio(pres[idx], dwpk[idx]);
            prof_vtmp[idx] = sharp::virtual_temperature(tmpk[idx], 
                prof_mixr[idx]);
            prof_theta[idx] = sharp::theta(pres[idx], tmpk[idx], 
                sharp::THETA_REF_PRESSURE);
            prof_theta_e[idx] = sharp::thetae(pres[idx], tmpk[idx], dwpk[idx]);
            float specific_humidity = sharp::specific_humidity(prof_mixr[idx]);
            prof_moist_static_energy[idx] = sharp::moist_static_energy(
                hght[idx], tmpk[idx], specific_humidity);
        }
    }

    /**
     * \author Amelia Urquhart - OU-SoM
     *
     * \brief Automatically determines the appropriate entrainment rate.
     *
     * While defaults are given, this function allows for configuration of
     * custom inflow layers and a left-mover option. This may allow for
     * more appropriate entrainment rates for situations such as elevated
     * convection.
     * 
     * It is highly recommended that the user runs this function for each 
     * LPL used by the lifter to ensure that the entrainment rate is as 
     * accurate as possible.
     *
     * \param   prof             (sharp::profile*)
     * \param   lpl              (sharp::LPL)
     * \param   inflow_bottom    (meters)
     * \param   inflow_top       (meters)
     * \param   left_mover       (bool)
     */
    void determine_entrainment_rate(LPL lpl,
                        const float inflow_bottom = 0.0f, 
                        const float inflow_top = 1000.0f,
                        const bool left_mover = false) {
        Parcel pcl;

        define_parcel(prof_pres, prof_tmpk, prof_dwpk, prof_mixr,
                            prof_theta, prof_theta_e, prof_NZ, pcl,
                            lpl);

        lifter_peters_et_al lifter;
        lifter.entr_rate = 0;

        if(ma_type == ascent_type::adiab_entr 
                || ma_type == ascent_type::adiab_nonentr) {
            lifter.ma_type = ascent_type::adiab_nonentr;
        }

        if(ma_type == ascent_type::pseudo_entr 
                || ma_type == ascent_type::pseudo_nonentr) {
            lifter.ma_type = ascent_type::pseudo_nonentr;
        }

        float* pcl_vtmpk_arr = new float[prof_NZ];
        pcl.lift_parcel(lifter, prof_pres, pcl_vtmpk_arr, (std::ptrdiff_t) prof_NZ);
        float* pcl_buoy_arr = new float[prof_NZ];
        sharp::buoyancy(pcl_vtmpk_arr, prof_vtmp, pcl_buoy_arr, (std::ptrdiff_t) prof_NZ);
        
        pcl.cape_cinh(prof_pres, prof_hght, pcl_buoy_arr, prof_NZ);

        entr_rate = compute_epsilon(prof_pres, prof_hght, prof_tmpk, 
            prof_moist_static_energy, prof_uwin, prof_vwin, prof_NZ, 
            &pcl, inflow_bottom, inflow_top, left_mover);
    }


    void define_parcel(const float pressure[], const float temperature[],
                   const float dewpoint[], const float wv_mixratio[],
                   const float theta_arr[], const float thetae[], const int N,
                   Parcel& pcl, LPL source) {
        pcl.source = source;

        if (source == LPL::SFC) {
            pcl.pres = pressure[0];
            pcl.tmpk = temperature[0];
            pcl.dwpk = dewpoint[0];
            return;
        } else if (source == LPL::FCST) {
            // TO-DO: Write the forecast_surface routine
            return;
        } else if (source == LPL::MU) {
            // Search for the most unstable parcel in the bottom
            // 400 hPa of the profile
            static constexpr float mu_depth = 40000.0f;  // 400 hPa in Pa
            PressureLayer mu_layer(pressure[0], pressure[0] - mu_depth);

            // layer_max returns the max, and will set the pressure
            // of the max via a pointer to a float.
            layer_max(mu_layer, pressure, thetae, N, &(pcl.pres));
            pcl.tmpk = interp_pressure(pcl.pres, pressure, temperature, N);
            pcl.dwpk = interp_pressure(pcl.pres, pressure, dewpoint, N);
            return;
        } else if (source == LPL::ML) {
            static constexpr float ml_depth = 10000.0;  // 100 hPa in Pa
            const float sfcpres = pressure[0];
            PressureLayer mix_layer(sfcpres, sfcpres - ml_depth);

            // get the mean attributes of the lowest 100 hPa
            const float mean_mixr = layer_mean(mix_layer, pressure, wv_mixratio, N);
            const float mean_thta = layer_mean(mix_layer, pressure, theta_arr, N);

            // set the parcel attributes
            pcl.pres = sfcpres;
            pcl.tmpk = theta(THETA_REF_PRESSURE, mean_thta, sfcpres);
            pcl.dwpk = temperature_at_mixratio(mean_mixr, sfcpres);
            return;
        } else if (source == LPL::USR) {
            // do nothing - its already been set!
            return;
        } else {
            // TO-DO: probably should raise an error or something
            return;
        }
    }

    // same as sharp::entrainment_rate, has to be copied here outside of convective.h to avoid circular include
    float compute_epsilon(const float pressure[], const float height[],
                        const float temperature[], const float mse_arr[],
                        const float u_wind[], const float v_wind[], const int N,
                        Parcel *pcl, const float inflow_bottom = 0.0f, 
                        const float inflow_top = 1000.0f,
                        const bool left_mover = false) {
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
                                {0.0, 6000.0}, left_mover, false);

        // get the mean inflow layer storm relative wind (default: 0-1km AGL)
        HeightLayer layer = {hsfc + inflow_bottom, hsfc + inflow_top};
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
        constexpr float sigma = 1.1;
        constexpr float alpha = 0.8;
        const float pitchfork = (VKSQ * (alpha * alpha) * (PI * PI) * L) /
                                (4 * PRANDTL * (sigma * sigma) * H);
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
    Parcel(const float pressure, const float temperature, const float dewpoint,
           const LPL lpl);

    /**
     * \author Kelton Halbert - NWS Storm Prediction Center
     *
     * \brief Lifts a sharp::Parcel to compute its virtual temperature
     *
     * Lifts a sharp::Parcel dry adiabatically from its sharp::LPL to its
     * LCL dry adiabatically, and then moist adiabatically from the
     * LCL to the top of the profile. The moist adiabat used is determined
     * bu the type of lifting functor passed to the function (i.e.
     * sharp::lifter_wobus or sharp::lifter_cm1).
     *
     * \param   liftpcl         Parcel lifting function/functor
     * \param   pressure_arr    Array of env pressure (Pa)
     * \param   pcl_vtmpk_arr   The array to fill with virtual temp (K)
     * \param   N               The length of the arrays
     */
    template <typename Lft>
    void lift_parcel(Lft& liftpcl, const float pressure_arr[],
                     float pcl_vtmpk_arr[], const std::ptrdiff_t N) {
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
            pcl_vtmpk_arr[k] = 0.0f;
        }

        // Virtual potential temperature (Theta-V)
        // is conserved for a parcels dry ascent to the LCL
        for (std::ptrdiff_t k = dry_idx.kbot; k < dry_idx.ktop + 1; ++k) {
            const float pcl_vtmp =
                theta(THETA_REF_PRESSURE, thetav_lcl, pressure_arr[k]);
            pcl_vtmpk_arr[k] = pcl_vtmp;
        }

        float pres_bot = pres_lcl;
        float tmpk_bot = tmpk_lcl;
        liftpcl.setup(pres_lcl, tmpk_lcl);

        // fill the array with the moist parcel buoyancy
        for (std::ptrdiff_t k = sat_idx.kbot; k < N; ++k) {
            // compute above-lcl buoyancy here
            const float pcl_pres = pressure_arr[k];
            const float pcl_tmpk = liftpcl(pres_bot, tmpk_bot, pcl_pres);

            if constexpr (!Lft::lift_from_lcl) {
                pres_bot = pcl_pres;
                tmpk_bot = pcl_tmpk;
            }

            const float pcl_vtmpk =
                liftpcl.parcel_virtual_temperature(pcl_pres, pcl_tmpk);
            pcl_vtmpk_arr[k] = pcl_vtmpk;
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
     * \param   mix_layer           sharp::PressureLayer or sharp::HeightLayer
     * \param   pressure            Array of pressure (Pa)
     * \param   height              Array of height (meters)
     * \param   pot_temperature     Array of potential temperature (K)
     * \param   wv_mixratio         Array of water vapor mixing ratio (unitless)
     * \param   N                   Length of arrays
     * \return sharp::Parcel with mixed-layer values
     */
    template <typename Lyr>
    static Parcel mixed_layer_parcel(Lyr& mix_layer, const float pressure[],
                                     [[maybe_unused]] const float height[],
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
     * \param   pcl_virtemp         Writeable array for parcel lifting calcs (K)
     * \param   buoy_arr            Writeable array for buoyancy calcs (m/s^2)
     * \param   N                   Length of arrays
     * \param   search_layer        sharp::PressureLayer or sharp::HeightLay
     * \param   lifter              The parcel moist adiabatic ascent function
     * \return The most unstable sharp::Parcel within the search layer
     */
    template <typename Lyr, typename Lft>
    static Parcel most_unstable_parcel(
        Lyr& search_layer, Lft& lifter, const float pressure[],
        const float height[], const float temperature[], const float virtemp[],
        const float dewpoint[], float pcl_virtemp[], float buoy_arr[],
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

            pcl.lift_parcel(lifter, pressure, pcl_virtemp, N);
            buoyancy(pcl_virtemp, virtemp, buoy_arr, N);
            pcl.cape_cinh(pressure, height, buoy_arr, N);
            if (pcl.cape > max_parcel.cape) max_parcel = pcl;
        }

        return max_parcel;
    }
};

}  // end namespace sharp

#endif  // SHARP_PARCEL_H
