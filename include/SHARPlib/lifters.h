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
#ifndef __SHARP_LIFTERS_H__
#define __SHARP_LIFTERS_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>


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
        Profile* profile;

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
                    new_pres,  this->qv, this->qt, this->profile, 
                    this->pressure_incr, this->entr_rate, this->ma_type, 
                    this->mixed_phase_warm, this->mixed_phase_cold);
            } else {
                pcl_tmpk = moist_adiabat_peters_et_al(pres, tmpk, new_pres, 
                    this->qv, this->qt, this->profile, this->pressure_incr, 
                    this->entr_rate, this->ma_type, this->mixed_phase_warm, 
                    this->mixed_phase_cold);
            }

            pressure = new_pres;
            temperature = pcl_tmpk;

            return density_temperature(pcl_tmpk, qv, qt);
        }

        /**
         * This must be run to give the lifter an environmental profile.
         * 
         * This must be used whether entrainment is turned on or not, otherwise
         * a segmentation fault will occur.
         */
        void set_profile(Profile* prof) {
            profile = prof;
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
        void determine_entrainment_rate(Profile* prof, LPL lpl,
                            const float inflow_bottom = 0.0f, 
                            const float inflow_top = 1000.0f,
                            const bool left_mover = false) {
            Parcel pcl;

            define_parcel(prof->pres, prof->tmpk, prof->dwpk, prof->mixr,
                                prof->theta, prof->theta_e, prof->NZ, pcl,
                                lpl);

            lifter_peters_et_al lifter;
            lifter.profile = prof;
            lifter.entr_rate = 0;

            if(ma_type == ascent_type::adiab_entr 
                    || ma_type == ascent_type::adiab_nonentr) {
                lifter.ma_type = ascent_type::adiab_nonentr;
            }

            if(ma_type == ascent_type::pseudo_entr 
                    || ma_type == ascent_type::pseudo_nonentr) {
                lifter.ma_type = ascent_type::pseudo_nonentr;
            }

            lift_parcel(lifter, prof->pres, prof->vtmp, prof->buoyancy,
                            prof->NZ, &pcl);
            cape_cinh(prof->pres, prof->hght, prof->buoyancy, prof->NZ, &pcl);

            entr_rate = compute_epsilon(prof->pres, prof->hght, prof->tmpk, 
                prof->moist_static_energy, prof->uwin, prof->vwin, prof->NZ, 
                &pcl, inflow_bottom, inflow_top, left_mover);
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
}

//
////////////  END FUNCTORS   ///////////

#endif  // __SHARP_LIFTERS_H__