
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
#ifndef __SHARP_PARCEL
#define __SHARP_PARCEL


#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/profile.h>
#include <iostream>

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
     * \param tmpc      Initial parcel temperature (degC)
     * \param new_pres  Final level of parcel after lift (hPa)
     */
    float operator()(float pres, float tmpc, float new_pres) const noexcept { 
        return wetlift(pres, tmpc, new_pres);
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
    SFC  = 1,

    /**
     * \brief Forecast Surface Parcel
     */
    FCST = 2, 

    /**
     * \brief Most Unstable Parcel
     */
    MU   = 3, // most unstable

    /**
     * \brief 100mb Mixed Layer Parcel
     */
    ML   = 4, // 100mb mixed layer

    /**
     * \brief User-defined Parcel
     */
    USR  = 5, // user-defined

    /**
     * \brief Mean Effective Inflow Layer Parcel
     */
    EIL  = 6, // Mean effective inflow layer
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
    float pres;

    /**
     * \brief Parcel starting temperature (degC)
     */
    float tmpc;

    /**
     * \brief Parcel starting dewpoint (degC)
     */
    float dwpc;

    /**
     * \brief Pressure at the Lifted Condensation Level (hPa)
     */
    float lcl_pressure;

    /**
     * \brief Pressure at the Level of Free Convection (hPa)
     */
    float lfc_pressure;

    /**
     * \brief Pressure at the parcel Equilibrium Level (hPa)
     */
    float eql_pressure;

    /**
     * \brief Pressure at the Maximum Parcel Level (hPa)
     */
    float mpl_pressure;

    /**
     * \brief Parcel Convective Available Potential Energy (J/kg) between the LFC and EL
     */
    float cape;

    /**
     * \brief Parcel Convective Inhibition (J/kg) between the LFC and EL
     */
    float cinh;

    /**
     * \brief Parcel Convective Available Potnential Energy (J/kg) 
     * for the total profile depth (all positive areas above LCL)
     */
    float total_cape;

    /**
     * \brief Parcel Convective Inhibition (J/kg) 
     * for the total profile depth (all negative areas below EL)
     */
    float total_cinh;

    /**
     * \brief The type of parcel this is
     */
    LPL source; 
    
    /**
     * \brief Parcel empty constructor that sets all values to sharp::MISSING
     */
    Parcel();
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
 * \brief Computes the CINH from the sharp::LPL to the LCL
 *
 * Accumulates the CINH from the sharp::LPL to the parcel's LCL in 10 hPa
 * increments. The virtual temperature correction is used to properly 
 * account for additional buoyancy from water vapor (see Doswell 1994). 
 * This routine assumes you have either defined the parcel using 
 * sharp::define_parcel or have defined the parcel attributes manually.
 *
 * \param prof      sharp::Profile
 * \param pcl       sharp::Parcel
 * \param pres_lcl  Pressure at LCL (hPa)
 * \param tmpc_lcl  Temperature at LCL (degC)
 * \return Convective Inhibition (CINH; J/kg) below the LCL
 */
float cinh_below_lcl(Profile* prof, Parcel* pcl, float pres_lcl, 
                     float tmpc_lcl) noexcept;
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
 */
template <typename Lft>
void lift_parcel(Lft liftpcl, Profile* prof, Parcel* pcl) noexcept {
    // Lift the parcel from the LPL to the LCL 
    float pres_lcl; 
    float tmpc_lcl;

    drylift(pcl->pres, pcl->tmpc, pcl->dwpc, pres_lcl, tmpc_lcl);
    pcl->lcl_pressure = pres_lcl;

    float thetav_lcl = theta(
            pres_lcl, 
            virtual_temperature(pcl->pres, tmpc_lcl, tmpc_lcl),
            1000.0
        ); 

    // define the dry and saturated lift layers
    PressureLayer dry_lyr = {pcl->pres, pcl->lcl_pressure};
    PressureLayer sat_lyr = {pcl->lcl_pressure, prof->pres[prof->NZ-1]};

    // the LayerIndex excludes the top and bottom for interpolation reasons
    LayerIndex dry_idx = get_layer_index(dry_lyr, prof->pres, prof->NZ);
    LayerIndex sat_idx = get_layer_index(sat_lyr, prof->pres, prof->NZ);

    // Fill the array with dry parcel buoyancy.
    // virtual potential temperature (Theta-V)
    // is conserved for a parcels dry ascent to the LCL
    for (int k = dry_idx.kbot; k <= dry_idx.ktop; ++k) {
        float pcl_pres = prof->pres[k];
        float env_vtmp = prof->vtmp[k];
        float pcl_vtmp = theta(1000.0, thetav_lcl, pcl_pres);
        prof->buoyancy[k] = buoyancy(pcl_vtmp, env_vtmp);
    }

    // fill the array with the moist parcel buoyancy
    for (int k = sat_idx.kbot; k <= sat_idx.ktop; ++k) {
        // compute above-lcl buoyancy here
        float pcl_pres = prof->pres[k];
        float pcl_tmpc = liftpcl(pres_lcl, tmpc_lcl, pcl_pres);
        // parcel is saturated, so temperature and dewpoint are same
        float pcl_vtmp = virtual_temperature(pcl_pres, pcl_tmpc, pcl_tmpc);
        float env_vtmp = prof->vtmp[k];

        prof->buoyancy[k] = buoyancy(pcl_vtmp, env_vtmp);
    }
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Lifts and integrates a parcel to compute CAPE, CINH, and the LFC/EL heights
 *
 * Lifts a parcel from its sharp::LPL to the top of the given sharp::Profile.
 * Integrates CAPE and CINH, and computes the parcel LCL, LFC, and EL. The
 * virtual temperature correction is used where possible - see Doswell 1994. 
 *
 * In edge cases where profiles are close to neutrally buoyant, there can be
 * multiple layers of buoyancy and therefore multiple LFC/EL candidates. This
 * algorithm selects the layer with the largest CAPE, and sets the LFC/EL as 
 * the levels bounding this layer. Total profile buoyancy is also accounted
 * for to support legacy behavior. 
 *
 * The method of computing the moist adiabats for parcel ascent can be passed
 * to this routine via a functor, meaning that other methods of computing moist
 * adiabats may be used without negatively affecting performance. 
 *
 * \param liftpcl   A functor that that provides a means to lift the parcel
 * \param prof      A sharp::Profile of atmospheric data
 * \param pcl       A sharp::Parcel with its sharp::LPL/attributes defined.
 */
template <typename Lifter>
void integrate_parcel(Lifter liftpcl, Profile* prof, Parcel* pcl) noexcept {

    // Lift the parcel from the LPL to the LCL 
    float pres_lcl; 
    float tmpc_lcl;
    drylift(pcl->pres, pcl->tmpc, pcl->dwpc, pres_lcl, tmpc_lcl);

	if (pres_lcl > prof->pres[0]) pres_lcl = prof->pres[0];
    pcl->lcl_pressure = pres_lcl;

    // get the CINH below the LCL in 10 hPa increments
    float cinh = cinh_below_lcl(prof, pcl, pres_lcl, tmpc_lcl);

    // define the parcel saturated lift layer to be
    // from the LCL to the top of the profile available
    PressureLayer sat_layer(pres_lcl, prof->pres[prof->NZ-1]);

    // excludes the indices that would correspond to the exact top and
    // bottom of this layer - default is that bottom and top is interpolated
    LayerIndex sat_index = get_layer_index(sat_layer, prof->pres, prof->NZ);

    float pbot = pres_lcl;
    float hbot = interp_pressure(pres_lcl, prof->pres, prof->hght, prof->NZ);
    float ptop = MISSING;
    float htop = MISSING;

    // _enb - environment bottom layer
    // _ent - environment top layer
    // _pcb - parcel bottom layer
    // _pct - parcel top layer

    //float tmpc_pcb = tmpc_lcl;
    float vtmp_enb = interp_pressure(pres_lcl, prof->pres, prof->vtmp, prof->NZ);
	float vtmp_pcb = virtual_temperature(pres_lcl, tmpc_lcl, tmpc_lcl);

    float tmpc_pct = MISSING;
    float vtmp_ent = MISSING;
	float vtmp_pct = MISSING;

    float lyre = 0.0;
    float lyre_last = 0.0;
    float cape = 0.0;
    float cape_old = 0.0;
	float cinh_old = 0.0;
    float lfc_pres = MISSING;
    float lfc_old = MISSING;
    float el_pres = MISSING;
    float el_old = MISSING;
    pcl->total_cape = 0.0;
    pcl->total_cinh = 0.0;

    // iterate from the LCL to the top of the profile
    // -- the +1 takes us to the last level since it is
    // excluded by get_layer_index
    for (int k = sat_index.kbot; k <= sat_index.ktop+1; k++) {
#ifndef NO_QC
        if (prof->tmpc[k] == MISSING) {
            continue;
        }
#endif
        ptop = prof->pres[k]; 
        htop = prof->hght[k]; 

        tmpc_pct = liftpcl(pres_lcl, tmpc_lcl, ptop);
        // parcel is saturated, so temperature and dewpoint are same
        vtmp_pct = virtual_temperature(ptop, tmpc_pct, tmpc_pct);
        vtmp_ent = prof->vtmp[k];

        float buoy_bot = buoyancy(vtmp_pcb, vtmp_enb);
        float buoy_top = buoyancy(vtmp_pct, vtmp_ent);

        float dz = htop - hbot;

        lyre_last = lyre;
        lyre = ( ( buoy_bot + buoy_top ) / 2.0 ) * dz;

        if (lyre > 0) {
            cape += lyre;
            pcl->total_cape += lyre;
        }
        else {
			// we only want to accumulate CINH
			// between the LCL and LFC. 
			if (pcl->lfc_pressure == MISSING) {
				cinh += lyre;
			}
			// this keeps track of total CINH in case of
			// multiple LFC heights. 
            pcl->total_cinh += lyre;
        }

        // check for the LFC
        if ((lyre > 0) && (lyre_last <=  0)) {
            // Set the LFC pressure to the bottom pressure layer
            lfc_pres = pbot;
			bool found = false;
            // if we can't find the LFC within 50 hPa of this 
			// level, stop searching because it's either the 
			// first level or further up the profile.
			for (lfc_pres = pbot; lfc_pres > pbot-50.0; lfc_pres-=5.0) {
				float env_vtmp = interp_pressure(lfc_pres, prof->pres, prof->vtmp, prof->NZ);
				float pcl_tmpc = liftpcl(ptop, tmpc_pct, lfc_pres);
				float pcl_vtmp = virtual_temperature(lfc_pres, pcl_tmpc, pcl_tmpc);
				if (pcl_vtmp > env_vtmp) {
					found = true;
					break;
				}
			}
			// reset to pbot if not found, since
			// this was the level that triggered 
			// the LFC condition to begin with.  
			if (!found) {
				lfc_pres = pbot;
				found = true;
			}

			//sanity check against the LCL
            if (lfc_pres > pres_lcl) lfc_pres = pres_lcl;

            // We've already found an LFC candidate - need
            // to keep track of a few things to determine which
            // LFC is the "right one".
            if (pcl->lfc_pressure != MISSING) {
                // store the old CAPE value here - we've already
                // added a new positive layer, so pull it back out
                // so we don't accidentally double count later
                cape_old = cape - lyre;
				cinh_old = cinh;
                lfc_old = pcl->lfc_pressure;
				el_old = pcl->eql_pressure;

                // reset the CAPE and integrate from the new LFC
                cape = lyre;
				cinh = pcl->total_cinh;
            }
            if (found) pcl->lfc_pressure = lfc_pres;
        } // end LFC check

        // check for the EL
        if (((lyre < 0) && (lyre_last >= 0)) || 
            ((ptop == prof->pres[prof->NZ-1]) && (lyre > 0))) {
            el_pres = pbot;
			bool found = false;
            // if we can't find the EL within 50 hPa 
			// of this level, stop searching 
			for (el_pres = pbot; el_pres > ptop-50.0; el_pres-=5.0) {
                if (el_pres < prof->pres[prof->NZ-1]) break;
				float env_vtmp = interp_pressure(el_pres, prof->pres, prof->vtmp, prof->NZ);
				float pcl_tmpc = liftpcl(ptop, tmpc_pct, el_pres);
				float pcl_vtmp = virtual_temperature(el_pres, pcl_tmpc, pcl_tmpc);
				if (pcl_vtmp < env_vtmp) {
					found = true;
					break;
				}
			}

			if ((!found) && (pbot < pcl->lfc_pressure)) {
				el_pres = pbot;
				found = true;
			}

			if (found) pcl->eql_pressure = el_pres;

            // check which CAPE layer should stay,
			// and keep the layer with the maximum
			// buoyancy in the profile. 
            bool swapped = false;
            if (cape_old > cape) {
                pcl->lfc_pressure = lfc_old;
                pcl->eql_pressure = el_old;
				cinh = cinh_old; 
                cape = cape_old;
                swapped = true;
            }

            if ((pcl->eql_pressure != MISSING) && (!swapped) && 
                (pcl->eql_pressure < pcl->lfc_pressure)) {
                // take off the last layer if positive and
                // then integrate up to the EL
                if (lyre > 0) {
                    cape -= lyre;
                }

                float el_htop = interp_pressure(pcl->eql_pressure, prof->pres, prof->hght, prof->NZ);
                float el_tmpc_pct = liftpcl(pres_lcl, tmpc_lcl, pcl->eql_pressure);
                // parcel is saturated, so temperature and dewpoint are same
                float el_vtmp_pct = virtual_temperature(pcl->eql_pressure, el_tmpc_pct, el_tmpc_pct);
                float el_vtmp_ent = interp_pressure(pcl->eql_pressure, prof->pres, prof->vtmp, prof->NZ);

                float el_buoy_bot = buoyancy(vtmp_pcb, vtmp_enb);
                float el_buoy_top = buoyancy(el_vtmp_pct, el_vtmp_ent);

                float el_dz = el_htop - hbot;

                float el_lyre = ( ( el_buoy_bot + el_buoy_top ) / 2.0 ) * el_dz;
                if ( el_lyre > 0 ) {
                    cape += el_lyre;
                }

                // no need to keep lifting/integrating past 200 hPa
                // if these conditions have been met. 
                if (ptop < 200) {
                    break;
                }
            }
        } // end EL check

        // set the top of the current layer to the
        // bottom of the next layer
        pbot = ptop;
        hbot = htop;
        vtmp_enb = vtmp_ent;
        //tmpc_pcb = tmpc_pct;
        vtmp_pcb = vtmp_pct;
    } // end profile iteration

	// truncate CAPE values
	// below 1 J/kg to zero
    if (cape < 1) {
        cape = 0;
        cinh = 0;
    }
    if (pcl->total_cape < 1) {
        pcl->total_cape = 0.0;
        pcl->total_cinh = 0.0;
    }
    pcl->cinh = cinh;
    pcl->cape = cape;
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Lifts a parcel using the Wobus method of computing moist adiabats 
 *
 *
 * Lifts a parcel from its sharp::LPL to the Equilibrium Level (EL) using the
 * Wobus method for computing moist adiabats. See sharp::wobf and sharp::wetlift
 * for implementation details. Attributes such as CAPE, CINH, and specific parcel
 * levels are stored within the sharp::Parcel. 
 *
 * \param prof      A sharp::Profile of atmospheric data
 * \param pcl       A sharp::Parcel with its sharp::LPL/attributes defined.
 */
void parcel_wobf(Profile* prof, Parcel* pcl) noexcept;


// TO-DO ECAPE: 
// Step 5 requires a vertical array of integrated MSE, that
// then gets integrated vertically in Step 6 to get NCAPE. Allocating
// memory within entrainment_cape would be unwise for gridded data,
// since calling new/malloc every grid cell is a bad idea for performance. 
// A reasonable solution could be including some vertical dummy/temproary 
// arrays within the Parcel or Profile objects, or include MSE as a vertical
// field that's always available. 

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Computes the NCAPE term used for ECAPE 
 *
 *
 * Computes the buoyancy dilution rate (NCAPE) used for evaluating
 * Entrainment CAPE (ECAPE).
 *
 * \param prof      A sharp::Profile of atmospheric data
 * \param pcl       A sharp::Parcel with its sharp::LPL/attributes defined.
 * \param mse_star  Vertical array of Saturation Moist Static Energy
 * \param mse_bar   Vertical array of integrated Moist Static Energy 
 */
float buoyancy_dilution(Profile* prof, Parcel *pcl, 
                        const float *mse_star, 
                        const float *mse_bar) noexcept;


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
