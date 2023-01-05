
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
#include <SHARPlib/utils.h>
#include <SHARPlib/profile.h>

// Trying to mentally sketch out how I want parcel lifting
// and CAPE/CINH integration to work. I want it to be modular
// such that different methods of computing moist adiabats can
// be supported. I want to separate the parcel lifting from the
// actual numerical integration. I want to be able to store and
// reuse the temperature/pressure traces of parcel lifting. Need 
// to support a "fast CAPE/CINH" for effective inflow calculations.
// 
//
// 1) Create the parcel struct
// 2) Define its starting attributes (MU/ML/SFC/etc)
// 3) Lift the parcel/compute temperature trace
// 4) Set LCL/LFC/EL values
// 5) Integrate CAPE/CINH
//
// So, the previous attempt at this workflow didn't work. In order to separate 
// the lifting of the parcel from the actual integration, it required storing
// the parcel virtual temperature values at corresponding profile pressure/height
// levels. Since the size of the profile is initially unknown, and the height of
// the LCL is unknown, the ultimate size of the parcel temperature, pressure, and
// height traces would have required heap allocations. This is clearly a terrible
// idea for a routine that is meant to be called within loops, whether it be the
// Effective Inflow Layer, or working with gridded data, allocating heap memory
// with every iteration is just incredibly stupid and slow. 
//
// So, after some thinking and research, I think a reasonable solution that
// gives flexibility for both parcel lifting routines and integration mechanisms
// is to use C++ functors and tempalte functions. A functor is essentially just
// an object that overloads the () operator. What's neat is it can be inlined
// by the compiler, and passed as a template argument to, say, the parcel
// lifting routine. 

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
    float operator()(float pres, float tmpc, float new_pres) const { 
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
     * \brief Parcel Convective Available Potential Energy (J/kg)
     */
    float cape;

    /**
     * \brief Parcel Convective Inhibition (J/kg)
     */
    float cinh;

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
void define_parcel(Profile* prof, Parcel* pcl, LPL source);


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
float cinh_below_lcl(Profile* prof, Parcel* pcl, float pres_lcl, float tmpc_lcl);


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Lifts and integrates a parcel to compute various indices
 *
 * Lifts a parcel from its sharp::LPL to the top of the given sharp::Profile.
 * Integrates CAPE and CINH, and computes the parcel LCL, LFC, and EL. The
 * virtual temperature correction is used where possible - see Doswell 1994. 
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
void integrate_parcel(Lifter liftpcl, Profile* prof, Parcel* pcl) {

    // Lift the parcel from the LPL to the LCL 
    float pres_lcl; 
    float tmpc_lcl;
    drylift(pcl->pres, pcl->tmpc, pcl->dwpc, pres_lcl, tmpc_lcl);

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

    float tmpc_pcb = tmpc_lcl;
    float vtmp_enb = interp_pressure(pres_lcl, prof->pres, prof->vtmp, prof->NZ);
	float vtmp_pcb = virtual_temperature(pres_lcl, tmpc_lcl, tmpc_lcl);

    float tmpc_pct = MISSING;
    float vtmp_ent = MISSING;
	float vtmp_pct = MISSING;

    float lyre = 0.0;
    float lyre_last = 0.0;
    float cape = 0.0;
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

        tmpc_pct = liftpcl(pbot, tmpc_pcb, ptop);
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
        }
        else {
            if (ptop > 500.0) cinh += lyre;
        }

        // check for the LFC
        // TO-DO: If there is a layer of positive buoyancy between
        // the LCL and the true LFC, we need to reset CAPE values to
        // not include that layer. I think. 
        if ((lyre >= 0) && (lyre_last <= 0)) {
            // Set the LFC pressure to the 
            float lfc_pres = pbot;
            while (interp_pressure(lfc_pres, prof->pres, prof->vtmp, prof->NZ) 
                    > virtual_temperature(lfc_pres,
                        liftpcl(ptop, tmpc_pct, lfc_pres), 
                        liftpcl(ptop, tmpc_pct, lfc_pres))) {
                lfc_pres -= 5;
            }

            pcl->lfc_pressure = lfc_pres;
        }

        // check for the EL
        if ((lyre <= 0) && (lyre_last >= 0)) {
            float el_pres = pbot;
            while (interp_pressure(el_pres, prof->pres, prof->vtmp, prof->NZ) 
                    < virtual_temperature(el_pres,
                        liftpcl(ptop, tmpc_pct, el_pres), 
                        liftpcl(ptop, tmpc_pct, el_pres))) {
                el_pres -= 5;
            }

            pcl->eql_pressure = el_pres;
            if ((pcl->eql_pressure < pcl->lfc_pressure) && (ptop < 500)) break;
        }

        // set the top of the current layer to the
        // bottom of the next layer
        pbot = ptop;
        hbot = htop;
        vtmp_enb = vtmp_ent;
        tmpc_pcb = tmpc_pct;
        vtmp_pcb = vtmp_pct;
    }

    if (cape == 0) cinh = 0;
    pcl->cinh = cinh;
    pcl->cape = cape;
}


template <typename Lifter>
void integrate_parcel_layer(Lifter liftpcl, Profile* prof, Parcel* pcl, PressureLayer layer) {
    // Lift the parcel from the LPL to the LCL 
    float pres_lcl; 
    float tmpc_lcl;
    drylift(pcl->pres, pcl->tmpc, pcl->dwpc, pres_lcl, tmpc_lcl);

    pcl->lcl_pressure = pres_lcl;

	// if the LCL is above the entirety of our layer, then we need to change things up
	if (pres_lcl < layer.ptop) {
		pcl->cape = 0;
		pcl->cinh = 0;
		return;
	}

    // define the parcel saturated lift layer to be
    // from the LCL to the top of the layer we are integrating 
    PressureLayer sat_layer(pres_lcl, layer.ptop);

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

    float tmpc_pcb = tmpc_lcl;
    float vtmp_enb = interp_pressure(pres_lcl, prof->pres, prof->vtmp, prof->NZ);
	float vtmp_pcb = virtual_temperature(pres_lcl, tmpc_lcl, tmpc_lcl);

    float tmpc_pct = MISSING;
    float vtmp_ent = MISSING;
	float vtmp_pct = MISSING;

	float buoy_bot = 0.0;
	float dz = 0.0;
	float buoy_top = 0.0;
    float lyre = 0.0;
    float lyre_last = 0.0;
    float cape = 0.0;
    // iterate from the LCL to the top of the layer. Excludes
	// layer.ptop so that we can interpoalte the last layer.  
    for (int k = sat_index.kbot; k <= sat_index.ktop; k++) {
#ifndef NO_QC
        if (prof->tmpc[k] == MISSING) {
            continue;
        }
#endif
        ptop = prof->pres[k]; 
        htop = prof->hght[k]; 

        tmpc_pct = liftpcl(pbot, tmpc_pcb, ptop);
        // parcel is saturated, so temperature and dewpoint are same
        vtmp_pct = virtual_temperature(ptop, tmpc_pct, tmpc_pct);
        vtmp_ent = prof->vtmp[k];

        buoy_bot = buoyancy(vtmp_pcb, vtmp_enb);
        buoy_top = buoyancy(vtmp_pct, vtmp_ent);

        dz = htop - hbot;

        lyre_last = lyre;
        lyre = ( ( buoy_bot + buoy_top ) / 2.0 ) * dz;

        if (lyre > 0) {
            cape += lyre;
        }

        // check for the LFC
		// TO-DO: This may not be needed for the buoyancy-layer
		// type stuff... HOWEVER... if there is a positive layer of
		// buoyancy above the LCL but below the true LFC, technically it
		// shouldn't be included in CAPE. So... going to leave this check
		// for the LFC so that when I inevitably have to debug this, its 
		// already there and good to go. 
        if ((lyre >= 0) && (lyre_last <= 0)) {
            // Set the LFC pressure to the 
            float lfc_pres = pbot;
            while (interp_pressure(lfc_pres, prof->pres, prof->vtmp, prof->NZ) 
                    > virtual_temperature(lfc_pres,
                        liftpcl(ptop, tmpc_pct, lfc_pres), 
                        liftpcl(ptop, tmpc_pct, lfc_pres))) {
                lfc_pres -= 5;
            }

            pcl->lfc_pressure = lfc_pres;
        }

        // set the top of the current layer to the
        // bottom of the next layer
        pbot = ptop;
        hbot = htop;
        vtmp_enb = vtmp_ent;
        tmpc_pcb = tmpc_pct;
        vtmp_pcb = vtmp_pct;
    }

	// now integrate the last layer using the interpolated
	// layer top
	ptop = layer.ptop; 
    htop = interp_pressure(ptop, prof->pres, prof->hght, prof->NZ);

	tmpc_pct = liftpcl(pbot, tmpc_pcb, ptop);
	// parcel is saturated, so temperature and dewpoint are same
	vtmp_pct = virtual_temperature(ptop, tmpc_pct, tmpc_pct);
	vtmp_ent = interp_pressure(ptop, prof->pres, prof->vtmp, prof->NZ);

	buoy_bot = buoyancy(vtmp_pcb, vtmp_enb);
	buoy_top = buoyancy(vtmp_pct, vtmp_ent);

	dz = htop - hbot;

	lyre_last = lyre;
	lyre = ( ( buoy_bot + buoy_top ) / 2.0 ) * dz;

	if (lyre > 0) {
		cape += lyre;
	}

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
void parcel_wobf(Profile* prof, Parcel* pcl);


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
