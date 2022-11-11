
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

#include "profile.h"
#include "thermo.h"
#include "utils.h"

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
 * \brief Lifts a parcel to compute CAPE, CINH, and various parcel levels.
 *
 * Lifts a parcel from its sharp::LPL to the top of the given sharp::Profile.
 * Integrates CAPE and CINH, and computes the parcel LCL, LFC, EL, and MPL. The
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
void lift_parcel(Lifter liftpcl, Profile* prof, Parcel* pcl) {
    // virtual temperature at LPL
    float vtmp_pcl = virtual_temperature(pcl->pres, pcl->tmpc, pcl->dwpc);

    // Lift the parcel from the LPL to the LCL 
    float pres_at_lcl; 
    float tmpc_at_lcl;
    drylift(pcl->pres, pcl->tmpc, pcl->dwpc, pres_at_lcl, tmpc_at_lcl);

    // define the parcel saturated lift layer to be
    // from the LCL to the top of the profile available
    PressureLayer sat_layer(pres_at_lcl, prof->pres[prof->NZ-1]);

    // excludes the indices that would correspond to the exact top and
    // bottom of this layer - default is that bottom and top is interpolated
    LayerIndex sat_index = get_layer_index(sat_layer, prof->pres, prof->NZ);

    // iterate from the LCL to the top of the profile
    float pbot = pres_at_lcl;
    float tbot = virtual_temperature(pres_at_lcl, tmpc_at_lcl, tmpc_at_lcl);
    float ptop = MISSING;
    float ttop = MISSING;
    for (int k = sat_index.kbot; k <= sat_index.ktop; k++) {
#ifndef NO_QC
        if ((prof->pres[k] == MISSING) || (prof->tmpc[k] == MISSING)) {
            continue;
        }
#endif
        ptop = prof->pres[k]; 
        ttop = liftpcl(pbot, tbot, ptop);
        vtmp_pcl = virtual_temperature(ptop, ttop, ttop);

        // set the top of the current layer to the
        // bottom of the next layer
        pbot = ptop;
        tbot = ttop;
    }

    // lift final level
    ttop = liftpcl(pbot, tbot, sat_layer.ptop);
    vtmp_pcl = virtual_temperature(ptop, ttop, ttop); 
}


} // end namespace sharp


namespace sharp::exper {


} // end namespace sharp::exper


#endif
