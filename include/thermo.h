/**
 * \file
 * \brief Routines used to compute thermodynamic attributes of vertical sounding profiles
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
#pragma once

namespace sharp {

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Estimate the saturated potential temperature
 *
 * Compute the Wobus function for moist parcel 
 * ascent given a temperature in degrees Celsius. 
 * Returns the Sat. Pot. Temperature of a parcel 
 * in Celsius.
 *
 * \param    temperature                     (degC)
 * \return   Sat. Pot. Temperature of Parcel (degC)
 */
float wobf(float temperature);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the vapor pressure
 *
 * Computes the vapor pressure of dry air at the
 * given temperature in degrees Celsius.
 *
 * \param    temperature    (degC)
 * \return   vapor_pressure (mb) 
 */
float vapor_pressure(float temperature);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the temperature of the LCL
 *
 * Computes the temperature of a parcels LCL in 
 * Celsius given the parcel temperature and 
 * dewpoint in Celsius. 
 *
 * \param    temperature     (degC)
 * \param    dewpoint        (degC)
 * \return   lcl_temperature (degC) 
 */
float lcl_temperature(float temperature, float dewpoint);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the temperature at a given mixing ratio and pressure level. 
 *
 * Computes the temperature in Celsius of air at the
 * given water vapor mixing ratio in g/kg and the 
 * air pressure in mb.
 *
 * \param    mixratio    (g/kg)
 * \param    pressure    (mb)
 * \return   temperature (degC) 
 */
float temperature_at_mixratio(float mixratio, float pressure);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Comute the pressure level given potential temperature and temperature. 
 *
 * Returns the pressure level in millibars of a parcel
 * given the potential temperature in Celsius and the
 * temperature of the parcel in Celsius.
 *
 * \param    potential_temperature (degC)
 * \param    temperature           (degC)
 * \return   pressure              (mb) 
 */
float theta_level(float potential_temperature, float temperature);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the potential temperature. 
 *
 * Returns the potential temperature in Celsius of
 * a parcel given its pressure in millibars and
 * temperature in Celsius. The final argument is
 * the reference level, which is usually 1000.0 mb.
 *
 * \param    pressure              (mb)
 * \param    temperature           (degC)
 * \param    ref_pressure          (mb)
 * \return   potential_temperature (degC) 
 */
float theta(float pressure, float temperature, float ref_pressure);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the water vapor mixing ratio. 
 *
 * Returns the water vapor mixing ratio in g/kg given 
 * the environmental pressure in millibars and a 
 * temperature (dry-bulb or dewpoint) in Celsius.  
 *
 * \param    pressure              (mb)
 * \param    temperature           (degC)
 * \return   mixratio              (g/kg) 
 */
float mixratio(float pressure, float temperature); 

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the virtual temperature. 
 *
 * Returns the virtual temperature in Celsius given the ambient
 * pressure in millibars, dry-bulb temperature in Celsius, and 
 * the dewpoint temperature in Celsius.
 *
 * \param    pressure              (mb)
 * \param    temperature           (degC)
 * \param    dewpoint              (degC)
 * \return   virtual_temperature   (degC) 
 */
float virtual_temperature(float pressure, float temperature, float dewpoint);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the new temperature of a saturated parcel lifted to a new level. 
 *
 * Returns the lifted temperature in Celsius of a saturated parcel
 * that is lifted to a new pressure level in millibars. Caution should
 * be exercised when using this routine with data at higher resolutions, such
 * as with full resolution radiosondes, as the default convergence criteria
 * can cause drift in the calculations, resulting in erroneous moist adiabats.
 *
 * \param    pressure              (mb)
 * \param    theta_sat             (degC)
 * \return   lifted_temperature    (degC) 
 */
float saturated_lift(float pressure, float theta_sat);

}
