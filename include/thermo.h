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
 * \brief Computes the difference between the wet-bulb potential temperatures for saturated and dry air given the temperature.
 *
 * The Wobus Function (wobf) is defined as the difference between
 * the wet-bulb potential temperature for saturated air (WBPTS)
 * and the wet-bulb potential temperature for dry air (WBPTD) given
 * the same temperature in Celsius.
 *
 * WOBF(T) := WBPTS - WBPTD
 *
 * Although WBPTS and WBPTD are functions of both pressure and
 * temperature, their difference is a function of temperature
 * only. The difference is also proportional to the heat imparted 
 * to a parcel.
 *
 * This function uses a polynomial approximation to the wobus function,
 * fitted to values in Table 78 of PP.319-322 of the Smithsonian Meteorological 
 * Tables by Roland List (6th Revised Edition). Herman Wobus, a mathematician
 * for the Navy Weather Research Facility in Norfolk, VA computed these 
 * coefficients a very long time ago, as he was retired as of the time of 
 * the documentation found on this routine written in 1981.
 *
 * It was shown by Robert Davies-Jones (2007) that the Wobus function has
 * a slight dependence on pressure, which results in errors of up to 1.2
 * degrees Kelvin in the temperature of a lifted parcel. 
 *
 * \param    temperature                     (degC)
 * \return   Sat. Pot. Temperature of Parcel (degC)
 */
float wobf(float temperature);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the vapor pressure over liquid water
 *
 * Computes the vapor pressure (or saturation vapor pressure) in
 * millibars over liquied water given the temperature in Celsius 
 * (or dewpoint temperature in Celsius).
 *
 * This function uses a polynomial fit approximated by Herman Wobus,
 * where the coefficients were chosen to fit the values in Table 94 
 * on PP. 351-353 of the Smithsonian Meteorological Tables by Roland
 * List (6th Edition). 
 *
 * The approximation is valid for -50 C < T < 100 C.
 *
 * \param    temperature    (degC)
 * \return   vapor_pressure (mb) 
 */
float vapor_pressure(float temperature);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the temperature of the LCL.
 *
 * Computes the temperature of a parcels LCL in 
 * Celsius given the parcel temperature and 
 * dewpoint in Celsius. 
 *
 * This is a third-order polynomial approximation written by Herman
 * Wobus, a mathematician for the Navy Weather Research Facility in 
 * Norfolk, VA. He was retired as of 1981, the time when the 
 * documentation on this function was written. The source data for
 * fitting the polynomial is unknown.
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
 * The formula is given by Table 1 on page 7 of Stipanuk (1973).
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
 * This function approximates the water vapor mixing ratio.
 * The formula is given by P. 302 of the Smithsonian
 * Meteorological Tables by Roland List (6th Edition).
 * The function also uses a correction factor (WFW) 
 * computed by Herman Wobus for the departure of the 
 * mixture of air and water vapor from the ideal gas
 * law. The correction forumla fits values in Table 89,
 * P. 340 of the Smithsonian Meteorological Tables, but 
 * only for pressures and temperatures normally 
 * encountered in Earth's atmosphere. 
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
 * \brief Compute the temperature along a moist adiabat (wet-bulb potential temperature) at a given pressure
 *
 * Compute the temperature at which the moist adiabat intersects a line 
 * of constant pressure on a Skew-T log-P diagram. The wet-bulb potential 
 * temperature, given by theta_sat, defines a moist adiabat in Celsius, and
 * the temperature at the given pressure level in millibars is returned.
 *
 * This function relies on the Wobus Function, and it was shown by Robert 
 * Davies-Jones (2007) that the Wobus function has a slight dependence on 
 * pressure, which results in errors of up to 1.2 degrees Kelvin in the 
 * temperature of a lifted parcel. 
 *
 * \param    pressure              (mb)
 * \param    theta_sat             (degC)
 * \return   lifted_temperature    (degC) 
 */
float saturated_lift(float pressure, float theta_sat);


/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 * \brief Compute the temperature of a parcel lifted moist adiabatically to a new level. 
 *
 * With a given parcel defined by a pressure and temperature (in millibars and Celsius), lift it moist adiabatically to a new pressure level (in millibars) and return the temperture of the parcel at that level. 
 *
 * This computation relies on the Wobus Function to compute the moist adiabats,
 * and it was shown by Robert Davies-Jones (2007) that the Wobus function has a 
 * slight dependence on pressure. This can result in errors of up to 1.2 
 * degrees Kelvin int he temperature of a lifted parcel. 
 *
 * \param    pressure              (mb)
 * \param    temperature           (degC)
 * \return   lifted_pressure       (degC) 
 */
float wetlift(float pressure, float temperature, float lifted_pressure);

}
