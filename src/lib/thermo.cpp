// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
#pragma once
#include <cmath>
#include "../include/constants.h"

namespace sharp {

/*
 * WOBF
 * Authored by John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * Compute the Wobus function for moist parcel 
 * ascent given a temperature in degrees Celsius. 
 * Returns the Sat. Pot. Temperature of a parcel 
 * in Celsius.
 *
 * @param    float temperature                     (degC)
 * @return   float Sat. Pot. Temperature of Parcel (degC)
 */
inline float wobf(float temperature) {
    float x;
    double pol;

    if (temperature == MISSING) return MISSING;

    x = temperature - 20.0;
    if (x <= 0.0) {
       pol = 1.0 + x * (-8.841660499999999e-03 + x * ( 1.4714143e-04 \
                 + x * (-9.671989000000001e-07 + x * (-3.2607217e-08 \
                 + x * (-3.8598073e-10)))));
       pol = pol * pol;
       return 15.13 / (pol * pol);
    }
    else {
       pol = x * (4.9618922e-07 + x * (-6.1059365e-09 + \
             x * (3.9401551e-11 + x * (-1.2588129e-13 + \
             x * (1.6688280e-16)))));
       pol = 1.0 + x * (3.6182989e-03 + x * (-1.3603273e-05 + pol));
       pol = pol * pol;
       return 29.93 / (pol * pol) + 0.96 * x - 14.8;
    }
}

/*
 * VAPPRES
 * Authored by John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * Computes the vapor pressure of dry air at the
 * given temperature in degrees Celsius.
 *
 * @param    float temperature    (degC)
 * @return   float vapor_pressure (mb) 
 */
inline float vappres(float temperature) {
	double pol;

    if (temperature == MISSING) return MISSING;

	pol = temperature * (1.1112018e-17 + temperature * (-3.0994571e-20));
	pol = temperature * (2.1874425e-13 + temperature * (-1.789232e-15 + pol));
	pol = temperature * (4.3884180e-09 + temperature * (-2.988388e-11 + pol));
	pol = temperature * (7.8736169e-05 + temperature * (-6.111796e-07 + pol));
	pol = .99999683e-00 + temperature * (-9.082695e-03 + pol);
	pol = (pol * pol);
	pol = (pol * pol);
	return 6.1078 / (pol * pol);
}

/*
 * LCLTEMP
 * Authored by John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * Computes the temperature of a parcels LCL in 
 * Celsius given the parcel temperature and 
 * dewpoint in Celsius. 
 *
 * @param    float temperature     (degC)
 * @param    float dewpoint        (degC)
 * @return   float lcl_temperature (degC) 
 */
inline float lcltemp(float temperature, float dewpoint) {
	float s, t, dlt;

    if ((temperature == MISSING) || (dewpoint == MISSING)) {
       return MISSING;
    } 

	s = temperature - dewpoint;
	dlt = s * (1.2185 + .001278 * temperature + s * 
          (-.00219 + 1.173E-05 * s - .0000052 * temperature));
	return temperature - dlt;
}

/*
 * TEMPERATURE_AT_MIXRATIO 
 * Authored by John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * Computes the temperature in Celsius of air at the
 * given water vapor mixing ratio in g/kg and the 
 * air pressure in mb.
 *
 * @param    float mixratio    (g/kg)
 * @param    float pressute    (mb)
 * @return   float temperature (degC) 
 */
inline float temperature_at_mixratio(float mixratio, float pressure) {
    if ((mixratio == MISSING) || (pressure == MISSING)) {
       return MISSING;
    } 

	double c1 = 0.0498646455;
	double c2 = 2.4082965;
	double c3 = 7.07475;
	double c4 = 38.9114;
	double c5 = 0.0915;
	double c6 = 1.2035;

	double x    = std::log10( mixratio * pressure / (622.0 + mixratio));
	double tmrk = std::pow(10.0, c1 * x + c2) - c3 + c4 * 
                  std::pow(std::pow(10.0, c5 * x) - c6, 2.0);

	return (float)(tmrk - ZEROCNK);
}

/*
 * THALVL 
 * Authored by John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * Returns the pressure level in millibars of a parcel
 * given the potential temperature in Celsius and the
 * temperature of the parcel in Celsius.
 *
 * @param    float potential_temperature (degC)
 * @param    float temperature           (degC)
 * @return   float pressure              (mb) 
 */
inline float thalvl(float potential_temperature, float temperature) {
    if ((potential_temperature == MISSING)
               || (temperature == MISSING)) {
       return MISSING;
    } 

	potential_temperature += ZEROCNK;
	temperature           += ZEROCNK;
	return 1000.0 / std::pow((potential_temperature/temperature), (1.0/ROCP));
}


/*
 * THETA 
 * Authored by John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * Returns the potential temperature in Celsius of
 * a parcel given its pressure in millibars and
 * temperature in Celsius. The final argument is
 * the reference level, which is usually 1000.0 mb.
 *
 * @param    float pressure              (mb)
 * @param    float temperature           (degC)
 * @param    float ref_pressure          (mb)
 * @return   float potential_temperature (degC) 
 */
inline float theta(float pressure, float temperature, float ref_pressure) {
    if ((temperature == MISSING) || (pressure == MISSING)
                             || (ref_pressure == MISSING)) {
       return MISSING;
    } 

	temperature += ZEROCNK;
	return (temperature * std::pow(ref_pressure / pressure, ROCP)) - ZEROCNK;
}

}
