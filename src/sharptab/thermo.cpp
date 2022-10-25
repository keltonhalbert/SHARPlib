// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
#include <cmath>
#include "constants.h"

namespace sharp {


float wobf(float temperature) {
    float x;
    double pol;

#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif

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


float vapor_pressure(float temperature) {
	double pol;

#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif

	pol = temperature * (1.1112018e-17 + temperature * (-3.0994571e-20));
	pol = temperature * (2.1874425e-13 + temperature * (-1.789232e-15 + pol));
	pol = temperature * (4.3884180e-09 + temperature * (-2.988388e-11 + pol));
	pol = temperature * (7.8736169e-05 + temperature * (-6.111796e-07 + pol));
	pol = .99999683e-00 + temperature * (-9.082695e-03 + pol);
	pol = (pol * pol);
	pol = (pol * pol);
	return 6.1078 / (pol * pol);
}


float lcl_temperature(float temperature, float dewpoint) {
	float s, t, dlt;

#ifndef NO_QC
    if ((temperature == MISSING) || (dewpoint == MISSING)) {
       return MISSING;
    } 
#endif

	s = temperature - dewpoint;
	dlt = s * (1.2185 + .001278 * temperature + s * 
          (-.00219 + 1.173E-05 * s - .0000052 * temperature));
	return temperature - dlt;
}


float temperature_at_mixratio(float wv_mixratio, float pressure) {
#ifndef NO_QC
    if ((wv_mixratio == MISSING) || (pressure == MISSING)) {
       return MISSING;
    } 
#endif

	constexpr double c1 = 0.0498646455;
	constexpr double c2 = 2.4082965;
	constexpr double c3 = 7.07475;
	constexpr double c4 = 38.9114;
	constexpr double c5 = 0.0915;
	constexpr double c6 = 1.2035;

	double x    = std::log10( wv_mixratio * pressure / (622.0 + wv_mixratio));
	double tmrk = std::pow(10.0, c1 * x + c2) - c3 + c4 * 
                  std::pow(std::pow(10.0, c5 * x) - c6, 2.0);

	return (float)(tmrk - ZEROCNK);
}


float theta_level(float potential_temperature, float temperature) {

#ifndef NO_QC
    if ((potential_temperature == MISSING)
               || (temperature == MISSING)) {
       return MISSING;
    } 
#endif

	potential_temperature += ZEROCNK;
	temperature           += ZEROCNK;
	return 1000.0 / std::pow((potential_temperature/temperature), (1.0/ROCP));
}


float theta(float pressure, float temperature, float ref_pressure) {

#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)
                             || (ref_pressure == MISSING)) {
       return MISSING;
    } 
#endif

	temperature += ZEROCNK;
	return (temperature * std::pow(ref_pressure / pressure, ROCP)) - ZEROCNK;
}


float mixratio(float pressure, float temperature) {

#ifndef NO_QC
    if (( temperature == MISSING ) || (pressure == MISSING)) {
        return MISSING;
    }
#endif

    // correction factor for the departure of the mixture
    // of air and water vapor from the ideal gas law
    float x = 0.02 * (temperature - 12.5 +7500.0 / pressure);
    float wfw = 1.0 + 0.0000045 * pressure + 0.0014 * x * x;
    float fwesw = wfw * vapor_pressure(temperature);
    return 621.97 * (fwesw / (pressure - fwesw));
}


float virtual_temperature(float pressure, float temperature, float dewpoint) {

    // Just an example of how we can remove missing data checks
    // for grid operations by compiling without the QC_DATA flag. 
#ifndef NO_QC
    if (dewpoint == MISSING) {
        return temperature;
    }
    else if ((pressure == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif

    constexpr double eps = 0.62197;
    double temperature_k = temperature + ZEROCNK;
    double w = 0.001 * (double)mixratio(pressure, dewpoint);
    return (float)(temperature_k * (1.0 + w / eps) / (1.0 + w) - ZEROCNK);
}


float saturated_lift(float pressure, float theta_sat) {

#ifndef NO_QC
    if ((pressure == MISSING) || (theta_sat == MISSING)) {
        return MISSING;
    }
#endif
    
    // this function is about lifting a parcel, so if we are lifting it below
    // 1000.0 hPa then something is probably wrong
	if ((std::fabs(pressure - 1000.0) - 0.001) <= 0.0) return theta_sat;

	float pwrp, t1, t2, woto, wotm, wot2, woe2, e1, e2, rate;
	float eor = 999;

    // TO-DO: Pass the 0.1 as an argument or 
    // define it as a constant, since this controls
    // the rate of convergence and therefore the quality
    // of the calculation
	while (std::fabs(eor) - 0.1 > 0.0) {
        
        // TO-DO: This section can be optimized because it is
        // only called once. There is no need to check this
        // iff condition multiple times to achieve the result.
        if (eor == 999) {                    /* First Pass */
	        pwrp = std::pow(pressure / 1000.0, ROCP);
            // get the temperature 
	        t1 = (theta_sat + ZEROCNK) * pwrp - ZEROCNK;
	        woto = wobf(t1);
	        wotm = wobf(theta_sat);
	        e1 = woto - wotm;
	        rate = 1;
        }
        else {  /* Successive Passes */
            rate = (t2 - t1) / (e2 - e1);
	        t1 = t2;
	        e1 = e2;
        }
        t2 = t1 - e1 * rate;
	    e2 = (t2 + ZEROCNK) / pwrp - ZEROCNK;
	    wot2 = wobf( t2 );
	    woe2 = wobf( e2 );
	    e2 = e2 + wot2 - woe2 - theta_sat;
	    eor = e2 * rate;
	}
	return t2 - eor;
}


float wetlift(float pressure, float temperature, float lifted_pressure) {

#ifndef NO_QC
	if ((temperature == MISSING) || (pressure == MISSING)
                           || (lifted_pressure = MISSING)) {
	  return MISSING;
    }
#endif

    // parcels potential temperature
	float pcl_theta = theta(pressure, temperature, 1000.0);
    // some Wobus voodo
	float woth = wobf(pcl_theta);
	float wott = wobf(temperature);
    // This is the wet bulb potential temperature
	float pcl_thetaw = pcl_theta - woth + wott;
    // get the temperature that crosses the moist adiabat at
    // this pressure level
	return saturated_lift(lifted_pressure, pcl_thetaw);
}


void drylift(float pressure, float temperature, float dewpoint,
             float& pressure_at_lcl, float& temperature_at_lcl) {

    pressure_at_lcl    = MISSING;
    temperature_at_lcl = MISSING;

#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) 
                              || (dewpoint == MISSING)) {
        return;
    }
#endif

    // theta is constant from parcel level to LCL
    float pcl_theta = theta(pressure, temperature, 1000.0);

    temperature_at_lcl = lcl_temperature(temperature, dewpoint);
    pressure_at_lcl = theta_level(pcl_theta, temperature_at_lcl);

    return;

}


float lifted(float pressure, float temperature, float dewpoint,
             float lifted_pressure) {

#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING)
                                 || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl, temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, lifted_pressure);
}


float wetbulb(float pressure, float temperature, float dewpoint) {

#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) 
                                 || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl, temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, pressure);
}


float theta_wetbulb(float pressure, float temperature, float dewpoint) {

#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING)
                                 || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl, temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, 1000.0);
}


float thetae(float pressure, float temperature, float dewpoint) {

#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING)
                                 || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl, temperature_at_lcl);
    // Lift a saturated parcel to 100 mb
    float lifted_temperature = wetlift(pressure_at_lcl, temperature_at_lcl, 100.0);
    // Return the potential temperature of the 100 hPa value
    return theta(100.0, lifted_temperature, 1000.0);

}

}
