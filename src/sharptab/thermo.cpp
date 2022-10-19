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

inline float vapor_pressure(float temperature) {
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

inline float lcl_temperature(float temperature, float dewpoint) {
	float s, t, dlt;

    if ((temperature == MISSING) || (dewpoint == MISSING)) {
       return MISSING;
    } 

	s = temperature - dewpoint;
	dlt = s * (1.2185 + .001278 * temperature + s * 
          (-.00219 + 1.173E-05 * s - .0000052 * temperature));
	return temperature - dlt;
}

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

inline float theta_level(float potential_temperature, float temperature) {
    if ((potential_temperature == MISSING)
               || (temperature == MISSING)) {
       return MISSING;
    } 

	potential_temperature += ZEROCNK;
	temperature           += ZEROCNK;
	return 1000.0 / std::pow((potential_temperature/temperature), (1.0/ROCP));
}


inline float theta(float pressure, float temperature, float ref_pressure) {
    if ((temperature == MISSING) || (pressure == MISSING)
                             || (ref_pressure == MISSING)) {
       return MISSING;
    } 

	temperature += ZEROCNK;
	return (temperature * std::pow(ref_pressure / pressure, ROCP)) - ZEROCNK;
}

}
