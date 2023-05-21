// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC.
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#include <cmath>
#include <iostream>

namespace sharp {

float wobf(float temperature) noexcept {
#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif
    float pol;
    float x = temperature - 20.0f;
    if (x <= 0.0f) {
       pol = 1.0f + x * (-8.841660499999999e-03f + x * ( 1.4714143e-04f \
                 + x * (-9.671989000000001e-07f + x * (-3.2607217e-08f \
                 + x * (-3.8598073e-10f)))));
       pol = pol * pol;
       return 15.13f / (pol * pol);
    }
    else {
       pol = x * (4.9618922e-07f + x * (-6.1059365e-09f + \
             x * (3.9401551e-11f + x * (-1.2588129e-13f + \
             x * (1.6688280e-16f)))));
       pol = 1.0f + x * (3.6182989e-03f + x * (-1.3603273e-05f + pol));
       pol = pol * pol;
       return 29.93f / (pol * pol) + 0.96f * x - 14.8f;
    }
}

float vapor_pressure(float temperature) noexcept {
#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif
    float pol;

    pol = temperature * (1.1112018e-17f + temperature * (-3.0994571e-20f));
    pol = temperature * (2.1874425e-13f + temperature * (-1.789232e-15f + pol));
    pol = temperature * (4.3884180e-09f + temperature * (-2.988388e-11f + pol));
    pol = temperature * (7.8736169e-05f + temperature * (-6.111796e-07f + pol));
    pol = .99999683e-00f + temperature * (-9.082695e-03f + pol);
    pol = (pol * pol);
    pol = (pol * pol);
    return 6.1078f / (pol * pol);
}

float lcl_temperature(float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif
    float s, dlt;

    s = temperature - dewpoint;
    dlt = s * (1.2185f + 0.001278f * temperature +
          s * (-0.00219f + 1.173E-05f * s - 0.0000052f * temperature));
    return temperature - dlt;
}

float temperature_at_mixratio(float wv_mixratio, float pressure) noexcept {
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

    double x = std::log10(wv_mixratio * pressure / (622.0 + wv_mixratio));
    double tmrk = std::pow(10.0, c1 * x + c2) - c3 +
                  c4 * std::pow(std::pow(10.0, c5 * x) - c6, 2.0);

    return (float)(tmrk - ZEROCNK);
}

float theta_level(float potential_temperature, float temperature) noexcept {
#ifndef NO_QC
    if ((potential_temperature == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif

    // convert to degrees kelvin
    potential_temperature += ZEROCNK;
    temperature += ZEROCNK;
    return 1000.0f /
           std::pow((potential_temperature / temperature), (1.0f / ROCP));
}

float theta(float pressure, float temperature, float ref_pressure) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING) ||
        (ref_pressure == MISSING)) {
        return MISSING;
    }
#endif

    // get temperature in kelvin
    temperature += ZEROCNK;
    return (temperature * std::pow(ref_pressure / pressure, ROCP)) - ZEROCNK;
}

float mixratio(float pressure, float temperature) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif

    // correction factor for the departure of the mixture
    // of air and water vapor from the ideal gas law
    float x = 0.02f * (temperature - 12.5f + 7500.0f / pressure);
    float wfw = 1.0f + 0.0000045f * pressure + 0.0014f * x * x;
    float fwesw = wfw * vapor_pressure(temperature);
    return 621.97f * (fwesw / (pressure - fwesw));
}

float virtual_temperature(float pressure, float temperature,
                          float dewpoint) noexcept {
#ifndef NO_QC
    if (dewpoint == MISSING) {
        return temperature;
    } else if ((pressure == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif

    constexpr float eps = 0.62197f;
    float temperature_k = temperature + ZEROCNK;
    float w = 0.001f * mixratio(pressure, dewpoint);
    return (temperature_k * ((1.0f + (w / eps)) / (1.0f + w))) - ZEROCNK;
}

float saturated_lift(float pressure, float theta_sat) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (theta_sat == MISSING)) {
        return MISSING;
    }
#endif

    // we do not want to go below the 1000.0 hPa reference level
    if ((std::fabs(pressure - 1000.0f) - 0.001f) <= 0.0f) return theta_sat;

    float eps = 0.001f;
    float pwrp = std::pow(pressure / 1000.0f, ROCP);
    // get the temperature
    float t1 = (theta_sat + ZEROCNK) * pwrp - ZEROCNK;
    float e1 = wobf(t1) - wobf(theta_sat);
    float rate = 1.0f;
    float eor = 999;
    float t2;
	// Testing the original showed that only
	// 5 or so iterations are needed, but
	// double that just in case. It'll exit
	// early if it converges anyway. 
	for (int iter = 0; iter < 10; ++iter) {
        t2 = t1 - e1 * rate;
        float e2 = (t2 + ZEROCNK) / pwrp - ZEROCNK;
        e2 = e2 + wobf(t2) - wobf(e2) - theta_sat;

        eor = e2 * rate;
        rate = (t2 - t1) / (e2 - e1);
        t1 = t2;
		e1 = e2;

        if (std::fabs(eor) - eps < 0.0f) {
            return t2 - eor;	
        }
    }
    return t2;
}

float wetlift(float pressure, float temperature,
              float lifted_pressure) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING) ||
        (lifted_pressure == MISSING)) {
        return MISSING;
    }
#endif

    // parcels potential temperature
    float pcl_theta = theta(pressure, temperature, 1000.0f);
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
             float& pressure_at_lcl, float& temperature_at_lcl) noexcept {
    // we do this before the QC check so that
    // these values are passed back as missing
    pressure_at_lcl = MISSING;
    temperature_at_lcl = MISSING;
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return;
    }
#endif

    // theta is constant from parcel level to LCL
    float pcl_theta = theta(pressure, temperature, 1000.0f);

    temperature_at_lcl = lcl_temperature(temperature, dewpoint);
    pressure_at_lcl = theta_level(pcl_theta, temperature_at_lcl);

	if (pressure_at_lcl > pressure) pressure_at_lcl = pressure;

    return;
}

float lifted(float pressure, float temperature, float dewpoint,
             float lifted_pressure) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, lifted_pressure);
}

float wetbulb(float pressure, float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, pressure);
}

float theta_wetbulb(float pressure, float temperature,
                    float dewpoint) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    return wetlift(pressure_at_lcl, temperature_at_lcl, 1000.0f);
}

float thetae(float pressure, float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pressure_at_lcl = MISSING;
    float temperature_at_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    // Lift a saturated parcel to 100 mb
    float lifted_temperature =
        wetlift(pressure_at_lcl, temperature_at_lcl, 100.0f);
    // Return the potential temperature of the 100 hPa value
    return theta(100.0f, lifted_temperature, 1000.0f);
}

float lapse_rate(HeightLayer layer_agl, const float* height,
                 const float* temperature, int num_levs) noexcept {
#ifndef NO_QC
    if ((layer_agl.bottom == MISSING) || (layer_agl.top == MISSING)) {
        return MISSING;
    }
#endif

    // convert from agl to msl
    layer_agl.bottom += height[0];
    layer_agl.top += height[0];

    // bounds check the height layer
    if (layer_agl.bottom < height[0]) {
        layer_agl.bottom = height[0];
    }
    if (layer_agl.top > height[num_levs - 1]) {
        layer_agl.top = height[num_levs - 1];
    }

    // lower and upper temperature
    float tmpc_l =
        interp_height(layer_agl.bottom, height, temperature, num_levs);
    float tmpc_u = interp_height(layer_agl.top, height, temperature, num_levs);
#ifndef NO_QC
    if ((tmpc_l == MISSING) || (tmpc_u == MISSING)) {
        return MISSING;
    }
#endif

    // dT/dz, positive (definition of lapse rate), in km
    float dz = layer_agl.top - layer_agl.bottom;
    return ((tmpc_u - tmpc_l) / dz) * -1000.0f;
}

float lapse_rate(PressureLayer layer, const float* pressure,
                 const float* height, const float* temperature,
                 int num_levs) noexcept {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    // bounds check the pressure layer
    if (layer.bottom > pressure[0]) {
        layer.bottom = pressure[0];
    }
    if (layer.top < pressure[num_levs - 1]) {
        layer.top = pressure[num_levs - 1];
    }

    HeightLayer h_layer =
        pressure_layer_to_height(layer, pressure, height, num_levs, true);

    return lapse_rate(h_layer, height, temperature, num_levs);
}

float buoyancy(float pcl_temperature, float env_temperature) noexcept {
    return GRAVITY * (pcl_temperature - env_temperature) /
           (env_temperature + ZEROCNK);
}

float moist_static_energy(float height_agl, float temperature,
                          float specific_humidity) noexcept {
#ifndef NO_QC
    if ((height_agl == MISSING) || (temperature == MISSING) ||
        (specific_humidity == MISSING)) {
        return MISSING;
    }
#endif

    return (CP_DRYAIR * temperature) + (LV * specific_humidity) +
           (GRAVITY * height_agl);
}

}  // end namespace sharp

namespace sharp::exper {


} // end namespace sharp::exper


