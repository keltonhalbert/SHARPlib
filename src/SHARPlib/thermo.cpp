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

namespace sharp {

// To-Do: Can I make this branchless?
float wobf(float temperature) noexcept {
#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif
    float pol;
    float x = temperature - ZEROCNK - 20.0f;
    if (x <= 0.0f) {
       pol = 1.0f + x * (-8.841660499999999e-03f + x * ( 1.4714143e-04f \
                 + x * (-9.671989000000001e-07f + x * (-3.2607217e-08f \
                 + x * (-3.8598073e-10f)))));
       pol = pol * pol;
       return (15.13f / (pol * pol))+ZEROCNK;
    }
    else {
       pol = x * (4.9618922e-07f + x * (-6.1059365e-09f + \
             x * (3.9401551e-11f + x * (-1.2588129e-13f + \
             x * (1.6688280e-16f)))));
       pol = 1.0f + x * (3.6182989e-03f + x * (-1.3603273e-05f + pol));
       pol = pol * pol;
       return (29.93f / (pol * pol) + 0.96f * x - 14.8f)+ZEROCNK;
    }
}

// To-Do: Implement Bolton (1980) formula
float vapor_pressure(float temperature) noexcept {
#ifndef NO_QC
    if (temperature == MISSING) return MISSING;
#endif
    float pol;
    float tmpc = temperature - ZEROCNK;

    pol = tmpc * (1.1112018e-17f + tmpc * (-3.0994571e-20f));
    pol = tmpc * (2.1874425e-13f + tmpc * (-1.789232e-15f + pol));
    pol = tmpc * (4.3884180e-09f + tmpc * (-2.988388e-11f + pol));
    pol = tmpc * (7.8736169e-05f + tmpc * (-6.111796e-07f + pol));
    pol = .99999683e-00f + tmpc * (-9.082695e-03f + pol);
    pol = (pol * pol);
    pol = (pol * pol);
    return (610.78f / (pol * pol)); 
}

float lcl_temperature(float temperature, float dewpoint) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (dewpoint == MISSING)) {
        return MISSING;
    }
#endif
    float s, dlt;
    float tmpc = temperature - ZEROCNK;

    s = temperature - dewpoint;
    dlt = s * (1.2185f + 0.001278f * tmpc +
          s * (-0.00219f + 1.173E-05f * s - 0.0000052f * tmpc));
    return (tmpc - dlt) + ZEROCNK;
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
    float pres = pressure * PA_TO_HPA;

    double x = std::log10(wv_mixratio * pres / (EPSILON + wv_mixratio));
    double tmrk = std::pow(10.0, c1 * x + c2) - c3 +
                  c4 * std::pow(std::pow(10.0, c5 * x) - c6, 2.0);

    return (float)tmrk;
}

float theta_level(float potential_temperature, float temperature) noexcept {
#ifndef NO_QC
    if ((potential_temperature == MISSING) || (temperature == MISSING)) {
        return MISSING;
    }
#endif
    return THETA_REF_PRESSURE /
           std::pow((potential_temperature / temperature), (1.0f / ROCP));
}

float theta(float pressure, float temperature, float ref_pressure) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING) ||
        (ref_pressure == MISSING)) {
        return MISSING;
    }
#endif
    return (temperature * std::pow(ref_pressure / pressure, ROCP));
}

// To-Do: Implement exact version
float mixratio(float pressure, float temperature) noexcept {
#ifndef NO_QC
    if ((temperature == MISSING) || (pressure == MISSING)) {
        return MISSING;
    }
#endif

    // correction factor for the departure of the mixture
    // of air and water vapor from the ideal gas law
    float tmpc = temperature - ZEROCNK;
    float pres = pressure * PA_TO_HPA;

    float x = 0.02f * (tmpc - 12.5f + 7500.0f / pres);
    float wfw = 1.0f + 0.0000045f * pres + 0.0014f * x * x;
    float fwesw = wfw * vapor_pressure(temperature) * PA_TO_HPA;
    return EPSILON * (fwesw / (pres - fwesw));
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
    float w = mixratio(pressure, dewpoint);
    return (temperature * ((1.0f + (w / EPSILON)) / (1.0f + w)));
}

float saturated_lift(float pressure, float theta_sat) noexcept {
#ifndef NO_QC
    if ((pressure == MISSING) || (theta_sat == MISSING)) {
        return MISSING;
    }
#endif

    float iter_step = 0.001f;
    if ((std::fabs(pressure - THETA_REF_PRESSURE) - iter_step) <= 0.0f)
        return theta_sat;

    float pwrp = std::pow(pressure / THETA_REF_PRESSURE, ROCP);
    // get the temperature
    float t1 = theta_sat * pwrp;
    float e1 = wobf(t1) - wobf(theta_sat);
    float rate = 1.0f;
    float eor = 999;
    float t2;
	int condition = 0;
	// Testing the original showed that only
	// 5 or so iterations are needed, but
	// double that just in case. It'll exit
	// early if it converges anyway. 
	for (int iter = 0; iter < 10; ++iter) {
        t2 = t1 - e1 * rate;
        float e2 = (t2) / pwrp;
        e2 = e2 + wobf(t2) - wobf(e2) - theta_sat;

        eor = e2 * rate;
        rate = (t2 - t1) / (e2 - e1);
        t1 = t2;
		e1 = e2;
		condition |= ((std::fabs(eor) - iter_step) < 0.0f);
		if (condition) break;
    }
    return t2 - eor;
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
    float pcl_theta = theta(pressure, temperature, THETA_REF_PRESSURE);
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
    float pcl_theta = theta(pressure, temperature, THETA_REF_PRESSURE);

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
    return wetlift(pressure_at_lcl, temperature_at_lcl, THETA_REF_PRESSURE);
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
    constexpr float lift_top = 10000.0f; // 100 hPa but its Pa

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pressure_at_lcl,
            temperature_at_lcl);
    // Lift a saturated parcel to 100 mb
    float lifted_temperature =
        wetlift(pressure_at_lcl, temperature_at_lcl, lift_top);
    // Return the potential temperature of the 100 hPa value
    return theta(lift_top, lifted_temperature, THETA_REF_PRESSURE);
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
           (env_temperature);
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


