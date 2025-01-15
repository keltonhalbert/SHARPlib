/**
 * \file
 * \brief Routines used to compute thermodynamic attributes of vertical<!--
 * --> sounding profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef SHARP_THERMO_H
#define SHARP_THERMO_H

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#include <cmath>
#include <cstddef>

namespace sharp {

/**
 * \brief Defines the kinds of moist adiabats used by sharp::moist_adiabat_cm1
 */
enum class adiabat : int {
    /**
     * \brief Pseudoadiabatic considering only liquid water
     */
    pseudo_liq = 1,

    /**
     * \brief Adiabatic considering only liquid water
     */
    adiab_liq = 2,

    /**
     * \brief Pseudoadiabatic considering liquid water + ice
     */
    pseudo_ice = 3,

    /**
     * \brief Adiabatic considering liquid water + ice
     */
    adiab_ice = 4,
    END,
};

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the difference between the wet-bulb potential<!--
 * --> temperatures for saturated and dry air given the temperature.
 *
 * The Wobus Function (wobf) is defined as the difference between
 * the wet-bulb potential temperature for saturated air (WBPTS)
 * and the wet-bulb potential temperature for dry air (WBPTD) given
 * an input temperature in Kelvin.
 *
 * WOBF(T) := WBPTS - WBPTD
 *
 * Although WBPTS and WBPTD are functions of both pressure and
 * temperature, it is assumed their difference is a function of
 * temperature only. The difference is also proportional to the
 * heat imparted to a parcel.
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
 * Kelvin in the temperature of a lifted parcel.
 *
 * \param   temperature     (K)
 *
 * \return  wobf            (K)
 */
[[nodiscard]] float wobf(float temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the vapor pressure over liquid water
 *
 * Computes the vapor pressure (or saturation vapor pressure) in
 * Pascals (Pa) over liquid water given the temperature in
 * Kelvin (or dewpoint temperature in Kelvin).
 * The air pressure is used as a minimum floor for extremely
 * cold temperatures and low pressures.
 *
 * This function uses the formulation by Bolton (1980), and is
 * accurate to 0.3% for the temperature range of -35C <= T <= 35 C
 *
 * \param   pressure       (Pa)
 * \param   temperature    (K)
 *
 * \return  vapor_pressure (Pa)
 */
[[nodiscard]] float vapor_pressure(float pressure, float temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the vapor pressure over ice
 *
 * Computes the vapor pressure (or saturation vapor pressure) in
 * Pascals (Pa) over ice given the temperature in
 * Kelvin (or dewpoint temperature in Kelvin).
 * The air pressure is used as a minimum floor for extremely
 * cold temperatures and low pressures.
 *
 * This function uses the formulation by Bolton (1980), and is
 * accurate to 0.3% for the temperature range of -35C <= T <= 35C
 *
 * \param   pressure       (Pa)
 * \param   temperature    (K)
 *
 * \return  vapor_pressure (Pa)
 */
[[nodiscard]] float vapor_pressure_ice(float pressure, float temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the temperature of the LCL.
 *
 * Computes the temperature of a parcels LCL in
 * degres Kelvin given the parcel temperature and
 * dewpoint in Kelvin.
 *
 * This is implemented as in Bolton (1980) eq 15, and is considered
 * to be within a 10th of a degree of the more exact iterative formula.
 *
 * \param    temperature     (K)
 * \param    dewpoint        (K)
 *
 * \return   lcl_temperature (K)
 */
[[nodiscard]] float lcl_temperature(const float temperature,
                                    const float dewpoint);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the temperature at a given water vapor mixing ratio<!--
 * --> and pressure level.
 *
 * Computes the temperature in Kelvin of air at the
 * given water vapor mixing ratio in kg/kg and the
 * air pressure in Pa.
 *
 * The is implemented as in Bolton (1980) eq 11, and is considered
 * to be accurate to 0.03 for -35C <= T <= 35C
 *
 * \param    wv_mixratio    (kg/kg)
 * \param    pressure       (Pa)
 *
 * \return   temperature    (K)
 */
[[nodiscard]] float temperature_at_mixratio(float wv_mixratio, float pressure);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Comute the pressure level given potential temperature and temperature.
 *
 * Returns the pressure level in Pascals (Pa) of a parcel
 * given the potential temperature in Kelvin and the
 * temperature of the parcel in Kelvin.
 *
 * \param    potential_temperature (K)
 * \param    temperature           (K)
 *
 * \return   pressure              (Pa)
 */
[[nodiscard]] float theta_level(float potential_temperature, float temperature);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the potential temperature.
 *
 * Returns the potential temperature in Kelvin of
 * a parcel given its pressure in Pascals and temperature
 * in Kelvin. The final argument is the reference
 * level, which is usually 100000.0 Pa, conveniently also
 * called sharp::THETA_REF_PRESSURE.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (K)
 * \param    ref_pressure          (Pa)
 *
 * \return   potential_temperature (K)
 */
[[nodiscard]] float theta(float pressure, float temperature,
                          float ref_pressure = THETA_REF_PRESSURE);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the water vapor mixing ratio from specific humidity.
 *
 * Returns the water vapor mixing ratio in kg/kg given the specific
 * humidity in kg/kg.
 *
 * \param	 q            (kg/kg)
 *
 * \return   mixratio   (kg/kg)
 */
[[nodiscard]] float mixratio(float q);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the water vapor mixing ratio.
 *
 * Returns the water vapor mixing ratio in kg/kg given
 * the environmental pressure in Pascals (Pa) and a
 * temperature (dry-bulb or dewpoint) in Kelvin.
 *
 * This is computed by calling sharp::vapor_pressure, which uses
 * the Bolton (1980) equations.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (K)
 *
 * \return   mixratio              (kg/kg)
 */
[[nodiscard]] float mixratio(float pressure, float temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the ice water mixing ratio.
 *
 * Returns the ice water mixing ratio in kg/kg given
 * the environmental pressure in Pascals (Pa) and a
 * temperature (dry-bulb or dewpoint) in Kelvin.
 *
 * This is computed by calling sharp::vapor_pressure_ice, which uses
 * the Bolton (1980) equations.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (K)
 *
 * \return   mixratio              (kg/kg)
 */
[[nodiscard]] float mixratio_ice(float pressure, float temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute specific humidity from a mixing ratio
 *
 * Returns the specific humidity given the mixing ratio.
 *
 * \param   rv                  (kg/kg)
 *
 * \return  specific_humidity   (unitless)
 */
[[nodiscard]] float specific_humidity(float rv);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the full virtual temperature.
 *
 * Returns the virtual temperature in Kelvin given the dry-bulb
 * temperature in Kelvin, the specific humidity of water vapor (qv) in
 * kg/kg, the specific humidity of liquid water (ql) in kg/kg, and the specific
 * himidity of ice water kg/kg.
 *
 * For convenience, the ql and qi terms have a default value of zero to
 * easily support returning just the virtual temperature from water vapor.
 * If you are okay with lower accuracy, you may pass mixing ratios intead
 * of specific humidities.
 *
 * \param   temperature             (K)
 * \param   qv                      (kg/kg)
 * \param   ql                      (kg/kg)
 * \param   qi                      (kg/kg)
 *
 * \return  virtual_temperature     (K)
 */
[[nodiscard]] float virtual_temperature(float temperature, float qv,
                                        float ql = 0.0f, float qi = 0.0f);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the temperature along a moist adiabat <!--
 * -->(wet-bulb potential temperature) at a given pressure
 *
 * Compute the temperature at which the moist adiabat intersects a line
 * of constant pressure on a Skew-T log-P diagram. The wet-bulb potential
 * temperature, given by theta_sat, defines a moist adiabat in Kelvin,
 * and the temperature at the given pressure level in Pascals (Pa) is returned.
 *
 * This function relies on the Wobus Function ( sharp::wobf ), and it was shown
 * by Robert Davies-Jones (2007) that the Wobus function has a slight
 * dependence on pressure, which results in errors of up to 1.2 Kelvin
 * in the temperature of a lifted parcel.
 *
 * \param   pressure            (Pa)
 * \param   theta_sat           (K)
 * \param   converge            (convergence criteria; default = 0.001f)
 *
 * \return  lifted_temperature  (K)
 */
[[nodiscard]] float saturated_lift(float pressure, float theta_sat,
                                   const float converge = 0.001f);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the temperature of a parcel lifted moist <!--
 * -->adiabatically to a new level.
 *
 * With a given parcel defined by a pressure and temperature (in Pascals and
 * Kelvin), lift it moist adiabatically to a new pressure level
 * (in Pascals) and return the temperture of the parcel at that level.
 *
 * This function relies on the Wobus Function ( sharp::wobf ), and it was shown
 * by Robert Davies-Jones (2007) that the Wobus function has a slight
 * dependence on pressure, which results in errors of up to 1.2 Kelvin
 * in the temperature of a lifted parcel.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (K)
 * \param    lifted_pressure       (Pa)
 *
 * \return   lifted_temperature    (K)
 */
[[nodiscard]] float wetlift(float pressure, float temperature,
                            float lifted_pressure);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute the temperature of a parcel lifted moist <!--
 * -->adiabatically to a new level.
 *
 *  With a given parcel defined by a pressure and temperature (in Pascals and
 *  Kelvin), lift it moist adiabatically to a new pressure level
 *  (in Pascals) and return the temperature of the parcel at that level.
 *
 *  This function is based on/ripped from George Byran's (NCAR) routine in CM1.
 *  It has several options for convergence criteria, pressure increments, and
 *  the type of moist adiabat/assumptions (i.e. pseudoadiabatic, adiabatic,
 *  liquid, liquid+ice). The sharp::adiabat defines the types of adiabats,
 *  and the float references to qv, ql, and qi are used to keep budgets of the
 *  moisture variables.
 *
 *  NOTE: rv_total should most likely be the water vapor mixing ratio either
 *  at the parcel's sharp::LPL, or the water vapor mixing ratio at the LCL
 *  (these are going to be the same value). Essentially, the total water vapor
 *  before condensation.
 *
 * \param   pressure        (Pa)
 * \param   temperature     (K)
 * \param   new_pressure    (Pa)
 * \param   rv_total        (kg/kg)
 * \param   rv              (kg/kg)
 * \param   rl              (kg/kg)
 * \param   ri              (kg/kg)
 * \param   pres_incr       (Pa)
 * \param   converge        (precision)
 * \param   ma_type         (sharp::adiabat)
 *
 * \return  pcl_temperature (K)
 */
[[nodiscard]] float moist_adiabat_cm1(float pressure, float temperature,
                                      float new_pressure, float& rv_total,
                                      float& rv, float& rl, float& ri,
                                      const float pres_incr,
                                      const float converge,
                                      const adiabat ma_type);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Lift a parcel dry adiabatically to its Lifted Condensation Level
 * (LCL).
 *
 * Given a parcel's initial pressure (Pa), temperature (K), and
 * dewpoint temperature (K), lift the parcel dry adiabatically to its
 * Lifted Condensation Level and store the resulting LCL pressure (Pa)
 * and LCL temperature (K) in the variables passed by reference to the
 * routine.
 *
 *
 * The LCL temperature is computed using an approximation. See the
 * sharp::lcl_temperature documentation for more information.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (K)
 * \param    dewpoint              (K)
 * \param    pressure_at_lcl       (Pa)
 *
 * \param    temperature_at_lcl    (K)
 */
void drylift(float pressure, float temperature, float dewpoint,
             float& pressure_at_lcl, float& temperature_at_lcl);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the wetbulb temperature.
 *
 * Compute the wet bulb temperature (K) given the pressure
 * (Pa), temperature (K), and dewpoint (K).
 *
 * First, it lifts a parcel with the given pressure, temperature, and
 * dewpoint temperature to its Lifted Condensation Level (LCL). To
 * compute the temperature and pressure of the LCL, an approximation
 * is used. See the sharp::lcl_temperature documentation for more
 * information.
 *
 * After the parcel has reached the LCL, the lifter passed in
 * lowers the parcel to its initial pressure level along a moist adiabat.
 *
 * \param   lifter                  (cm1, wobf)
 * \param   pressure                (Pa)
 * \param   temperature             (K)
 * \param   dewpoint                (K)
 *
 * \return  wetbulb_temperature     (K)
 */
template <typename Lft>
[[nodiscard]] float wetbulb(Lft lifter, float pressure, float temperature,
                            float dewpoint) {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pres_lcl = MISSING;
    float tmpk_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pres_lcl, tmpk_lcl);
    lifter.setup(pres_lcl, tmpk_lcl);

    return lifter(pres_lcl, tmpk_lcl, pressure);
}

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the wetbulb potential temperature.
 *
 * Compute the wet bulb potential temperature (K) given
 * the pressure (Pa), temperature (K), and dewpoint
 * (K).
 *
 * First, it lifts a parcel with the given pressure, temperature, and
 * dewpoint temperature to its Lifted Condensation Level (LCL). To
 * compute the temperature and pressure of the LCL, an approximation
 * is used. See the sharp::lcl_temperature documentation for more
 * information.
 *
 * After the parcel has reached the LCL, the lifter passed in
 * lowers the parcel to the standard reference pressure level
 * (1000.0 hPa) along a moist adiabat.
 *
 * \param   lifter                          (wobf, cm1)
 * \param   pressure                        (Pa)
 * \param   temperature                     (K)
 * \param   dewpoint                        (K)
 *
 * \return  wetbulb_potential_temperature   (K)
 */
template <typename Lft>
[[nodiscard]] float theta_wetbulb(Lft lifter, float pressure, float temperature,
                                  float dewpoint) {
#ifndef NO_QC
    if ((pressure == MISSING) || (temperature == MISSING) ||
        (dewpoint == MISSING)) {
        return MISSING;
    }
#endif

    float pres_lcl = MISSING;
    float tmpk_lcl = MISSING;

    // pressure_at_lcl and temperature_at_lcl are passed by reference,
    // so the values are changed by the drylift routine
    drylift(pressure, temperature, dewpoint, pres_lcl, tmpk_lcl);
    lifter.setup(pres_lcl, tmpk_lcl);

    return lifter(pres_lcl, tmpk_lcl, THETA_REF_PRESSURE);
}

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the equivalent potential temperature.
 *
 * Compute the equivalent potential temperature (K) given
 * the pressure (Pa), temperature (K), and dewpoint
 * (K). The equivalent potential temperature is computed as in
 * Bolton 1980. It was found in Davies-Jones 2009 that this is
 * the most accurate non-iterative formulation of ThetaE.
 *
 * Bolton 1980:
 * https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2
 * Davies-Jones 2009: https://doi.org/10.1175/2009MWR2774.1
 *
 * \param    pressure                           (Pa)
 * \param    temperature                        (K)
 * \param    dewpoint                           (K)
 *
 * \return   equivalent_potential_temperature   (K)
 */
[[nodiscard]] float thetae(float pressure, float temperature, float dewpoint);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief compute the lapse rate over the given sharp::HeightLayer (AGL)
 *
 * Computes the lapse rate over a given sharp::HeightLayer above-ground-level.
 * This routine handles converting AGL to MSL by adding the surface height
 * value to the layer.
 *
 * \param   layer_agl               (meters AGL)
 * \param   height                  (meters MSL)
 * \param   temperature             (K)
 * \param   N                       (length of arrays)
 *
 * \return  Temperature Lapse Rate  (K/km)
 */
[[nodiscard]] float lapse_rate(HeightLayer layer_agl, const float height[],
                               const float temperature[],
                               const std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief compute the lapse rate over the given sharp::PressureLayer
 *
 * Computes the lapse rate over a given sharp::PressureLayer. This routine
 * handles converting the pressure layer into a height layer, and
 * then calles the sharp::HeightLayer implementation of this routine.
 *
 * \param   layer                   (Pa)
 * \param   pressure                (Pa)
 * \param   height                  (meters MSL)
 * \param   temperature             (K)
 * \param   N		                    (length of arrays)
 *
 * \return  Temperature Lapse Rate  (K/km)
 */
[[nodiscard]] float lapse_rate(PressureLayer layer, const float pressure[],
                               const float height[], const float temperature[],
                               const std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Find the maximum lapse rate within a specified layer/depth (AGL).
 *
 * Given a layer of the atmosphere (e.g. 2 - 6 km), find the maximum lapse rate
 * over the provided depth (e.g. 2 km) within that given layer. Returns the
 * maximum lapse rate, with the optional ability to return the layer that lapse
 * rate exists within.
 *
 * \param   layer_agl       (meters AGL)
 * \param   depth           (meters)
 * \param   height          (meters MSL)
 * \param   temperature     (K)
 * \param   N               (length of arrays)
 * \param   max_lyr         (The sharp::HeightLayer of the max lapse rate)
 *
 * \return  Maximum Temperature Lapse Rate  (K/km)
 */
[[nodiscard]] float lapse_rate_max(HeightLayer layer_agl, const float depth,
                                   const float height[],
                                   const float temperature[],
                                   const std::ptrdiff_t N,
                                   HeightLayer* max_lyr = nullptr);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Find the maximum lapse rate within a specified layer/depth (Pres)
 *
 * Given a mayer of the atmisphere (e.g. 800 hPa - 500 hPa), find the maximum
 * lapse rate over the provided depth (e.g. 100 hPa) within that given layer.
 * Returns the maximum lapse rate, with the optional ability to return the layer
 * that lapse rate exists within.
 *
 * \param   layer           (Pa)
 * \param   depth           (Pa)
 * \param   pressure        (Pa)
 * \param   height          (meters MSL)
 * \param   temperature     (K)
 * \param   N		            (length of arrays)
 * \param   max_lyr         (The sharp::PressureLayer of the max lapse rate)
 *
 * \return  Maximum Temperature Lapse Rate  (K/km)
 */
[[nodiscard]] float lapse_rate_max(PressureLayer layer, const float depth,
                                   const float pressure[], const float height[],
                                   const float temperature[],
                                   const std::ptrdiff_t N,
                                   PressureLayer* max_lyr = nullptr);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief compute buoyancy given parcel and environment temperatures
 *
 * \param   pcl_temperature     (K)
 * \param   env_temperature	    (K)
 *
 * \return  buoyancy            (m/s^2)
 */
[[nodiscard]] float buoyancy(const float pcl_temperature,
                             const float env_temperature);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute buoyancy given arrays of parcel & environment temperatures
 *
 * \param pcl_temperature   (K)
 * \param env_temperature   (K)
 * \param buoy_arr          Output array to fill with buoyancy
 * \param N                 Length of arrays
 */
void buoyancy(const float pcl_temperature[], const float env_temperature[],
              float buoy_arr[], std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Compute moist static energy.
 *
 * \param   height_agl          (meters)
 * \param   temperature         (K)
 * \param   specific_humidity   (kg/kg)
 *
 * \return  moist static energy ()
 */
[[nodiscard]] float moist_static_energy(float height_agl, float temperature,
                                        float specific_humidity);

[[nodiscard]] float buoyancy_dilution_potential(const float temperature,
                                                const float mse_bar,
                                                const float saturation_mse);

}  // end namespace sharp

#endif  // SHARP_THERMP_H
