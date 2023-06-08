/**
 * \file
 * \brief Routines used to compute thermodynamic attributes of vertical<!--
 * --> sounding profiles
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
#ifndef __SHARP_THERMO_H__
#define __SHARP_THERMO_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

#include <algorithm>
#include <cmath>

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
 * the same temperature in Celsius.
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
 * degrees Kelvin in the temperature of a lifted parcel.
 *
 * \param   temperature     (degK)
 *
 * \return  wobf            (degK)
 */
[[nodiscard]] float wobf(float temperature) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the vapor pressure over liquid water
 *
 * Computes the vapor pressure (or saturation vapor pressure) in
 * Pascals (Pa) over liquid water given the temperature in 
 * degrees Kelvin (or dewpoint temperature in degrees Kelvin).
 * The air pressure is used as a minimum floor for extremely
 * cold temperatures and low pressures.
 *
 * This function uses the formulation by Bolton (1980), and is
 * accurate to 0.3% for the temperature range of -35C <= T <= 35 C
 *
 * \param   pressure       (Pa)
 * \param   temperature    (degK)
 *
 * \return  vapor_pressure (Pa)
 */
[[nodiscard]] float vapor_pressure(float pressure, float temperature) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the vapor pressure over ice 
 *
 * Computes the vapor pressure (or saturation vapor pressure) in
 * Pascals (Pa) over ice given the temperature in 
 * degrees Kelvin (or dewpoint temperature in degrees Kelvin).
 * The air pressure is used as a minimum floor for extremely
 * cold temperatures and low pressures.
 *
 * This function uses the formulation by Bolton (1980), and is
 * accurate to 0.3% for the temperature range of -35C <= T <= 35C
 *
 * \param   pressure       (Pa)
 * \param   temperature    (degK)
 *
 * \return  vapor_pressure (Pa)
 */
[[nodiscard]] float vapor_pressure_ice(float pressure,
                                       float temperature) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the temperature of the LCL.
 *
 * Computes the temperature of a parcels LCL in
 * degres Kelvin given the parcel temperature and
 * dewpoint in degrees Kelvin.
 *
 * This is implemented as in Bolton (1980) eq 15, and is considered
 * to be within a 10th of a degree of the more exact iterative formula.
 *
 * \param    temperature     (degK)
 * \param    dewpoint        (degK)
 *
 * \return   lcl_temperature (degK)
 */
[[nodiscard]] float lcl_temperature(float temperature, float dewpoint) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the temperature at a given water vapor mixing ratio<!--
 * --> and pressure level.
 *
 * Computes the temperature in degrees Kelvin of air at the
 * given water vapor mixing ratio in kg/kg and the
 * air pressure in Pa.
 *
 * The is implemented as in Bolton (1980) eq 11, and is considered
 * to be accurate to 0.03 for -35C <= T <= 35C
 *
 * \param    wv_mixratio    (g/kg)
 * \param    pressure       (mb)
 *
 * \return   temperature    (degC)
 */
[[nodiscard]] float temperature_at_mixratio(float wv_mixratio,
                                            float pressure) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Comute the pressure level given potential temperature and temperature.
 *
 * Returns the pressure level in Pascals (Pa) of a parcel
 * given the potential temperature in degrees Kelvin and the
 * temperature of the parcel in degrees Kelvin.
 *
 * \param    potential_temperature (degK)
 * \param    temperature           (degK)
 *
 * \return   pressure              (Pa)
 */
[[nodiscard]] float theta_level(float potential_temperature,
                                float temperature) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the potential temperature.
 *
 * Returns the potential temperature in degrees Kelvin of
 * a parcel given its pressure in Pascals and temperature 
 * in degrees Kelvin. The final argument is the reference 
 * level, which is usually 100000.0 Pa, conveniently also
 * called sharp::THETA_REF_PRESSURE.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (degK)
 * \param    ref_pressure          (Pa)
 *
 * \return   potential_temperature (degK)
 */
[[nodiscard]] float theta(float pressure, float temperature,
                          float ref_pressure=THETA_REF_PRESSURE) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the water vapor mixing ratio.
 *
 * Returns the water vapor mixing ratio in kg/kg given
 * the environmental pressure in Pascals (Pa) and a
 * temperature (dry-bulb or dewpoint) in degrees Kelvin.
 *
 * This is computed by calling sharp::vapor_pressure, which uses
 * the Bolton (1980) equations.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (degK)
 *
 * \return   mixratio              (kg/kg)
 */
[[nodiscard]] float mixratio(float pressure, float temperature) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the ice water mixing ratio.
 *
 * Returns the ice water mixing ratio in kg/kg given
 * the environmental pressure in Pascals (Pa) and a
 * temperature (dry-bulb or dewpoint) in degrees Kelvin.
 *
 * This is computed by calling sharp::vapor_pressure_ice, which uses
 * the Bolton (1980) equations.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (degK)
 *
 * \return   mixratio              (kg/kg)
 */
[[nodiscard]] float mixratio_ice(float pressure, float temperature) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute specific humidity from a mixing ratio
 *
 * Returns the specific humidity given the mixing ratio.
 *
 * \param   q                   (kg/kg)
 *
 * \return  specific_humidity   (unitless)
 */
[[nodiscard]] float specific_humidity(float q) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the full virtual temperature.
 *
 * Returns the virtual temperature in degrees Kelvin given the dry-bulb 
 * temperature in degrees Kelvin, the water vapor mixing ratio (qv) in 
 * kg/kg, the liquid water mixing ratio (ql) in kg/kg, and the ice water
 * mixing ratio in kg/kg. 
 *
 * For convenience, the ql and qi terms have a default value of zero to 
 * easily support returning just the virtual temperature from water vapor. 
 *
 * \param   temperature             (degK)
 * \param   qv                      (kg/kg)
 * \param   ql                      (kg/kg)
 * \param   qi                      (kg/kg)
 *
 * \return  virtual_temperature     (degK)
 */
[[nodiscard]] float virtual_temperature(float temperature, float qv,
                                        float ql = 0.0f,
                                        float qi = 0.0f) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the temperature along a moist adiabat <!--
 * -->(wet-bulb potential temperature) at a given pressure
 *
 * Compute the temperature at which the moist adiabat intersects a line
 * of constant pressure on a Skew-T log-P diagram. The wet-bulb potential
 * temperature, given by theta_sat, defines a moist adiabat in degrees Kelvin, 
 * and the temperature at the given pressure level in Pascals (Pa) is returned.
 *
 * This function relies on the Wobus Function ( sharp::wobf ), and it was shown
 * by Robert Davies-Jones (2007) that the Wobus function has a slight
 * dependence on pressure, which results in errors of up to 1.2 degrees Kelvin
 * in the temperature of a lifted parcel.
 *
 * \param   pressure            (Pa)
 * \param   theta_sat           (degK)
 * \param   converge            (convergence criteria; default = 0.001f)
 *
 * \return  lifted_temperature  (degK)
 */
[[nodiscard]] float saturated_lift(float pressure, float theta_sat,
                                   const float converge = 0.001f) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the temperature of a parcel lifted moist <!--
 * -->adiabatically to a new level.
 *
 * With a given parcel defined by a pressure and temperature (in Pascals and
 * degrees Kelvin), lift it moist adiabatically to a new pressure level 
 * (in Pascals) and return the temperture of the parcel at that level.
 *
 * This function relies on the Wobus Function ( sharp::wobf ), and it was shown
 * by Robert Davies-Jones (2007) that the Wobus function has a slight
 * dependence on pressure, which results in errors of up to 1.2 degrees Kelvin
 * in the temperature of a lifted parcel.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (degK)
 * \param    lifted_pressure       (Pa)
 *
 * \return   lifted_temperature    (degK)
 */
[[nodiscard]] float wetlift(float pressure, float temperature,
                            float lifted_pressure) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the temperature of a parcel lifted moist <!--
 * -->adiabatically to a new level.
 *
 *  With a given parcel defined by a pressure and temperature (in Pascals and
 *  degrees Kelvin), lift it moist adiabatically to a new pressure level 
 *  (in Pascals) and return the temperature of the parcel at that level.
 *
 *  This function is based on/ripped from George Byran's (NCAR) routine in CM1.
 *  It has several options for convergence criteria, pressure increments, and
 *  the type of moist adiabat/assumptions (i.e. pseudoadiabatic, adiabatic, 
 *  liquid, liquid+ice). The sharp::adiabat defines the types of adiabats, 
 *  and the float references to qv, ql, and qi are used to keep budgets of the
 *  moisture variables. 
 *
 *  NOTE: qv_total should most likely be the water vapor mixing ratio either
 *  at the parcel's sharp::LPL, or the water vapor mixing ratio at the LCL
 *  (these are going to be the same value). Essentially, the total water vapor
 *  before condensation. 
 *
 * \param   pressure        (Pa)
 * \param   temperature     (degK)
 * \param   new_pressure    (Pa)
 * \param   qv_total        (kg/kg)
 * \param   qv              (kg/kg)
 * \param   ql              (kg/kg)
 * \param   qi              (kg/kg)
 * \param   pres_incr       (Pa)
 * \param   converge        (precision)
 * \param   ma_type         (sharp::adiabat)
 *
 * \return  pcl_temperature (degK)
 */
[[nodiscard]] float moist_adiabat_cm1(float pressure, float temperature,
                                      float new_pressure, float& qv_total,
                                      float& qv, float& ql, float& qi,
                                      const float pres_incr,
                                      const float converge,
                                      const adiabat ma_type);

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Lift a parcel dry adiabatically to its Lifted Condensation Level
 * (LCL).
 *
 * Given a parcel's initial pressure (Pa), temperature (degK), and
 * dewpoint temperature (degK), lift the parcel dry adiabatically to its
 * Lifted Condensation Level and store the resulting LCL pressure (Pa)
 * and LCL temperature (degK) in the variables passed by reference to the
 * routine.
 *
 *
 * The LCL temperature is computed using an approximation. See the
 * sharp::lcl_temperature documentation for more information.
 *
 * \param    pressure              (Pa)
 * \param    temperature           (degK)
 * \param    dewpoint              (degK)
 * \param    pressure_at_lcl       (Pa)
 *
 * \param    temperature_at_lcl    (degK)
 */
void drylift(float pressure, float temperature, float dewpoint,
             float& pressure_at_lcl, float& temperature_at_lcl) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the temperature of an unsaturated parcel lifted to <!--
 * -->a given pressure level.
 *
 * This routine computes the temperature required to derive the Lifted Index
 * for a particulat pressure level. Given a parcel's initial pressure
 * (Pa), temperature (degK), and dewpoint (degK), it first lifts a
 * parcel to its LCL, and then continues to lift it moist adiabatically to the
 * given lifted pressure level (Pa).
 *
 * The LCL temperature is computed using an approximation. See the
 * sharp::lcl_temperature documentation for more information.
 * The moist adiabatic ascent is done by calling sharp::wetlift, which
 * relies on the Wobus Function ( sharp::wobf ). There are inherent
 * estimation errors, so see documentation to learn more.
 *
 * \param   pressure                    (Pa)
 * \param   temperature                 (degK)
 * \param   dewpoint                    (degK)
 * \param   lifted_pressure             (Pa)
 *
 * \return  lifted_index_temperature    (degK)
 */
[[nodiscard]] float lifted(float pressure, float temperature, float dewpoint,
                           float lifted_pressure) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the wetbulb temperature.
 *
 * Compute the wet bulb temperature (degK) given the pressure
 * (Pa), temperature (degK), and dewpoint (degK).
 *
 * First, it lifts a parcel with the given pressure, temperature, and
 * dewpoint temperature to its Lifted Condensation Level (LCL). To
 * compute the temperature and pressure of the LCL, an approximation
 * is used. See the sharp::lcl_temperature documentation for more
 * information.
 *
 * After the parcel has reached the LCL, the sharp::wetlift routine
 * lowers the parcel to its initial pressure level along a moist adiabat.
 * The sharp::wetlift routine relies on the Wobus Function ( sharp::wobf ),
 * which is an approximation with some inherent errors. See the
 * sharp::wetlift documentation for more information.
 *
 * \param   pressure                (Pa)
 * \param   temperature             (degK)
 * \param   dewpoint                (degK)
 *
 * \return  wetbulb_temperature     (degK)
 */
[[nodiscard]] float wetbulb(float pressure, float temperature,
                            float dewpoint) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the wetbulb potential temperature.
 *
 * Compute the wet bulb potential temperature (degK) given
 * the pressure (Pa), temperature (degK), and dewpoint
 * (degK).
 *
 * First, it lifts a parcel with the given pressure, temperature, and
 * dewpoint temperature to its Lifted Condensation Level (LCL). To
 * compute the temperature and pressure of the LCL, an approximation
 * is used. See the sharp::lcl_temperature documentation for more
 * information.
 *
 * After the parcel has reached the LCL, the sharp::wetlift routine
 * lowers the parcel to the standard reference pressure level
 * (1000.0 hPa) along a moist adiabat. The sharp::wetlift routine relies
 * on the Wobus Function ( sharp::wobf ), which is an approximation with
 * some inherent errors. See the sharp::wetlift documentation for more
 * information.
 *
 * \param   pressure                        (Pa)
 * \param   temperature                     (degK)
 * \param   dewpoint                        (degK)
 *
 * \return  wetbulb_potential_temperature   (degK)
 */
[[nodiscard]] float theta_wetbulb(float pressure, float temperature,
                                  float dewpoint) noexcept;

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Compute the equivalent potential temperature.
 *
 * Compute the equivalent potential temperature (degK) given
 * the pressure (Pa), temperature (degK), and dewpoint
 * (degK).
 *
 * First, it lifts a parcel with the given pressure, temperature, and
 * dewpoint temperature to its Lifted Condensation Level (LCL). To
 * compute the temperature and pressure of the LCL, an approximation
 * is used. See the sharp::lcl_temperature documentation for more
 * information.
 *
 * After the parcel has reached the LCL, the sharp::wetlift routine
 * lifts the parcel to 100 hPa along a moist adiabat. Finally, the
 * parcel is then lowered dry adiabatically to the standard reference
 * pressure level of 1000.0 hPa. The sharp::wetlift routine relies on
 * the Wobus Function ( sharp::wobf ), which is an approximation with
 * some inherent errors. See the sharp::wetlift documentation for
 * more information.
 *
 * \param    pressure                           (Pa)
 * \param    temperature                        (degK)
 * \param    dewpoint                           (degK)
 *
 * \return   equivalent_potential_temperature   (degK)
 */
[[nodiscard]] float thetae(float pressure, float temperature,
                           float dewpoint) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief compute the lapse rate over the given sharp::HeightLayer (AGL)
 *
 * Computes the lapse rate over a given sharp::HeightLayer above-ground-level.
 * This routine handles converting AGL to MSL by adding the surface height
 * value to the layer.
 *
 * \param   layer_agl               (meters AGL)
 * \param   height                  (meters MSL)
 * \param   temperature             (degK)
 * \param   N                       (length of arrays)
 *
 * \return  Temperature Lapse Rate  (degK/km)
 */
[[nodiscard]] float lapse_rate(HeightLayer layer_agl, const float height[],
                               const float temperature[], const int N) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief compute the lapse rate over the given sharp::PressureLayer 
 *
 * Computes the lapse rate over a given sharp::PressureLayer. This routine
 * handles converting the pressure layer into a height layer, and
 * then calles the sharp::HeightLayer implementation of this routine.
 *
 * \param   layer                   (meters AGL)
 * \param   pressure                (hPa)
 * \param   height                  (meters MSL)
 * \param   temperature             (degK)
 * \param   N		                (length of arrays)
 *
 * \return  Temperature Lapse Rate  (degK/km)
 */
[[nodiscard]] float lapse_rate(PressureLayer layer, const float pressure[],
                               const float height[], const float temperature[],
                               const int N) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief compute buoyancy given parcel and environment temperatures
 *
 * \param   pcl_temperature     (degK)
 * \param   env_temperature	    (degK)
 *
 * \return  buoyancy            (m/s^2)
 */
[[nodiscard]] float buoyancy(float pcl_temperature,
                             float env_temperature) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute moist static energy.
 *
 * \param   height_agl          (meters)
 * \param   temperature         (degK)
 * \param   specific_humidity   (kg/kg)
 *
 * \return  moist static energy ()
 */
[[nodiscard]] float moist_static_energy(float height_agl, float temperature,
                                        float specific_humidity) noexcept;

[[nodiscard]] float buoyancy_dilution_potential(const float temperature,
                                                const float mse_bar,
                                                const float saturation_mse);

}  // end namespace sharp

#endif // __SHARP_THERMP_H__
