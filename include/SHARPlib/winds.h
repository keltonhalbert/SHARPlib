/**
 * \file
 * \brief Routines used to compute kinematic attributes and indices<!--
 * --> of vertical sounding profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef __SHARP_WINDS_H__
#define __SHARP_WINDS_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named floats that represent a wind vector.
 */
struct WindVector {
    /**
     * \brief The magntidue of a wind vector
     */
    float speed;

    /**
     * \brief The direction (degrees from North) of a wind vector.
     */
    float direction;
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two names floats that represent<!--
 * --> the components of a wind vector.
 */
struct WindComponents {
    /**
     * \brief The zonal (U) component of a wind vector.
     */
    float u;
    /**
     * \brief The meridional (V) component of a wind vector.
     */
    float v;
};

/**
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the zonal (U) wind component from a wind vector.
 *
 * Given the wind speed and direction, compute and return the zonal
 * (U) wind component.
 *
 * \param   wind_speed      (m/s)
 * \param   wind_direction  (degrees from North)
 *
 * \return  u_component     (m/s)
 */
[[nodiscard]] float u_component(float wind_speed, float wind_direction);

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the meridional (V) wind component from a wind vector.
 *
 * Given the wind speed and direction, compute and return the meridional
 * (V) wind component.
 *
 * \param   wind_speed      (m/s)
 * \param   wind_direction  (degrees from North)
 *
 * \return  v_component     (m/s)
 */
[[nodiscard]] float v_component(float wind_speed, float wind_direction);

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the direction (degrees from North) a vector points to.
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the direction (degrees from North) the vector points to.
 *
 * \param   u_comp          (m/s)
 * \param   v_comp          (m/s)
 *
 * \return  wind_direction  (degrees from North)
 */
[[nodiscard]] float vector_angle(float u_comp, float v_comp);

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Computes the magnitude of a vector given its components.
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the magnitude (m/s) of the vector.
 *
 * \param   u_comp      (m/s)
 * \param   v_comp      (m/s)
 *
 * \return  wind_speed  (m/s)
 */
[[nodiscard]] float vector_magnitude(float u_comp, float v_comp);

/**
 *
 * \author John Hart - NSSFC KCMO / NWSSPC OUN
 *
 * \brief Precisely computes the magnitude of a vector given its components.
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the magnitude (distance / time) of the vector. Instead of using
 * a standard sqrt(u*u + v*v) approach, this uses the c++ Standard Template
 * Library implementation of std::hypot. Using std::hypot is ~20 times slower,
 * but ensures a higher accuracy and checks for overflow/underflow when using
 * extremely large or extremely small values. Chances are, you are fine using
 * the imprecise version as long as values are in the range of
 * 1e-100 < values < 1e100.
 *
 * \param   u_comp      (m/s)
 * \param   v_comp      (m/s)
 *
 * \return  wind_speed  (m/s)
 */
[[nodiscard]] float vector_magnitude_precise(float u_comp, float v_comp);

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute wind speed and direction<!--
 * --> from components
 *
 * Given the zonal (U) and meridional (V) components of a vector, compute
 * and return the wind speed (m/s) and direction (degrees from
 * North) from the components as a struct.
 *
 * \param   u_comp  (m/s)
 * \param   v_comp  (m/s)
 *
 * \return  {wind_speed, wind_direction}
 */
[[nodiscard]] WindVector components_to_vector(float u_comp, float v_comp);

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute wind speed and direction<!--
 * --> from components
 *
 * Given the components of a vector via sharp::WindComponents, compute
 * and return the wind speed (m/s) and direction (degrees from
 * North) and return them as sharp::WindVector.
 *
 * \param   comp    {u_comp, v_comp}
 *
 * \return  {wind_speed, wind_direction}
 */
[[nodiscard]] WindVector components_to_vector(WindComponents comp);

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute the U and V components<!--
 * --> of a wind vector.
 *
 * Given the wind speed and direction of a vector, compute and return the
 * zonal and meridional vector components as a struct.
 *
 * \param   wind_speed      (m/s)
 * \param   wind_direction  (degrees from North)
 *
 * \return  {u_comp, v_comp}
 */
[[nodiscard]] WindComponents vector_to_components(float wind_speed,
                                                  float wind_direction);

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Conveniance function to compute the U and V components<!--
 * --> of a wind vector.
 *
 * Given the wind speed and direction of a vector via sharp::WindVector,
 * compute and return the zonal and meridional vector components as a struct.
 *
 * \param   vect    {speed, direction}
 *
 * \return  {u_comp, v_comp}
 */
[[nodiscard]] WindComponents vector_to_components(WindVector vect);

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compute the mean wind over the given sharp::PressureLayer
 *
 * Computes the mean wind over the given sharp::PressureLayer
 * and arrays of pressure and wind components with a length of N.
 *
 *
 * \param   layer       {bottom, top}
 * \param   pres        (Pa)
 * \param   u_wind      (m/s)
 * \param   v_wind      (m/s)
 * \param   N           (length of arrays)
 * \param   weighted    (whether to weight by pressure level)
 *
 * \return  {mean_u, mean_v}
 */
[[nodiscard]] WindComponents mean_wind(PressureLayer layer, const float pres[],
                                       const float u_wind[],
                                       const float v_wind[],
                                       const std::ptrdiff_t N,
                                       const bool weighted);

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Compue the wind shear for a given layer
 *
 * Computes the U and V components of the wind shear over a
 * layer given the vertical sounding arrays
 * of pressure/height, u_wind, v_wind, and their length.
 *
 * This generic template handles the logic
 * at compike time on whether or not this is a sharp::HeightLayer
 * or a sharp::PressureLayer, and calls the interpolation routines
 * accordingly
 *
 * \param   layer   {bottom, top}
 * \param   coord   (meters or Pa)
 * \param   u_wind  (m/s)
 * \param   v_wind  (m/s)
 * \param   N       (length of arrays)
 *
 * \return  {shear_u, shear_v}
 */
template <typename L>
[[nodiscard]] constexpr WindComponents wind_shear(L layer, const float coord[],
                                                  const float u_wind[],
                                                  const float v_wind[],
                                                  const std::ptrdiff_t N) {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING))
        return {MISSING, MISSING};
#endif

    float u_bot = MISSING;
    float u_top = MISSING;
    float v_bot = MISSING;
    float v_top = MISSING;
    // This if statement disappears at compile time depending
    // on which kind of layer is passed
    if constexpr (layer.coord == LayerCoordinate::pressure) {
        if (layer.bottom > coord[0]) layer.bottom = coord[0];
        if (layer.top < coord[N - 1]) layer.top = coord[N - 1];

        u_bot = interp_pressure(layer.bottom, coord, u_wind, N);
        u_top = interp_pressure(layer.top, coord, u_wind, N);

        v_bot = interp_pressure(layer.bottom, coord, v_wind, N);
        v_top = interp_pressure(layer.top, coord, v_wind, N);
    } else {
        // AGL to MSL
        layer.bottom += coord[0];
        layer.top += coord[0];

        if (layer.bottom < coord[0]) layer.bottom = coord[0];
        if (layer.top > coord[N - 1]) layer.top = coord[N - 1];

        u_bot = interp_height(layer.bottom, coord, u_wind, N);
        u_top = interp_height(layer.top, coord, u_wind, N);

        v_bot = interp_height(layer.bottom, coord, v_wind, N);
        v_top = interp_height(layer.top, coord, v_wind, N);
    }

#ifndef NO_QC
    if ((u_bot == MISSING) || (v_bot == MISSING) || (u_top == MISSING) ||
        (v_top == MISSING)) {
        return {MISSING, MISSING};
    }
#endif

    return {u_top - u_bot, v_top - v_bot};
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Computes the Storm Relative Helicity (SRH) over a given layer
 *
 * Computes the Storm Relative Helicity (SRH) over a given layer
 * using storm motion vector components stored in sharp::WindComponents.
 * This tempalte is the generalized form that will fill in the appropriate
 * interpolation calls and coordinate arrays depending on whether a
 * sharp::PressureLayer or sharp::HeightLayer gets passed to this
 * template. The other API instances are just static instantiations
 * of this template.
 *
 * This integration occurs over the given arrays of coord, u_wind, and v_wind,
 * with N elements in each. The integration only uses interpolation
 * for the top and bottom of the specified layer.
 *
 * If using a sharp::HeightLayer, units are expected in above-ground-level
 * (AGL), and gets converted to above-mean-sea-level (MSL) during the
 * computation.
 *
 * \param   layer           {bottom, top}
 * \param   storm_motion    {storm_u, storm_v}
 * \param   coord           (meters or Pa)
 * \param   u_wind          (m/s)
 * \param   v_wind          (m/s)
 * \param   N               (length of arrays)
 *
 * \return  storm_relative_helicity
 */
template <typename L>
[[nodiscard]] float helicity(L layer, WindComponents storm_motion,
                             const float coord[], const float u_wind[],
                             const float v_wind[], const std::ptrdiff_t N) {
#ifndef NO_QC
    if ((storm_motion.u == MISSING) || (storm_motion.v == MISSING)) {
        return MISSING;
    }
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    float sru_bot;
    float srv_bot;
    if constexpr (layer.coord == LayerCoordinate::height) {
        // get the height in MSL by adding the surface height
        layer.bottom += coord[0];
        layer.top += coord[0];

        // Get the interpolated bottom of the layer to integrate
        // and convert to storm-relative winds
        sru_bot = interp_height(layer.bottom, coord, u_wind, N);
        srv_bot = interp_height(layer.bottom, coord, v_wind, N);
    } else {
        // Get the interpolated bottom of the layer to integrate
        // and convert to storm-relative winds
        sru_bot = interp_pressure(layer.bottom, coord, u_wind, N);
        srv_bot = interp_pressure(layer.bottom, coord, v_wind, N);
    }
    sru_bot -= storm_motion.u;
    srv_bot -= storm_motion.v;

    // Get the vertical array indices corresponding to our layer,
    // while also bounds checking our search. Indices exclude the
    // top and bottom layers that will be interpolated.
    LayerIndex layer_idx = get_layer_index(layer, coord, N);

    // will get set in first loop iter
    float sru_top;
    float srv_top;
    float layer_helicity = 0.0;
    for (std::ptrdiff_t k = layer_idx.kbot; k <= layer_idx.ktop; ++k) {
#ifndef NO_QC
        if ((u_wind[k] == MISSING) || (v_wind[k] == MISSING)) {
            continue;
        }
#endif
        // top of layer storm relative winds
        sru_top = u_wind[k] - storm_motion.u;
        srv_top = v_wind[k] - storm_motion.v;

        // integrate layer
        layer_helicity += (sru_top * srv_bot) - (sru_bot * srv_top);
        // set the top to be the bottom
        // of the next layer
        sru_bot = sru_top;
        srv_bot = srv_top;
    }

    // Get the interpolated top of the layer to integrate
    // and convert to storm-relative winds
    if constexpr (layer.coord == LayerCoordinate::height) {
        sru_top = interp_height(layer.top, coord, u_wind, N);
        srv_top = interp_height(layer.top, coord, v_wind, N);
    } else {
        sru_top = interp_pressure(layer.top, coord, u_wind, N);
        srv_top = interp_pressure(layer.top, coord, v_wind, N);
    }
    sru_top -= storm_motion.u;
    srv_top -= storm_motion.v;

    // integrate the final layer
    layer_helicity += (sru_top * srv_bot) - (sru_bot * srv_top);

    return layer_helicity;
}

}  // end namespace sharp

#endif  // __SHARP_WINDS_H__
