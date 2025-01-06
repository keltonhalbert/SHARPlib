/**
 * \file
 * \brief Data structures and functions for performing operations over <!--
 * --> atmospheric layers.
 *
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2022-11-02
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */
#ifndef __SHARP_LAYERS_H__
#define __SHARP_LAYERS_H__

#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>

#include <cmath>
#include <functional>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Enum defining the coordinate system of a layer
 */
enum class LayerCoordinate {
    /**
     * \brief Height coordinate
     */
    height = 0,

    /**
     * \brief Pressure coordinate
     */
    pressure = 1,

    /**
     * \brief End value for range checking parameters
     */
    END,
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A simple structure of named floats that represent a height layer.
 */
struct HeightLayer {
    /**
     * \brief The bottom of the height layer (meters)
     */
    float bottom;

    /**
     * \brief The top of the height layer (meters)
     */
    float top;

    /**
     * \brief The height interval with which to iterate over the layer (meters)
     */
    float delta;

    /**
     * \brief The coordinate system of the layer
     */
    static constexpr LayerCoordinate coord = LayerCoordinate::height;

    /**
     * \brief Constructs a sharp::HeightLayer
     *
     * \param   bot     (bottom of layer, meters)
     * \param   top     (top of layer, meters)
     * \param   delta   (height increment, meters)
     *
     */
    HeightLayer(float bot, float top, float delta = 100.0);
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A simple structure of two named floats that represent a <!--
 * --> pressure layer.
 */
struct PressureLayer {
    /**
     * \brief The bottom of the pressure layer (Pa)
     */
    float bottom;

    /**
     * \brief The top of the pressure layer (Pa)
     */
    float top;

    /**
     * \brief The pressure interval with which to iterate over the layer (Pa)
     */
    float delta;

    /**
     * \brief The coordinate system of the layer
     */
    static constexpr LayerCoordinate coord = LayerCoordinate::pressure;

    /**
     * \brief Constructs a sharp::PressureLayer
     *
     * \param   bot     (bottom of layer, Pa)
     * \param   top     (top of layer, Pa)
     * \param   delta   (pressure increment, Pa)
     *
     */
    PressureLayer(float bot, float top, float delta = -1000.0);
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief A simple structure of two named integer indices for the top <!--
 * -->and bottom of a layer
 */
struct LayerIndex {
    /**
     * \brief The array index of the bottom of the layer
     */
    std::ptrdiff_t kbot;

    /**
     * \brief The array index of the top of the layer
     */
    std::ptrdiff_t ktop;
};

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Returns the array indices corresponding to the given layer.
 *
 * This template function encapsulates the algorithm that both bounds
 * checks layer operations and returns the array indices to the corresponding
 * coordinate array, excluding the top and bottom of the layer. Specifically,
 * the LayerIndex is [bottom, top] exclusive, and can be interpreted as the
 * interior range the layer bounds. This behaviour is due to the fact many
 * algorithms used will get the interpolated top and bottom values of the
 * layer, meaning for looping purposes we only want the inner range of
 * indices.
 *
 * NOTE: If layer.bottom or layer.top are out of bounds, this function will
 * truncate the layer to the coordinate range of data provided by coord[]
 * in an attempt to gracefully continue and produce a result.
 * This will modify the layer and is why it is passed as a reference.
 * If you do not wish to have this behavior, it is up to the user to ensure
 * they assign meaningful values to the layers and not request data out of
 * bounds.
 *
 * \param   layer           {bottom, top}
 * \param   coord           (height or pressure)
 * \param   N               (length of array)
 * \param   bottom_comp     (function comparing the layer bottom to coord array)
 * \param   top_comp        (function comparing the layer top to coord array)
 *
 * \return  {kbot, ktop}
 */
template <typename L, typename Cb, typename Ct>
[[nodiscard]] constexpr LayerIndex get_layer_index(L& layer,
                                                   const float coord[],
                                                   const std::ptrdiff_t N,
                                                   const Cb bottom_comp,
                                                   const Ct top_comp) {
    // bounds check out search!
    if (bottom_comp(layer.bottom, coord[0])) {
        layer.bottom = coord[0];
    }
    if (top_comp(layer.top, coord[N - 1])) {
        layer.top = coord[N - 1];
    }

    // whether pressure or height coordiantes, the bottom
    // comparitor passed to the function will determine
    // how to search
    std::ptrdiff_t lower_idx = lower_bound(coord, N, layer.bottom, bottom_comp);
    std::ptrdiff_t upper_idx = upper_bound(coord, N, layer.top, bottom_comp);

    // If the layer top is in between two levels, this check ensures
    // that our index is below the top for interpolation reasons
    upper_idx -= ((top_comp(coord[upper_idx], layer.top)) & (upper_idx > 0));

    return {lower_idx, upper_idx};
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Finds the array indices corresponding to the <!--
 * -->given sharp::PressureLayer.
 *
 * Returns the array indices corresponding to the given sharp::PressureLayer,
 * and performs bounds checking on the layer. As part of the bounds checking,
 * the sharp::PressureLayer is modified if the bottom or top of the layer
 * is out of bounds, which is why it gets passed as a reference.
 *
 * If the exact values of the top and bottom of the layer are present,
 * their indices are ignored. The default behavior is that the bottom
 * and top of a layer is computed by interpolation by default, since
 * it may or may not be present in the native data.
 *
 * \param   layer       {bottom, top}
 * \param   pressure    (Pa)
 * \param   N           (length of array)
 *
 * \return  {kbot, ktop}
 */
[[nodiscard]] LayerIndex get_layer_index(PressureLayer& layer,
                                         const float pressure[],
                                         const std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Finds the array indices corresponding to the <!--
 * -->given sharp::HeightLayer.
 *
 * Returns the array indices corresponding to the given sharp::HeightLayer,
 * and performs bounds checking on the layer. As part of the bounds checking,
 * the sharp::HeightLayer is modified if the bottom or top of the layer
 * is out of bounds, which is why it gets passed as a reference.
 *
 * If the exact values of the top and bottom of the layer are present,
 * their indices are ignored. The default behavior is that the bottom
 * and top of a layer is computed by interpolation by default, since
 * it may or may not be present in the native data.
 *
 *
 * \param   layer   {bottom, top}
 * \param   height  (meters)
 * \param   N       (length of array)
 *
 * \return  {kbot, ktop}
 */
[[nodiscard]] LayerIndex get_layer_index(HeightLayer& layer,
                                         const float height[],
                                         const std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Converts a sharp::HeightLayer to a sharp::PressureLayer
 *
 * Converts a sharp::HeightLayer to a sharp::PressureLayer via
 * interpolation, with a flag to signal whether the input layer is
 * in meters AGL or meters MSL. If for some strange reason you
 * provide a HeightLayer that is out of the bounds of height[],
 * then the bottom and top of the output layer will be set to
 * sharp::MISSING.
 *
 * \param   layer       (meters)
 * \param   pressure    (Pa)
 * \param   height      (meters)
 * \param   N           (Length of arrays)
 * \param   isAGL       Whether the input layer is AGL / MSL (default: false)
 *
 * \return  {bottom, top}
 */
[[nodiscard]] PressureLayer height_layer_to_pressure(HeightLayer layer,
                                                     const float pressure[],
                                                     const float height[],
                                                     const std::ptrdiff_t N,
                                                     const bool isAGL = false);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Converts a sharp::PressureLayer to a sharp::HeightLayer
 *
 * Converts a sharp::PressureLayer to a sharp::HeightLayer via
 * interpolation, with the obtion of returning the layer in meters
 * AGL or MSL. If for some strange reason you provide a PressureLayer
 * that is out of the bounds of pressure[], then the bottom and top
 * of the output layer will be set to sharp::MISSING.
 *
 * \param   layer       (Pa)
 * \param   pressure    (Pa)
 * \param   height      (meters)
 * \param   N           (Length of arrays)
 * \param   toAGL       Flag whether to return meters AGL or MSL (default:
 * false)
 *
 * \return  {bottom, top}
 */
[[nodiscard]] HeightLayer pressure_layer_to_height(PressureLayer layer,
                                                   const float pressure[],
                                                   const float height[],
                                                   const std::ptrdiff_t N,
                                                   const bool toAGL = false);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Template for min/max value over sharp::PressureLayer <!--
 * -->or sharp::HeightLayer
 *
 * This tempalte function contains the algorithm for searching for either
 * the minimum or maximum value over a given sharp::PressureLayer or
 * sharp::HeightLayer. This function is wrapped by the sharp::min_value and
 * sharp::max_value functions, which pass the appropriate comparitor to the
 * template.
 *
 * If lvl_min_or_max is not a nullptr, then the pointer will be
 * dereferenced and filled with the pressure or height of the maximum/minum
 * value.
 *
 * \param   layer           (sharp::PressureLayer or sharp::HeightLayer)
 * \param   coord_arr       (pressure or height)
 * \param   data_arr        (data array to find max on)
 * \param   N               (length of arrays)
 * \param   lvl_min_or_max  (level of min/max val)
 * \param   comp            Comparitor (i.e. std::less or std::greater)
 *
 * \return  minmax_value
 */
template <typename L, typename C>
[[nodiscard]] constexpr float layer_minmax(L layer, const float coord_arr[],
                                           const float data_arr[],
                                           const std::ptrdiff_t N,
                                           float* lvl_min_or_max,
                                           const C comp) {
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    LayerIndex layer_idx = get_layer_index(layer, coord_arr, N);

    float min_or_max = MISSING;
    float top_val = MISSING;
    if constexpr (layer.coord == LayerCoordinate::pressure) {
        min_or_max = interp_pressure(layer.bottom, coord_arr, data_arr, N);
        top_val = interp_pressure(layer.top, coord_arr, data_arr, N);
    } else {
        min_or_max = interp_height(layer.bottom, coord_arr, data_arr, N);
        top_val = interp_height(layer.top, coord_arr, data_arr, N);
    }

    float coord_lvl = layer.bottom;
    for (std::ptrdiff_t k = layer_idx.kbot; k < layer_idx.ktop + 1; ++k) {
        const float val = data_arr[k];
        if (comp(val, min_or_max)) {
            min_or_max = val;
            coord_lvl = coord_arr[k];
        }
    }

    if (comp(top_val, min_or_max)) {
        min_or_max = top_val;
        coord_lvl = layer.top;
    }

    if (lvl_min_or_max) *lvl_min_or_max = coord_lvl;
    return min_or_max;
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Returns the minimum value in the given layer.
 *
 * Returns the minimum value observed within the given data array
 * over the given sharp::PressureLayer or sharp::HeightLayer.
 * The function bounds checks the layer by calling
 * sharp::get_layer_index.
 *
 * If lvl_of_min is not a nullptr, then the pointer will be
 * dereferenced and filled with the coordinate of the minimum
 * value.
 *
 * \param   layer       (sharp::PressureLayer or sharp::HeightLayer)
 * \param   coord_arr   (coordinate units; Pa or meters)
 * \param   data_arr    (data array to find min on)
 * \param   N           (length of arrays)
 * \param   lvl_of_min  (level of min val)
 *
 * \return  layer_min
 *
 */
template <typename L>
constexpr float layer_min(L layer, const float coord_arr[],
                          const float data_arr[], const std::ptrdiff_t N,
                          float* lvl_of_min = nullptr) {
    constexpr auto comp = std::less<float>();
    return layer_minmax(layer, coord_arr, data_arr, N, lvl_of_min, comp);
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Returns the maximim value in the given layer.
 *
 * Returns the maximum value observed within the given data array
 * over the given sharp::PressureLayer or sharp::HeightLayer.
 * The function bounds checks the layer by calling
 * sharp::get_layer_index.
 *
 * If lvl_of_max is not a nullptr, then the pointer will be
 * dereferenced and filled with the coordinate of the maximum
 * value.
 *
 * \param   layer           (sharp::PressureLayer or sharp::HeightLayer)
 * \param   coord_arr       (coordinate units; Pa or meters)
 * \param   data_arr        (data array to find max on)
 * \param   N               (length of arrays)
 * \param   lvl_of_max      (level of max val)
 *
 * \return  layer_max
 */
template <typename L>
constexpr float layer_max(L layer, const float coord_arr[],
                          const float data_arr[], const std::ptrdiff_t N,
                          float* lvl_of_max = nullptr) {
    constexpr auto comp = std::greater<float>();
    return layer_minmax(layer, coord_arr, data_arr, N, lvl_of_max, comp);
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Returns a trapezoidal integration of the given layer.
 *
 * Returns a trapezoidal integration of the given data array over the
 * given sharp::PressureLayer or sharp::HeightLayer. There is an additional
 * default argument that determines whether this is a weighted average
 * or not. The sign of the integration may be passed as well, i.e.
 * integrating only positive or negative values, by passing a 1 or -1 to
 * integ_sign.
 *
 * \param   layer           (sharp::PressureLayer or sharp::HeightLayer)
 * \param   var_array       (data array to integrate)
 * \param   coord_array     (coordinate array to integrate over)
 * \param   N               (length of arrays)
 * \param   weighted        (bool; default=false)
 * \param   integ_sign	    (int; default=0)
 *
 * \return  integrated_value
 */
template <typename T, typename L>
[[nodiscard]] constexpr T integrate_layer_trapz(L layer, const T var_array[],
                                                const T coord_array[],
                                                const std::ptrdiff_t N,
                                                const int integ_sign = 0,
                                                const bool weighted = false) {
    T var_lyr_bottom;
    T coord_lyr_bottom = layer.bottom;

    T var_lyr_top;
    T coord_lyr_top = layer.top;

    T integrated = 0.0;
    T weights = 0.0;

    const bool isign = std::signbit(integ_sign);

    // using constexpr means that this if statement optimizes
    // away at compile time since the layer coordinate is known
    LayerIndex idx = get_layer_index(layer, coord_array, N);
    if constexpr (layer.coord == LayerCoordinate::height) {
        var_lyr_bottom = interp_height(layer.bottom, coord_array, var_array, N);
        var_lyr_top = interp_height(layer.top, coord_array, var_array, N);
    } else {
        var_lyr_bottom =
            interp_pressure(layer.bottom, coord_array, var_array, N);
        var_lyr_top = interp_pressure(layer.top, coord_array, var_array, N);
    }

    for (std::ptrdiff_t k = idx.kbot; k < idx.ktop; ++k) {
#ifndef NO_QC
        if (var_array[k] == MISSING) {
            continue;
        }
#endif

        T coord_bottom = coord_array[k];
        T var_bottom = var_array[k];

        T coord_top = coord_array[k + 1];
        T var_top = var_array[k + 1];

        T layer_avg = _integ_trapz(var_top, var_bottom, coord_top, coord_bottom,
                                   weights, weighted);

        T cond = ((!integ_sign) | (isign == std::signbit(layer_avg)));
        integrated += cond * layer_avg;
    }

    // interpolated bottom of layer
    T layer_avg =
        _integ_trapz(var_array[idx.kbot], var_lyr_bottom, coord_array[idx.kbot],
                     coord_lyr_bottom, weights, weighted);
    T cond = ((!integ_sign) | (isign == std::signbit(layer_avg)));
    integrated += cond * layer_avg;

    // interpolated top of layer
    layer_avg = _integ_trapz(var_lyr_top, var_array[idx.ktop], coord_lyr_top,
                             coord_array[idx.ktop], weights, weighted);
    cond = ((!integ_sign) | (isign == std::signbit(layer_avg)));
    integrated += cond * layer_avg;

    if constexpr (layer.coord == LayerCoordinate::pressure) {
        integrated *= -1.0f;
        weights *= -1.0f;
    }

    if (weighted) integrated /= weights;
    return integrated;
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes the mass-weighted mean value of a field over <!--
 * -->a given pressure layer.
 *
 * Computes the mass-weighted mean value of given arrays of data
 * and corresponding pressure coordinates over the given sharp::PressureLayer.
 *
 * \param   layer       (sharp::PressureLayer)
 * \param   pressure    (vertical pressure array; Pa)
 * \param   data_arr    (The data for which to compute a mean)
 * \param   N           (length of pressure and data arrays)
 *
 * \return  layer_mean
 */
[[nodiscard]] float layer_mean(PressureLayer layer, const float pressure[],
                               const float data_arr[], const std::ptrdiff_t N);

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Computes the mass-weighted mean value of a field over <!--
 * -->a given height layer.
 *
 * Computes the mass-weighted mean value of given arrays of data
 * and corresponding height coordinates over the given sharp::HeightLayer.
 * This is really just a fancy wrapper around the implementation that uses
 * sharp::PressureLayer.
 *
 * \param   layer       (sharp::HeightLayer)
 * \param   height      (vertical height array; meters)
 * \param   pressure    (vertical pressure array; Pa)
 * \param   data_arr    (The data for which to compute a mean)
 * \param   N           (length of pressure and data arrays)
 * \param 	isAGL 		(whether or not intput is in AGL or MSL)
 *
 * \return  layer_mean
 */
[[nodiscard]] float layer_mean(HeightLayer layer, const float height[],
                               const float pressure[], const float data_arr[],
                               const std::ptrdiff_t N,
                               const bool isAGL = false);

}  // end namespace sharp

#endif  // __SHARP_LAYERS_H__
