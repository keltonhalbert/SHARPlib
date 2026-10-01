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
#ifndef SHARP_LAYERS_H
#define SHARP_LAYERS_H

#include <SHARPlib/algorithms.h>
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>

#include <cmath>
#include <cstddef>
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
     * \brief Construct an empty sharp::HeightLayer
     *
     * Sets the top and bottom to sharp::MISSING
     */
    HeightLayer() {
        bottom = MISSING;
        top = MISSING;
        delta = 100.0;
    }

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
     * \brief Construct an empty sharp::PressureLayer
     *
     * Sets the top and bottom to sharp::MISSING
     */
    PressureLayer() {
        bottom = MISSING;
        top = MISSING;
        delta = -1000.0;
    }

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
 * interpolation, with the option of returning the layer in meters
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
 * QC builds skip MISSING and NaN data and return MISSING for a layer
 * wholly outside the profile. sharp::layer_min gives the details.
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

#ifndef NO_QC
    // Clipping to the profile inverts a correctly ordered layer exactly when
    // the layer lies wholly outside the profile. Clipping never moves the
    // endpoint nearest the profile. For a layer above the profile, that is
    // the bottom, and clipping moves the top onto the last level. For a
    // layer below the profile, it is the top.
    const bool outside = (layer.coord == LayerCoordinate::pressure)
                             ? (layer.bottom < layer.top)
                             : (layer.bottom > layer.top);
    if (outside) {
        if (lvl_min_or_max) {
            const bool above = (layer.top == coord_arr[N - 1]);
            *lvl_min_or_max = above ? layer.bottom : layer.top;
        }
        return MISSING;
    }
#endif

    float min_or_max = MISSING;
    float top_val = MISSING;
    if constexpr (layer.coord == LayerCoordinate::pressure) {
        min_or_max = interp_pressure(layer.bottom, coord_arr, data_arr, N);
        top_val = interp_pressure(layer.top, coord_arr, data_arr, N);
    } else {
        min_or_max = interp_height(layer.bottom, coord_arr, data_arr, N);
        top_val = interp_height(layer.top, coord_arr, data_arr, N);
    }

    // QC builds skip MISSING and NaN values. A MISSING min_or_max means no
    // value yet, which happens when the bottom endpoint has no valid level
    // to interpolate from.
    const auto replaces = [&](const float val) {
#ifndef NO_QC
        if ((val == MISSING) || std::isnan(val)) return false;
        if (min_or_max == MISSING) return true;
#endif
        return comp(val, min_or_max);
    };

    float coord_lvl = layer.bottom;
    for (std::ptrdiff_t k = layer_idx.kbot; k < layer_idx.ktop + 1; ++k) {
        const float val = data_arr[k];
        if (replaces(val)) {
            min_or_max = val;
            coord_lvl = coord_arr[k];
        }
    }

    if (replaces(top_val)) {
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
 * QC builds, the default, skip levels whose data is MISSING or NaN. They
 * interpolate the layer bottom and top across missing levels, as
 * sharp::interp_height and sharp::interp_pressure do, and skip an endpoint
 * that has no valid level on one side of it. A layer with no valid data
 * returns MISSING. A layer that lies wholly outside the profile also
 * returns MISSING, and lvl_of_min is set to the layer's endpoint nearest
 * the profile. Builds with NO_QC skip these checks.
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
 * QC builds, the default, skip levels whose data is MISSING or NaN. They
 * interpolate the layer bottom and top across missing levels, as
 * sharp::interp_height and sharp::interp_pressure do, and skip an endpoint
 * that has no valid level on one side of it. A layer with no valid data
 * returns MISSING. A layer that lies wholly outside the profile also
 * returns MISSING, and lvl_of_max is set to the layer's endpoint nearest
 * the profile. Builds with NO_QC skip these checks.
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
template <typename L>
[[nodiscard]] float integrate_layer_trapz(L layer, const float var_array[],
                                          const float coord_array[],
                                          const std::ptrdiff_t N,
                                          const int integ_sign = 0,
                                          const bool weighted = false) {
    float integrated = 0.0;
    float weights = 0.0;

    const bool isign = std::signbit(integ_sign);

    LayerIndex idx = get_layer_index(layer, coord_array, N);

    float var_lyr_bottom;
    float var_lyr_top;

    auto include = [integ_sign, isign](float layer_avg) -> float {
        return (float)((!integ_sign) | (isign == std::signbit(layer_avg)));
    };

    if constexpr (layer.coord == LayerCoordinate::height) {
        var_lyr_bottom = interp_height(layer.bottom, coord_array, var_array, N);
        var_lyr_top = interp_height(layer.top, coord_array, var_array, N);
    } else {
        var_lyr_bottom =
            interp_pressure(layer.bottom, coord_array, var_array, N);
        var_lyr_top = interp_pressure(layer.top, coord_array, var_array, N);
    }

    // When the layer falls entirely between two native levels
    // with no interior points, handle as a single trapezoid.
    // NOTE: weights is accumulated by _integ_trapz through its
    // reference parameter whenever weighted=true. The same weights
    // variable must be passed to every _integ_trapz call in this
    // function for the final division to be correct.
    if (idx.ktop < idx.kbot) {
        if ((var_lyr_bottom == MISSING) || (var_lyr_top == MISSING)) {
            return MISSING;
        }
        float layer_avg = _integ_trapz(var_lyr_top, var_lyr_bottom, layer.top,
                                       layer.bottom, weights, weighted);

        integrated += include(layer_avg) * layer_avg;
    } else {
        // Start from the interpolated layer bottom and walk
        // through native levels to the interpolated top.
        // When a level is MISSING, the next valid level forms a
        // wider trapezoid from the last valid point, bridging
        // the gap via linear interpolation.
        float coord_prev = layer.bottom;
        float var_prev = var_lyr_bottom;
        bool have_prev = (var_lyr_bottom != MISSING);
        bool any_segment = false;

        for (std::ptrdiff_t k = idx.kbot; k <= idx.ktop; ++k) {
#ifndef NO_QC
            if (var_array[k] == MISSING) {
                continue;
            }
#endif
            if (have_prev) {
                float layer_avg =
                    _integ_trapz(var_array[k], var_prev, coord_array[k],
                                 coord_prev, weights, weighted);

                integrated += include(layer_avg) * layer_avg;
                any_segment = true;
            }

            coord_prev = coord_array[k];
            var_prev = var_array[k];
            have_prev = true;
        }

        // Final segment: last valid native level to the
        // interpolated layer top
        if (have_prev && (var_lyr_top != MISSING)) {
            float layer_avg = _integ_trapz(var_lyr_top, var_prev, layer.top,
                                           coord_prev, weights, weighted);

            integrated += include(layer_avg) * layer_avg;
            any_segment = true;
        }
        if (!any_segment) return MISSING;
    }

    if constexpr (layer.coord == LayerCoordinate::pressure) {
        integrated *= -1.0f;
        weights *= -1.0f;
    }

    if (weighted) {
        if (weights == 0.0f) return MISSING;
        integrated /= weights;
    }
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
 * A layer that extends past the profile is clipped to it. A layer that lies
 * wholly outside the profile, or touches it at only one level, has no mean
 * and returns sharp::MISSING.
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
 * A layer that extends past the profile is clipped to it. A layer that lies
 * wholly outside the profile, or touches it at only one level, has no mean
 * and returns sharp::MISSING.
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

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Walks a height profile once and reports each layer above or <!--
 * -->below a threshold.
 *
 * Walks the profile from the bottom up and calls on_layer once per layer,
 * in order of increasing height. Each call gets the layer bounds, the side
 * of the threshold (above = true where the data exceed it), and two areas
 * between the data and the threshold, in data units times meters:
 *
 * - pos_area: the integral of max(value - threshold, 0) dz over the layer
 * - neg_area: the integral of min(value - threshold, 0) dz (never positive)
 *
 * Layer boundaries are linearly interpolated threshold crossings. Areas
 * are trapezoids, split exactly at the crossings.
 *
 * Levels exactly at the threshold never start a new layer. They continue
 * the current one and add their depth to it, and a crossing is placed
 * where the data leave the threshold toward the other side. Levels at the
 * threshold below the first level off it take that level's side. If every
 * valid level is at the threshold, or there is only one valid level, the
 * walker reports one layer spanning the valid levels, with above = false
 * and zero areas. With no valid levels it reports nothing.
 *
 * min_depth and min_area merge noise into the surrounding layers. A run is
 * a raw stretch of data between crossings, or between a crossing and the
 * end of the profile. It is significant if its depth is at least min_depth
 * and the absolute value of its area is at least min_area. Only the run's
 * own depth and area count, never what merges into it. The merge rule:
 *
 * 1. Consecutive insignificant runs form a zone.
 * 2. A zone between significant runs on the same side is absorbed, and the
 *    layer continues through it.
 * 3. A zone between significant runs on opposite sides goes entirely to
 *    the side whose runs make up more of the zone's depth. A tie goes to
 *    the lower layer. Swapping above and below gives the mirror result,
 *    but reversing the profile vertically does not.
 * 4. A zone at the bottom or top of the profile joins the adjacent
 *    significant layer. If no run is significant, the whole column is one
 *    layer, on the side with more total run depth. A tie goes to the side
 *    of the lowest run.
 * 5. Merging moves boundaries only. A merged layer's pos_area and neg_area
 *    are the sums over all of its runs.
 * 6. With min_depth = min_area = 0, every run is significant, and the
 *    layers are exactly the raw threshold crossings.
 *
 * The walk is a single forward pass with constant state that reads each
 * level once. The walker itself never allocates or throws, though data_at
 * and on_layer may. A layer's top depends on the zone above it, so a
 * layer is reported once the next significant run is confirmed, which
 * happens as soon as that run's own depth and area reach the thresholds.
 * With both thresholds at 0, that is the first level past the layer's top
 * crossing. A layer at the top of the profile is reported at the end of
 * the walk.
 *
 * Unless NO_QC is defined, levels whose height or value is sharp::MISSING
 * or NaN are skipped, and the walk joins their valid neighbors with a
 * straight line.
 *
 * Heights must be strictly increasing. This is not checked.
 *
 * \tparam  Accessor    callable as float(std::ptrdiff_t k)
 * \tparam  Callback    callable as bool(const sharp::HeightLayer& layer,
 *                      bool above, float pos_area, float neg_area)
 *
 * \param   height      (meters)
 * \param   data_at     (returns the data value at level k)
 * \param   N           (length of the profile)
 * \param   threshold   (data units)
 * \param   min_depth   (meters; 0 disables)
 * \param   min_area    (data units times meters; 0 disables)
 * \param   on_layer    (called once per layer; returning false stops the
 *                      walk)
 *
 * \return  The number of layers reported, including the one whose
 *          callback stopped the walk
 */
template <typename Accessor, typename Callback>
std::ptrdiff_t for_each_threshold_layer(const float height[], Accessor data_at,
                                        const std::ptrdiff_t N,
                                        const float threshold,
                                        const float min_depth,
                                        const float min_area,
                                        Callback on_layer) {
    // Reads the next valid level as its height z and its departure d from
    // the threshold. Each level is read once, and the data only where the
    // height is valid.
    std::ptrdiff_t k = 0;
    auto next_level = [&](float& z, float& d) -> bool {
        for (; k < N; ++k) {
            const float z_k = height[k];
#ifndef NO_QC
            if ((z_k == MISSING) || std::isnan(z_k)) continue;
#endif
            const float val_k = data_at(k);
#ifndef NO_QC
            if ((val_k == MISSING) || std::isnan(val_k)) continue;
#endif
            ++k;
            z = z_k;
            d = val_k - threshold;
            return true;
        }
        return false;
    };

    std::ptrdiff_t reported = 0;
    auto report = [&](const float bottom, const float top, const bool above,
                      const float pos_area, const float neg_area) -> bool {
        HeightLayer layer;
        layer.bottom = bottom;
        layer.top = top;
        ++reported;
        return on_layer(layer, above, pos_area, neg_area);
    };

    float z_prev = 0.0f;
    float d_prev = 0.0f;
    if (!next_level(z_prev, d_prev)) return 0;
    const float z_first = z_prev;

    // The first run takes the side of the first level off the threshold,
    // and the levels at the threshold below it join that run. A lone valid
    // level, or a column entirely at the threshold, has no side.
    float z = 0.0f;
    float d = 0.0f;
    for (;;) {
        if (!next_level(z, d)) {
            report(z_first, z_prev, false, 0.0f, 0.0f);
            return reported;
        }
        if ((d_prev != 0.0f) || (d != 0.0f)) break;
        z_prev = z;
    }

    // The run in progress, and whether its own depth and area have reached
    // the thresholds yet.
    bool above = (d_prev != 0.0f) ? (d_prev > 0.0f) : (d > 0.0f);
    const bool first_above = above;
    float run_bot = z_first;
    float run_area = 0.0f;
    bool run_sig = false;

    // The open layer: the latest significant layer, whose top waits on the
    // zone above it. It exists once any run is significant.
    bool has_open = false;
    bool open_above = false;
    float open_bot = 0.0f;
    float open_pos = 0.0f;
    float open_neg = 0.0f;

    // The zone: finished insignificant runs above the open layer, or above
    // the column bottom before any run is significant.
    float zone_bot = z_first;
    float zone_depth_above = 0.0f;
    float zone_depth_below = 0.0f;
    float zone_pos = 0.0f;
    float zone_neg = 0.0f;

    auto significant = [&](const float depth, const float area) -> bool {
        return (depth >= min_depth) && (std::fabs(area) >= min_area);
    };

    // The run in progress just became significant. Give the zone below it
    // to a layer, and report the open layer if this run starts a new one.
    // Returns false if the callback stopped the walk.
    auto confirm_run = [&]() -> bool {
        run_sig = true;
        bool keep_going = true;
        if (!has_open) {
            // A zone at the bottom of the profile joins the first
            // significant layer.
            has_open = true;
            open_above = above;
            open_bot = z_first;
            open_pos = zone_pos;
            open_neg = zone_neg;
        } else if (open_above == above) {
            // A zone between significant runs on the same side is absorbed.
            open_pos += zone_pos;
            open_neg += zone_neg;
        } else {
            // A zone between opposite sides goes to the side with more of
            // its depth, and a tie goes to the lower layer. With an empty
            // zone, both branches put the boundary at run_bot == zone_bot.
            const float depth_open =
                (open_above) ? zone_depth_above : zone_depth_below;
            const float depth_run =
                (open_above) ? zone_depth_below : zone_depth_above;
            if (depth_open >= depth_run) {
                keep_going = report(open_bot, run_bot, open_above,
                                    open_pos + zone_pos, open_neg + zone_neg);
                open_bot = run_bot;
                open_pos = 0.0f;
                open_neg = 0.0f;
            } else {
                keep_going =
                    report(open_bot, zone_bot, open_above, open_pos, open_neg);
                open_bot = zone_bot;
                open_pos = zone_pos;
                open_neg = zone_neg;
            }
            open_above = above;
        }
        zone_depth_above = 0.0f;
        zone_depth_below = 0.0f;
        zone_pos = 0.0f;
        zone_neg = 0.0f;
        return keep_going;
    };

    // Ends the run in progress at top, crediting it to the open layer if it
    // is significant and to the zone if not. Returns false if the callback
    // stopped the walk.
    auto finish_run = [&](const float top) -> bool {
        const float depth = top - run_bot;
        if (!run_sig && significant(depth, run_area)) {
            if (!confirm_run()) return false;
        }
        if (run_sig) {
            if (above) {
                open_pos += run_area;
            } else {
                open_neg += run_area;
            }
            zone_bot = top;
        } else if (above) {
            zone_depth_above += depth;
            zone_pos += run_area;
        } else {
            zone_depth_below += depth;
            zone_neg += run_area;
        }
        run_sig = false;
        return true;
    };

    // _integ_trapz accumulates weights through a reference. Unused here.
    float weights = 0.0f;
    do {
        if ((above) ? (d < 0.0f) : (d > 0.0f)) {
            // A crossing between the two levels ends the run there.
            const float z_cross = lerp(z_prev, z, d_prev / (d_prev - d));
            run_area += _integ_trapz(0.0f, d_prev, z_cross, z_prev, weights);
            if (!finish_run(z_cross)) return reported;
            above = !above;
            run_bot = z_cross;
            run_area = _integ_trapz(d, 0.0f, z, z_cross, weights);
        } else {
            run_area += _integ_trapz(d, d_prev, z, z_prev, weights);
        }
        z_prev = z;
        d_prev = d;

        // Confirm the run as soon as it is significant, so the layer below
        // it is reported without waiting for this run to end.
        if (!run_sig && significant(z_prev - run_bot, run_area)) {
            if (!confirm_run()) return reported;
        }
    } while (next_level(z, d));

    if (!finish_run(z_prev)) return reported;
    if (has_open) {
        // A zone at the top of the profile joins the last significant layer.
        report(open_bot, z_prev, open_above, open_pos + zone_pos,
               open_neg + zone_neg);
    } else {
        // No run is significant: one layer on the side with more depth, and
        // a tie goes to the side of the lowest run.
        const bool column_above = (zone_depth_above == zone_depth_below)
                                      ? first_above
                                      : (zone_depth_above > zone_depth_below);
        report(z_first, z_prev, column_above, zone_pos, zone_neg);
    }
    return reported;
}

/// @cond DOXYGEN_IGNORE

extern template float layer_min<PressureLayer>(PressureLayer layer,
                                               const float coord_arr[],
                                               const float data_arr[],
                                               const std::ptrdiff_t N,
                                               float* lvl_of_min);

extern template float layer_min<HeightLayer>(HeightLayer layer,
                                             const float coord_arr[],
                                             const float data_arr[],
                                             const std::ptrdiff_t N,
                                             float* lvl_of_min);

extern template float layer_max<PressureLayer>(PressureLayer layer,
                                               const float coord_arr[],
                                               const float data_arr[],
                                               const std::ptrdiff_t N,
                                               float* lvl_of_max);

extern template float layer_max<HeightLayer>(HeightLayer layer,
                                             const float coord_arr[],
                                             const float data_arr[],
                                             const std::ptrdiff_t N,
                                             float* lvl_of_max);

extern template float integrate_layer_trapz<PressureLayer>(
    PressureLayer layer, const float var_array[], const float coord_array[],
    const std::ptrdiff_t N, const int integ_sign, const bool weighted);

extern template float integrate_layer_trapz<HeightLayer>(
    HeightLayer layer, const float var_array[], const float coord_array[],
    const std::ptrdiff_t N, const int integ_sign, const bool weighted);

/// @endcond

}  // end namespace sharp

#endif  // SHARP_LAYERS_H
