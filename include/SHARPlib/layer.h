/**
 * \file
 * \brief Data structures and functions for performing operations over atmospheric layers.  
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-11-02
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */
#ifndef __SHARP_LAYERS
#define __SHARP_LAYERS

#include <functional>

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>

namespace sharp {

/*
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Enum defining the coordinate system of a layer
 *
 */
enum class LayerCoordinate {
	height = 0,
	pressure = 1,
};

/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named floats that represent a height layer. 
 *
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

    HeightLayer(float bot, float top, float delta=100.0);
};



/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named floats that represent a pressure layer. 
 *
 */
struct PressureLayer {

    /**
     * \brief The bottom of the pressure layer (hPa)
     */
    float bottom;

    /**
     * \brief The top of the pressure layer (hPa)
     */
    float top;

    /**
     * \brief The pressure interval with which to iterate over the layer (hPa)
     */
    float delta;

	/**
	 * \brief The coordinate system of the layer
	 */
	static constexpr LayerCoordinate coord = LayerCoordinate::pressure;

    PressureLayer(float bot, float top, float delta=-10);
};


/**
 *
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A simple structure of two named integer indices for the top and bottom of a layer
 *
 */
struct LayerIndex {
    /**
     * \brief The array index of the bottom of the layer
     */
    int kbot;

    /**
     * \brief The array index of the top of the layer
     */
    int ktop;
};


/*
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Returns the array indices corresponding to the given layer.
 *
 * This template function encapsulates the algorithm that both bounds
 * checks layer operations and returns the array indices to the corresponding
 * coordinate array, excluding the top and bottom of the layer. This behaviour
 * is due to the fact many algorithms used will get the interpolated top and 
 * bottom values of the layer, meaning for looping purposes we only want the 
 * inner range of indices. 
 *
 * \param layer	 (sharp::HeightLayer or sharp::PressureLayer) {bottom, top}
 * \param coord  (height or pressure)
 * \param N      (length of array)
 * \param bottom_comp (function comparing the layer bottom to the coordinate array)
 * \param top_comp    (function comparing the layer top to the coordinate array)
 *
 * \return       sharp::LayerIndex {kbot, ktop}
 */
template <typename L, typename Cb, typename Ct>
constexpr LayerIndex get_layer_index(L& layer, const float* coord, int N,
									 const Cb bottom_comp, const Ct top_comp) noexcept {

	if (bottom_comp(layer.bottom, coord[0])) {
		layer.bottom = coord[0];
	}

	if (top_comp(layer.top, coord[N-1])) {
		layer.top = coord[N-1];
	}

	// whether pressure or height coordiantes, the bottom
	// comparitor passed to the function will determine 
	// how to search
	int lower_idx = lower_bound(coord, N, layer.bottom, bottom_comp);
	int upper_idx = upper_bound(coord, N, layer.top, bottom_comp);

	if ((bottom_comp(coord[lower_idx], layer.bottom)) && (lower_idx < N-1)) {
		++lower_idx;
	}

	if ((top_comp(coord[upper_idx], layer.top)) && (upper_idx > 0)) {
		--upper_idx;
	}

	return {lower_idx, upper_idx};
}


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Finds the array indices corresponding to the given layer. 
 *
 * Returns the array indices corresponding to the given sharp::PressureLayer,
 * and performs bounds checking on the layer. As part of the bounds checking, 
 * the sharp::PressureLayer is modified if the bottom or top of the layer
 * is out of bounds, which is why it gets passed as a reference.
 *
 * If the exact values of the top and bottom of the layer are present,
 * their indices are ignored. The default behavior is that the bottom
 * and top of a layer is computed by interpolation by default, since
 * it may or may not be present. 
 *
 * \param layer     (sharp::PressureLayer)  {bottom, top}
 * \param pressure  (hPa)
 * \param num_levs  (length of array)
 * \return          sharp::LayerIndex       {kbot, ktop}
 *
 */
LayerIndex get_layer_index(PressureLayer& layer, 
          			       const float* pressure, 
                    	   int num_levs) noexcept; 


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Finds the array indices corresponding to the given layer. 
 *
 * Returns the array indices corresponding to the given sharp::HeightLayer,
 * and performs bounds checking on the layer. As part of the bounds checking, 
 * the sharp::HeightLayer is modified if the bottom or top of the layer
 * is out of bounds, which is why it gets passed as a reference.
 *
 * If the exact values of the top and bottom of the layer are present,
 * their indices are ignored. The default behavior is that the bottom
 * and top of a layer is computed by interpolation by default, since
 * it may or may not be present. 
 *
 *
 * \param layer     (sharp::HeightLayer)    {bottom, top}
 * \param height    (meters)
 * \param num_levs  (length of array)
 * \return          sharp::LayerIndex       {kbot, ktop}
 *
 */
LayerIndex get_layer_index(HeightLayer& layer, 
                           const float* height, 
                           int num_levs) noexcept;


/*
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Converts a sharp::HeightLayer to a sharp::PressureLayer
 *
 * Converts a sharp::HeightLayer to a sharp::PressureLayer via 
 * interpolation, with a flag to signal whether the input layer is 
 * in meters AGL or meters MSL. 
 *
 * \param layer     sharp::HeightLayer to convert (meters)
 * \param pressure  Vertical array of pressure (hPa)
 * \param height    Vertical array of height (meters)
 * \param num_levs  Length of arrays
 * \param isAGL     Flag whether the input layer is AGL or MSL (default: false)
 *
 * \return  sharp::PressureLayer
 */
PressureLayer height_layer_to_pressure(HeightLayer layer, 
                const float* pressure, const float* height, 
                int num_levs, bool isAGL=false) noexcept;


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Converts a sharp::PressureLayer to a sharp::HeightLayer
 *
 * Converts a sharp::PressureLayer to a sharp::HeightLayer via
 * interpolation, with the obtion of returning the layer in meters
 * AGL or MSL. 
 *
 * \param layer     sharp::PressureLayer to convert (hPa)
 * \param pressure  Vertical array of pressure (hPa)
 * \param height    Vertical array of height (meters)
 * \param num_levs  Length of arrays
 * \param toAGL     Flag whether to return meters AGL or MSL (default: false)
 *
 * \return sharp::HeightLayer
 */
HeightLayer pressure_layer_to_height(PressureLayer layer, 
                const float* pressure, const float* height, 
                int num_levs, bool toAGL=false) noexcept;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Template for min/max value over sharp::PressureLayer or sharp::HeightLayer 
 *
 * This tempalte function contains the algorithm for searching for either
 * the minimum or maximum value over a given sharp::PressureLayer or sharp::HeightLayer. 
 * This function is wrapped by the sharp::min_value and sharp::max_value
 * functions, which pass the appropriate comparitor to the template. 
 *
 * If lvl_min_or_max is not a nullptr, then the pointer will be 
 * dereferenced and filled with the pressure or height of the maximum/minum
 * value. 
 *
 * \param layer          (sharp::PressureLayer or sharp::HeightLayer)  {bottom, top}
 * \param coord_arr      (pressure or height)
 * \param data_arr       (data array to find max on)
 * \param N              (length of arrays)
 * \param lvl_min_or_max (level of min/max val; hpa)
 * \param comp           Comparitor (i.e. std::less or std::greater)
 * \return minmax_value
 *
 */
template <typename L, typename C>
constexpr float layer_minmax(L layer, const float* coord_arr, 
                   const float* data_arr, int N, float* lvl_min_or_max,
                   const C comp) noexcept {
    
#ifndef NO_QC
    if ((layer.bottom == MISSING) || (layer.top == MISSING)) {
        return MISSING;
    }
#endif

    LayerIndex layer_idx = get_layer_index(layer, coord_arr, N);

	float min_or_max = MISSING;
	if (layer.coord == LayerCoordinate::pressure) {
    	min_or_max = interp_pressure(layer.bottom, coord_arr, data_arr, N);
	}
	else {
		min_or_max = interp_height(layer.bottom, coord_arr, data_arr, N);
	}

    if (lvl_min_or_max)
        *lvl_min_or_max = layer.bottom;

    for (int k = layer_idx.kbot; k <= layer_idx.ktop; ++k) {
        float val = data_arr[k];    
        if (comp(val, min_or_max)) {
            min_or_max = val;
            if (lvl_min_or_max)
                *lvl_min_or_max = coord_arr[k];
        } 
    }

	float top_val = MISSING;
	if (layer.coord == LayerCoordinate::pressure) {
    	top_val = interp_pressure(layer.top, coord_arr, data_arr, N);
	}
	else {
		top_val = interp_height(layer.top, coord_arr, data_arr, N);
	}

    if (comp(top_val, min_or_max)) {
        min_or_max = top_val;
        if (lvl_min_or_max)
            *lvl_min_or_max = layer.top;
    }

    return min_or_max;
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Returns the minimum value in the given pressure layer. 
 *
 * Returns the minimum value observed within the given data array
 * over the given sharp::PressureLayer. The function bounds checks
 * the sharp::PressureLayer by calling sharp::get_layer_index.
 *
 * If pres_of_min is not a nullptr, then the pointer will be 
 * dereferenced and filled with the pressure of the minimum
 * value. 
 *
 * \param layer         (sharp::PressureLayer)  {bottom, top}
 * \param pressure      (hPa)
 * \param data_arr      (data array to find min on)
 * \param N             (length of arrays)
 * \param pres_of_min   (Pressure level of min val; hPa)
 * \return layer_min
 *
 */
template <typename L>
constexpr float layer_min(L layer, const float* coord_arr,
				const float* data_arr, int N,
				float* lvl_of_min=nullptr) noexcept {
	constexpr auto comp = std::less<float>();
	return layer_minmax(layer, coord_arr, data_arr, N, lvl_of_min, comp);
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Returns the maximim value in the given pressure layer. 
 *
 * Returns the maximum value observed within the given data array
 * over the given sharp::PressureLayer. The function bounds checks
 * the sharp::PressureLayer by calling sharp::get_layer_index.
 *
 * If pres_of_max is not a nullptr, then the pointer will be 
 * dereferenced and filled with the pressure of the maximum
 * value. 
 *
 * \param layer         (sharp::PressureLayer)  {bottom, top}
 * \param pressure      (hPa)
 * \param data_arr      (data array to find max on)
 * \param N             (length of arrays)
 * \param pres_of_max   (Pressure level of max val; hpa)
 * \return layer_max
 *
 */
template <typename L>
constexpr float layer_max(L layer, const float* coord_arr,
				const float* data_arr, int N,
				float* lvl_of_max=nullptr) noexcept {
	constexpr auto comp = std::greater<float>();
	return layer_minmax(layer, coord_arr, data_arr, N, lvl_of_max, comp);
}	


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Computes the mass-weighted mean value of a field over a given pressure layer.
 *
 * Computes the mean value of a given array of data and corresponding 
 * pressure coordinates over the given sharp::PressureLayer.
 *
 * \param layer     (sharp::PressureLayer)  {bottom, top}
 * \param pressure  (vertical pressure array; hPa)
 * \param data_arr  (The data for which to compute a mean)
 * \param num_levs  (length of pressure and data arrays)
 * \return layer_mean 
 *
 */
float layer_mean(PressureLayer layer,   const float* pressure,
                 const float* data_arr, int num_levs) noexcept;


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Computes the mass-weighted mean value of a field over a given height layer.
 *
 * Computes the mean value of a given array of data and corresponding 
 * height coordinates over the given sharp::HeightLayer. This is really
 * just a fancy wrapper around the implementation that uses 
 * sharp::PressureLayer. 
 *
 * \param layer     (sharp::PressureLayer)  {bottom, top}
 * \param height    (vertical height array; meters)
 * \param pressure  (vertical pressure array; hPa)
 * \param data_arr  (The data for which to compute a mean)
 * \param num_levs  (length of pressure and data arrays)
 * \return layer_mean 
 *
 */
float layer_mean(HeightLayer layer, const float* height, const float* pressure,
                 const float* data_arr, int num_levs) noexcept;


} // end namespace sharp


namespace sharp::exper {


} // end namespce sharp::exper


#endif
