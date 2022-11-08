/**
 * \file
 * \brief Utilities and structures useful for grouping objects, conducting QC, and things that otherwise don't have a home.  
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
#ifndef __SHARP_UTILS
#define __SHARP_UTILS

namespace sharp {

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
    float zbot;

    /**
     * \brief The top of the height layer (meters)
     */
    float ztop;

    /**
     * \brief The height interval with which to iterate over the layer (meters) 
     */
    float dz = 100.0;

    HeightLayer(float bot, float top, float delta);
    HeightLayer(float bot, float top);
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
    float pbot;

    /**
     * \brief The top of the pressure layer (hPa)
     */
    float ptop;

    /**
     * \brief The pressure interval with which to iterate over the layer (hPa)
     */
    float dp = 1;
    PressureLayer(float bot, float top, float delta);
    PressureLayer(float bot, float top);
};


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * \brief A simple structure of two named integer indices for the top and bottom of a layer
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
 * \param layer     (sharp::PressureLayer)  {pbot, ptop}
 * \param pressure  (hPa)
 * \param num_levs  (length of array)
 * \return          sharp::LayerIndex       {kbot, ktop}
 *
 */
LayerIndex get_layer_index(PressureLayer& layer, 
                           const float* pressure, int num_levs); 


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
 * \param layer     (sharp::HeightLayer)    {zbot, ztop}
 * \param height    (meters)
 * \param num_levs  (length of array)
 * \return          sharp::LayerIndex       {kbot, ktop}
 *
 */
LayerIndex get_layer_index(HeightLayer& layer, 
                           const float* height, int num_levs);


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
 * \param layer         (sharp::PressureLayer)  {pbot, ptop}
 * \param pressure      (hPa)
 * \param data_arr      (data array to find max on)
 * \param num_levs      (length of arrays)
 * \param pres_of_max   (Pressure level of max val; hpa)
 * \return max_value
 *
 */
float max_value(PressureLayer layer,   const float* pressure,
                const float* data_arr, int num_levs,
                float* pres_of_max); 


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Returns the maximim value in the given height layer.  
 *
 * Returns the maximum value observed within the given data array
 * over the given sharp::HeightLayer. The function bounds checks
 * the sharp::HeightLayer by calling sharp::get_layer_index.
 *
 * If hght_of_max is not a nullptr, then the pointer will be 
 * dereferenced and filled with the height of the maximum
 * value. 
 *
 * \param layer         (sharp::HeightLayer)  {zbot, ztop}
 * \param height        (meters)
 * \param data_arr      (data array to find max on)
 * \param num_levs      (length of arrays)
 * \param hght_of_max   (Height level of max val; meters)
 * \return max_value
 *
 */
float max_value(HeightLayer layer,     const float* height,
                const float* data_arr, int num_levs,
                float* hght_of_max);


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
 * \param layer         (sharp::PressureLayer)  {pbot, ptop}
 * \param pressure      (hPa)
 * \param data_arr      (data array to find min on)
 * \param num_levs      (length of arrays)
 * \param pres_of_min   (Pressure level of min val; hPa)
 * \return min_value
 *
 */
float min_value(PressureLayer layer,   const float* pressure,
                const float* data_arr, int num_levs,
                float* pres_of_min);


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Returns the minimum value in the given height layer.  
 *
 * Returns the minimum value observed within the given data array
 * over the given sharp::HeightLayer. The function bounds checks
 * the sharp::HeightLayer by calling sharp::get_layer_indices.
 *
 * If hght_of_min is not a nullptr, then the pointer will be 
 * dereferenced and filled with the height of the minimum
 * value. 
 *
 * \param layer         (sharp::HeightLayer)  {zbot, ztop}
 * \param height        (meters)
 * \param data_arr      (data array to find min on)
 * \param num_levs      (length of arrays)
 * \param hght_of_min   (Height level of min val; meters)
 * \return min_value
 *
 */
float min_value(HeightLayer layer,     const float* height,
                const float* data_arr, int num_levs,
                float* hght_of_min); 

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \breif Computes the mean value of a field over a given pressure layer.
 *
 * Computes the mean value of a given array of data and corresponding 
 * pressure coordinates over the given sharp::PressureLayer.
 *
 * \param layer     (sharp::PressureLayer)  {pbot, ptop}
 * \param pressure  (vertical pressure array; hPa)
 * \param data_arr  (The data for which to compute a mean)
 * \param num_levs  (length of pressure and data arrays)
 * \return mean_value
 *
 */
float mean_value(PressureLayer layer,   const float* pressure,
                 const float* data_arr, int num_levs);


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \breif Computes the mean value of a field over a given height layer.
 *
 * Computes the mean value of a given array of data and corresponding 
 * height coordinates over the given sharp::HeightLayer.
 *
 * \param layer     (sharp::PressureLayer)  {pbot, ptop}
 * \param pressure  (vertical pressure array; hPa)
 * \param data_arr  (The data for which to compute a mean)
 * \param num_levs  (length of pressure and data arrays)
 * \return mean_value
 *
 */
float mean_value(HeightLayer layer,     const float* height,
                 const float* data_arr, int num_levs);


} // end namespace sharp


namespace sharp::exper {


} // end namespce sharp::exper


#endif
