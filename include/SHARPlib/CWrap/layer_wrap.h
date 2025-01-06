/**
 * \file
 * \brief Header defining the C wrapper interface to layer based routines
 * \author
 *   Kelton Halbert                     \n
 *   Email: kelton.halbert@noaa.gov     \n
 * \date    2023-05-19
 *
 * Written for the NWS Storm Prediction Center,
 * based on NSHARP routines written by John Hart
 * Rich Thompson.
 */

#ifndef SHARP_LAYER_WRAP_H
#define SHARP_LAYER_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

struct sharp_PressureLayer {
    void* obj;
};

struct sharp_HeightLayer {
    void* obj;
};

struct sharp_LayerIndex {
    void* obj;
};

typedef struct sharp_PressureLayer sharp_PressureLayer_t;
typedef struct sharp_HeightLayer sharp_HeightLayer_t;
typedef struct sharp_LayerIndex sharp_LayerIndex_t;

sharp_PressureLayer_t* sharp_PressureLayer_create(float bottom, float top);
sharp_HeightLayer_t* sharp_HeightLayer_create(float bottom, float top);
sharp_LayerIndex_t* sharp_LayerIndex_create();

void sharp_PressureLayer_delete(sharp_PressureLayer_t* l);
void sharp_HeightLayer_delete(sharp_HeightLayer_t* l);
void sharp_LayerIndex_delete(sharp_LayerIndex_t* i);

float sharp_PressureLayer_get_bottom(sharp_PressureLayer_t* l);
float sharp_PressureLayer_get_top(sharp_PressureLayer_t* l);

float sharp_HeightLayer_get_bottom(sharp_HeightLayer_t* l);
float sharp_HeightLayer_get_top(sharp_HeightLayer_t* l);

int sharp_LayerIndex_get_bottom(sharp_LayerIndex_t* i);
int sharp_LayerIndex_get_top(sharp_LayerIndex_t* i);

void sharp_get_PressureLayer_index(sharp_PressureLayer_t* lyr,
                                   sharp_LayerIndex_t* idx,
                                   const float* pressure, const int N);

void sharp_get_HeightLayer_index(sharp_HeightLayer_t* lyr,
                                 sharp_LayerIndex_t* idx, const float* height,
                                 const int N);

void sharp_HeightLayer_to_PressureLayer(sharp_HeightLayer_t* hlyr,
                                        sharp_PressureLayer_t* plyr,
                                        const float* pressure,
                                        const float* height, const int N,
                                        int isAGL);

void sharp_PressureLayer_to_HeightLayer(sharp_PressureLayer_t* plyr,
                                        sharp_HeightLayer_t* hlyr,
                                        const float* pressure,
                                        const float* height, const int N,
                                        int toAGL);

float sharp_PressureLayer_min(sharp_PressureLayer_t* plyr, const float* pres,
                              const float* data, const int N,
                              float* lvl_of_min);

float sharp_PressureLayer_max(sharp_PressureLayer_t* plyr, const float* pres,
                              const float* data, const int N,
                              float* lvl_of_max);

float sharp_HeightLayer_min(sharp_HeightLayer_t* hlyr, const float* hght,
                            const float* data, const int N, float* lvl_of_min);

float sharp_HeightLayer_max(sharp_HeightLayer_t* hlyr, const float* hght,
                            const float* data, const int N, float* lvl_of_max);

float sharp_PressureLayer_mean(sharp_PressureLayer_t* plyr,
                               const float* pressure, const float* data,
                               const int N);

float sharp_HeightLayer_mean(sharp_HeightLayer_t* hlyr, const float* height,
                             const float* pressure, const float* data,
                             const int N, int isAGL);

float sharp_PressureLayer_integrate(sharp_PressureLayer_t* plyr,
                                    const float* data, const float* pressure,
                                    const int N, const int integ_sign,
                                    const int weighted);

float sharp_HeightLayer_integrate(sharp_HeightLayer_t* hlyr, const float* data,
                                  const float* height, const int N,
                                  const int integ_sign, const int weighted);

#ifdef cplusplus
}
#endif

#endif
