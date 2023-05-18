#ifndef __SHARP_LAYER_WRAP_H__
#define __SHARP_LAYER_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif


struct sharp_PressureLayer {
    void *obj;
};

struct sharp_HeightLayer {
    void *obj;
};

struct sharp_LayerIndex {
    void *obj;
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
                                   const float* pressure, int NZ);

void sharp_get_HeightLayer_index(sharp_HeightLayer_t* lyr,
                                 sharp_LayerIndex_t* idx, 
                                 const float* height, int NZ);

void sharp_HeightLayer_to_PressureLayer(sharp_HeightLayer_t* hlyr,
                                        sharp_PressureLayer_t* plyr,
                                        const float* pressure,
                                        const float* height, int NZ,
                                        int isAGL);

void sharp_PressureLayer_to_HeightLayer(sharp_PressureLayer_t* plyr,
                                        sharp_HeightLayer_t* hlyr,
                                        const float* pressure,
                                        const float* height, int NZ,
                                        int toAGL);

float sharp_PressureLayer_min(sharp_PressureLayer_t* plyr, const float* pres,
                              const float* data, int NZ, float* lvl_of_min);

float sharp_PressureLayer_max(sharp_PressureLayer_t* plyr, const float* pres,
                              const float* data, int NZ, float* lvl_of_max);

float sharp_HeightLayer_min(sharp_HeightLayer_t* hlyr, const float* hght,
                            const float* data, int NZ, float* lvl_of_min);

float sharp_HeightLayer_max(sharp_HeightLayer_t* hlyr, const float* hght,
                            const float* data, int NZ, float* lvl_of_max);

float sharp_PressureLayer_mean(sharp_PressureLayer_t* plyr,
                               const float* pressure, const float* data,
                               int NZ);

float sharp_HeightLayer_mean(sharp_HeightLayer_t* hlyr, const float* height,
                             const float* pressure, const float* data, int NZ,
                             int isAGL);

float sharp_PressureLayer_integrate(sharp_PressureLayer_t* plyr,
                                    const float* data, const float* pressure,
                                    int NZ, int integ_sign, int weighted);

float sharp_HeightLayer_integrate(sharp_HeightLayer_t* hlyr,
                                  const float* data, const float* height,
                                  int NZ, int integ_sign, int weighted);

#ifdef __cplusplus
}
#endif

#endif
