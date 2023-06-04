/**
 * \file
 * \brief Header defining the C wrapper interface to kinematic routines
 * \author
 *   Kelton Halbert                     \n
 *   Email: kelton.halbert@noaa.gov     \n
 *   License: Apache 2.0                \n
 * \date    2023-05-19
 *
 * Written for the NWS Storm Prediction Center,
 * based on NSHARP routines written by John Hart 
 * Rich Thompson.
 */

#ifndef __SHARP_WINDS_WRAP_H__
#define __SHARP_WINDS_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <SHARPlib/CWrap/layer_wrap.h>

struct sharp_WindVector {
    void *obj;
};

struct sharp_WindComponents {
    void *obj;
};

typedef struct sharp_WindVector sharp_WindVector_t;
typedef struct sharp_WindComponents sharp_WindComponents_t;

sharp_WindVector_t* sharp_WindVector_create(float speed, float direction);
sharp_WindComponents_t* sharp_WindComponents_create(float u, float v);

void sharp_WindVector_delete(sharp_WindVector_t* v);
void sharp_WindComponents_delete(sharp_WindComponents_t* c);

float sharp_WindVector_get_speed(sharp_WindVector_t* v);
float sharp_WindVector_get_direction(sharp_WindVector_t* v);

float sharp_WindComponents_get_u(sharp_WindComponents_t* c);
float sharp_WindComponents_get_v(sharp_WindComponents_t* c);

float sharp_u_component(float wind_speed, float wind_direction);
float sharp_v_component(float wind_speed, float wind_direction);
float sharp_vector_angle(float u_comp, float v_comp);
float sharp_vector_magnitude(float u_comp, float v_comp);

void sharp_components_to_vector(sharp_WindComponents_t* comp,
                                sharp_WindVector_t* vec);
void sharp_vector_to_components(sharp_WindVector_t* vec,
                                sharp_WindComponents_t* comp);

void sharp_mean_wind(sharp_PressureLayer_t* plyr, sharp_WindComponents_t* cmp,
                     const float* pres, const float* u_wind,
                     const float* v_wind, const int N, const int weighted);

void sharp_PressureLayer_wind_shear(sharp_PressureLayer_t* plyr,
                                    sharp_WindComponents_t* cmp,
                                    const float* pres, const float* u_wind,
                                    const float* v_wind, const int N);

void sharp_HeightLayer_wind_shear(sharp_HeightLayer_t* hlyr,
                                  sharp_WindComponents_t* cmp,
                                  const float* hght, const float* u_wind,
                                  const float* v_wind, const int N);

float sharp_HeightLayer_helicity(sharp_HeightLayer_t* hlyr,
                                 sharp_WindComponents_t* storm_motion,
                                 const float* height, const float* u_wind,
                                 const float* v_wind, const int N);

float sharp_PressureLayer_helicity(sharp_PressureLayer_t* plyr,
                                   sharp_WindComponents_t* storm_motion,
                                   const float* pressure, const float* u_wind,
                                   const float* v_wind, const int N);

#ifdef __cplusplus // end of extern C
}
#endif

#endif // __SHARP_WINDS_WRAP_H__
