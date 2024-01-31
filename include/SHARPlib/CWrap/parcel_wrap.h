/**
 * \file
 * \brief Header defining the C wrapper interface to parcel routines
 * \author
 *   Kelton Halbert                     \n
 *   Email: kelton.halbert@noaa.gov     \n
 * \date    2023-05-19
 *
 * Written for the NWS Storm Prediction Center,
 * based on NSHARP routines written by John Hart
 * Rich Thompson.
 */

#ifndef __SHARP_PARCEL_WRAP_H__
#define __SHARP_PARCEL_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <SHARPlib/CWrap/profile_wrap.h>

struct sharp_Parcel {
    void* obj;
};

typedef struct sharp_Parcel sharp_Parcel_t;

sharp_Parcel_t* sharp_Parcel_create();
void sharp_Parcel_delete(sharp_Parcel_t* pcl);

float sharp_Parcel_get_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_tmpk(sharp_Parcel_t* pcl);
float sharp_Parcel_get_dwpk(sharp_Parcel_t* pcl);
float sharp_Parcel_get_lcl_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_lfc_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_el_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_cape(sharp_Parcel_t* pcl);
float sharp_Parcel_get_cinh(sharp_Parcel_t* pcl);
int sharp_Parcel_get_lpl(sharp_Parcel_t* pcl);

void sharp_define_parcel(const float* pressure, const float* temperature,
                         const float* dewpoint, const float* wv_mixratio,
                         const float* theta, const float* thetae, const int N,
                         sharp_Parcel_t* pcl, int source);

void sharp_define_custom_parcel(sharp_Parcel_t* pcl, float pres, float tmpk,
                                float dwpk);

void sharp_lift_parcel_wobf(const float* pressure,
                            const float* virtual_temperature, float* buoyancy,
                            const int N, sharp_Parcel_t* pcl);

void sharp_find_lfc_el(sharp_Parcel_t* pcl, const float* pres,
                       const float* hght, const float* buoy, const int N);

void sharp_cape_cinh(const float* pressure, const float* height,
                     const float* buoyancy, const int N, sharp_Parcel_t* pcl);

#ifdef __cplusplus
}
#endif

#endif  // __SHARP_PARCEL_WRAP_H__
