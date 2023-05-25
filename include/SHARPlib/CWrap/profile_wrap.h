/**
 * \file
 * \brief Header defining the C wrapper interface to the Profile structure 
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

#ifndef __SHARP_PROFILE_WRAP_H__
#define __SHARP_PROFILE_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

struct sharp_Profile {
    void *obj;
};
typedef struct sharp_Profile sharp_Profile_t;

sharp_Profile_t* sharp_Profile_create(int NZ, int sounding_type); 
void sharp_Profile_delete(sharp_Profile_t* p);

float* sharp_Profile_get_pres_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_hght_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_tmpc_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_dwpc_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_mixr_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_relh_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_vtmp_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_wspd_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_wdir_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_uwin_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_vwin_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_vvel_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_theta_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_thetae_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_mse_ptr(sharp_Profile_t* p);
float* sharp_Profile_get_buoy_ptr(sharp_Profile_t* p);

#ifdef __cplusplus
}
#endif

#endif // __SHARP_PROFILE_WRAP_H__
