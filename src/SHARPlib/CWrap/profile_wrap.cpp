#include <SHARPlib/CWrap/profile_wrap.h>
#include <SHARPlib/profile.h>
#include <stdlib.h>

sharp_Profile_t* sharp_Profile_create(int NZ, int sounding_type) {
    sharp_Profile_t* p;
    sharp::Profile* prof;
    sharp::Source src;
    src = static_cast<sharp::Source>(sounding_type);

    p = (sharp_Profile_t*)malloc(sizeof(*p));
    prof = new sharp::Profile(NZ, src);
    p->obj = prof;
    return p;
}

void sharp_Profile_delete(sharp_Profile_t* p) {
    if (p == NULL) return;
    delete static_cast<sharp::Profile*>(p->obj);
    free(p);
}

float* sharp_Profile_get_pres_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->pres;
}

float* sharp_Profile_get_hght_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->hght;
}

float* sharp_Profile_get_tmpk_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->tmpk;
}

float* sharp_Profile_get_dwpk_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->dwpk;
}

float* sharp_Profile_get_mixr_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->mixr;
}

float* sharp_Profile_get_vtmp_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->vtmp;
}

float* sharp_Profile_get_wspd_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->wspd;
}

float* sharp_Profile_get_wdir_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->wdir;
}

float* sharp_Profile_get_uwin_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->uwin;
}

float* sharp_Profile_get_vwin_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->vwin;
}

float* sharp_Profile_get_vvel_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->vvel;
}

float* sharp_Profile_get_theta_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->theta;
}

float* sharp_Profile_get_thetae_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->theta_e;
}

float* sharp_Profile_get_mse_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->moist_static_energy;
}

float* sharp_Profile_get_buoy_ptr(sharp_Profile_t* p) {
    if (p == NULL) return NULL;
    return static_cast<sharp::Profile*>(p->obj)->buoyancy;
}
