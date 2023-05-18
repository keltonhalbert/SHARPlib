#ifndef __SHARP_PARCEL_WRAP_H__
#define __SHARP_PARCEL_WRAP_H__

#ifdef __cplusplus
extern "C" {
#endif

struct sharp_Parcel {
    void *obj;
};

typedef struct sharp_Parcel sharp_Parcel_t;

sharp_Parcel_t* sharp_Parcel_create();
void sharp_Parcel_delete(sharp_Parcel_t* pcl);

float sharp_Parcel_get_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_tmpc(sharp_Parcel_t* pcl);
float sharp_Parcel_get_dwpc(sharp_Parcel_t* pcl);
float sharp_Parcel_get_lcl_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_lfc_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_el_pres(sharp_Parcel_t* pcl);
float sharp_Parcel_get_cape(sharp_Parcel_t* pcl);
float sharp_Parcel_get_cinh(sharp_Parcel_t* pcl);
int sharp_Parcel_get_lpl(sharp_Parcel_t* pcl);

#ifdef __cplusplus
}
#endif

#endif // __SHARP_PARCEL_WRAP_H__
