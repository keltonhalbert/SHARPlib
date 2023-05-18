#include <SHARPlib/CWrap/parcel_wrap.h>
#include <SHARPlib/parcel.h>
#include <stdlib.h>

sharp_Parcel_t* sharp_Parcel_create() {
    sharp_Parcel_t* p;
    sharp::Parcel* pcl;

    p = (sharp_Parcel_t*)malloc(sizeof(*p));
    pcl = new sharp::Parcel();
    p->obj = pcl;
    return p;
}

void sharp_Parcel_delete(sharp_Parcel_t* p) {
    if (p == NULL) return;
    delete static_cast<sharp::Parcel*>(p->obj);
    free(p);
}

float sharp_Parcel_get_pres(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->pres;
}

float sharp_Parcel_get_tmpc(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->tmpc;
}

float sharp_Parcel_get_dwpc(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->dwpc;
}

float sharp_Parcel_get_lcl_pres(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->lcl_pressure;
}

float sharp_Parcel_get_lfc_pres(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->lfc_pressure;
}

float sharp_Parcel_get_el_pres(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->eql_pressure;
}

float sharp_Parcel_get_cape(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->cape;
}

float sharp_Parcel_get_cinh(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<sharp::Parcel*>(p->obj)->cinh;
}

int sharp_Parcel_get_lpl(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<int>(static_cast<sharp::Parcel*>(p->obj)->source);
}

