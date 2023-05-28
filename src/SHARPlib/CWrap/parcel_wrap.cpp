#include <SHARPlib/CWrap/parcel_wrap.h>
#include <SHARPlib/constants.h>
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
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->pres;
}

float sharp_Parcel_get_tmpk(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->tmpk;
}

float sharp_Parcel_get_dwpk(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->dwpk;
}

float sharp_Parcel_get_lcl_pres(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->lcl_pressure;
}

float sharp_Parcel_get_lfc_pres(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->lfc_pressure;
}

float sharp_Parcel_get_el_pres(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->eql_pressure;
}

float sharp_Parcel_get_cape(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->cape;
}

float sharp_Parcel_get_cinh(sharp_Parcel_t* p) {
	if (p == NULL) return sharp::MISSING;
	return static_cast<sharp::Parcel*>(p->obj)->cinh;
}

int sharp_Parcel_get_lpl(sharp_Parcel_t* p) {
	if (p == NULL) return 0;
	return static_cast<int>(static_cast<sharp::Parcel*>(p->obj)->source);
}

void sharp_define_parcel(sharp_Profile_t* prof, sharp_Parcel_t* pcl,
                         int source) {
    if ((prof == NULL) || (pcl == NULL)) return;
    sharp::Profile* pf = static_cast<sharp::Profile*>(prof->obj);
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    sharp::LPL src = static_cast<sharp::LPL>(source);
    sharp::define_parcel(pf, pc, src); 
}

void sharp_define_custom_parcel(sharp_Parcel_t* pcl, float pres, float tmpk,
                                float dwpk) {
    if (pcl == NULL) return;
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    pc->pres = pres;
    pc->tmpk = tmpk;
    pc->dwpk = dwpk;
}

void sharp_lift_parcel_wobf(sharp_Profile_t* prof, sharp_Parcel_t* pcl) {
    if ((prof == NULL) || (pcl == NULL)) return;
    sharp::Profile* pf = static_cast<sharp::Profile*>(prof->obj);
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    constexpr sharp::lifter_wobus lifter;
    sharp::lift_parcel(lifter, pf, pc);
}

void sharp_find_lfc_el(sharp_Parcel_t* pcl, const float* pres,
                       const float* hght, const float* buoy, const int NZ) {
    if (pcl == NULL) return;
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    sharp::find_lfc_el(pc, pres, hght, buoy, NZ);
}

void sharp_cape_cinh(sharp_Profile_t* prof, sharp_Parcel_t* pcl) {
    if ((prof == NULL) || (pcl == NULL)) return;
    sharp::Profile* pf = static_cast<sharp::Profile*>(prof->obj);
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    sharp::cape_cinh(pf, pc);
}

void sharp_parcel_wobf(sharp_Profile_t* prof, sharp_Parcel_t* pcl) {
    if ((prof == NULL) || (pcl == NULL)) return;
    sharp::Profile* pf = static_cast<sharp::Profile*>(prof->obj);
    sharp::Parcel* pc = static_cast<sharp::Parcel*>(pcl->obj);
    sharp::parcel_wobf(pf, pc);
}

