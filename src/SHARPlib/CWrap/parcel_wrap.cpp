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
