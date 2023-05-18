#include <SHARPlib/CWrap/profile_wrap.h> 
#include <SHARPlib/CWrap/parcel_wrap.h>
#include <SHARPlib/CWrap/layer_wrap.h>
#include <SHARPlib/CWrap/winds_wrap.h>
#include <stdio.h>

int main(int argc, char** argv) {

    int NZ = 40;

    sharp_Profile_t* prof = sharp_Profile_create(NZ, 0);
    sharp_Parcel_t* pcl = sharp_Parcel_create();
    sharp_PressureLayer_t* plyr = sharp_PressureLayer_create(1000.0, 500.0);
    sharp_HeightLayer_t* hlyr = sharp_HeightLayer_create(0.0, 3000.0);
    sharp_LayerIndex_t* lyridx = sharp_LayerIndex_create();
    sharp_WindComponents_t* comp = sharp_WindComponents_create();
    sharp_WindVector_t* vec = sharp_WindVector_create();

    printf("pbot: %f hPa ptop: %f hPa\n", sharp_PressureLayer_get_bottom(plyr), sharp_PressureLayer_get_top(plyr));
    printf("hbot: %f meters htop: %f meters\n", sharp_HeightLayer_get_bottom(hlyr), sharp_HeightLayer_get_top(hlyr));

    float* pres = sharp_Profile_get_pres_ptr(prof);
    int k;
    float P0 = 1000.0;
    float dp = 25.0;
    for (k = 0; k < NZ; ++k) {
        pres[k] = P0 - k*dp;
    }

    float* pres2 = sharp_Profile_get_pres_ptr(prof);
    for (k = 0; k < NZ; ++k) {
        printf("pres[%d] = %f\n", k, pres2[k]);
    }

	printf("pcl pres: %f pcl tmpc: %f pcl dwpc: %f\n", 
			sharp_Parcel_get_pres(pcl),
			sharp_Parcel_get_tmpc(pcl),
			sharp_Parcel_get_dwpc(pcl)
		  );
	printf("pcl lcl pres: %f pcl lfc pres: %f pcl el pres: %f\n",
			sharp_Parcel_get_lcl_pres(pcl),
			sharp_Parcel_get_lfc_pres(pcl),
			sharp_Parcel_get_el_pres(pcl)
		  );
	printf("pcl CAPE: %f pcl CINH: %f\n",
			sharp_Parcel_get_cape(pcl),
			sharp_Parcel_get_cinh(pcl)
		  );

    sharp_Profile_delete(prof);
    sharp_Parcel_delete(pcl);
    sharp_PressureLayer_delete(plyr);
    sharp_HeightLayer_delete(hlyr);
    sharp_LayerIndex_delete(lyridx);
    sharp_WindComponents_delete(comp);
    sharp_WindVector_delete(vec);
    return 0;
}
