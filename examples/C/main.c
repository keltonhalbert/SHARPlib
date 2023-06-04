#include <SHARPlib/CWrap/interp_wrap.h>
#include <SHARPlib/CWrap/layer_wrap.h>
#include <SHARPlib/CWrap/parcel_wrap.h>
#include <SHARPlib/CWrap/thermo_wrap.h>
#include <SHARPlib/CWrap/winds_wrap.h>
#include <SHARPlib/CWrap/params/convective_wrap.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

/**
 * Reads the sounding file and counts the number of
 * lines that correspond to data levels and returns
 * that number. This is so that arrays can be sized
 * properly for the data available.
 */
int get_num_levs_in_file(const char* filename) {
    FILE* file = fopen(filename, "r");
    char line[256];  // Assuming a maximum line length of 255 characters
    int i;

    if (file == NULL) {
        printf("Failed to open the file.\n");
        return 0;  // Exit the function
    }
    while (fgets(line, 256, file) != 0) {
        if (strstr(line, "%RAW%") != 0) {
            i = 0;
            while ((fgets(line, 256, file) != 0) &&
                   (strstr(line, "%END%") == 0)) {
                ++i;
            }
        }
    }
    fclose(file);
    return i;
}

/**
 * Reads the sounding file and fills the arrays with the data. These
 * arrays should have already been allocated to the appropriate size
 * before entering this routine.
 */
void read_snd_file(const char* filename, float* pres, float* hght, float* tmpk,
                   float* dwpk, float* wspd, float* wdir, float* omeg) {
    FILE* file = fopen(filename, "r");

    if (file == NULL) {
        printf("Failed to open the file.\n");
        return;  // Exit the function
    }

    char line[256];  // Assuming a maximum line length of 255 characters
    int i, j;

    while (fgets(line, 256, file) != 0) {
        if (strstr(line, "%RAW%") != 0) {
            i = 0;
            while ((fgets(line, 256, file) != 0) &&
                   (strstr(line, "%END%") == 0)) {
                j = sscanf(line, "%f,%f,%f,%f,%f,%f,%f", &(pres[i]), &(hght[i]),
                           &(tmpk[i]), &(dwpk[i]), &(wdir[i]), &(wspd[i]),
                           &(omeg[i]));
                pres[i] = pres[i] * 100.0;
                tmpk[i] = tmpk[i] + 273.15;
                dwpk[i] = dwpk[i] + 273.15;
				wspd[i] = wspd[i] * 0.544444f;
                ++i;
            }
        }
    }

    fclose(file);  // Close the file
}

/**
 * Loads the sounding data from the file into the Profile object.
 *
 * This function will call get_num_levs_in_file and read_snd_file
 * in order to allocate and fill the arrays. Then, the remaining
 * Sounding variables, such as mixing ratio and virtual temprature,
 * are computed and filled.
 */
void load_sounding(const char* filename, float* pres, float* hght, float* tmpk,
                   float* dwpk, float* mixr, float* vtmp, float* theta,
                   float* thetae, float* mse, float* wdir, float* wspd,
                   float* uwin, float* vwin, float* omeg, const int N) {

    // get the base data arrays and fill them with the data from the file
    read_snd_file(filename, pres, hght, tmpk, dwpk, wspd, wdir, omeg);

    // fill the other arrays with computed data
    for (int i = 0; i < N; ++i) {
        // wind speeds in these files are in knots,
        // but we want them in m/s instead
        wspd[i] = wspd[i];

        float p = pres[i];
        float h = hght[i];
        float t = tmpk[i];
        float d = dwpk[i];
        float spd = wspd[i];
        float dir = wdir[i];

        mixr[i] = sharp_mixratio(p, d);
        vtmp[i] = sharp_virtual_temperature(t, mixr[i], 0.0, 0.0);
        theta[i] = sharp_theta(p, t, 100000.0);
        thetae[i] = sharp_thetae(p, t, d);
        uwin[i] = sharp_u_component(spd, dir);
        vwin[i] = sharp_v_component(spd, dir);

        // need to convert mixing ratio from g/kg to kg/kg
        float spfh = sharp_specific_humidity(mixr[i]); 

        float height_agl = hght[i] - hght[0];
        mse[i] = sharp_moist_static_energy(height_agl, t, spfh);
    }
}

int main(int argc, char** argv) {
    // performance counters
    clock_t start, end;
    double cpu_time_used;

	// filename and number of vertical levels in the file
    const char* filename = "../../data/test_snds/20160524_2302_EF3_37.57_-100.13_108_613967.snd";
    int NZ = get_num_levs_in_file(filename);

	float* pres = malloc(NZ*sizeof(float)); 
    float* hght = malloc(NZ*sizeof(float)); 
    float* tmpk = malloc(NZ*sizeof(float)); 
	float* vtmp = malloc(NZ*sizeof(float)); 
    float* dwpk = malloc(NZ*sizeof(float)); 
    float* wdir = malloc(NZ*sizeof(float)); 
    float* wspd = malloc(NZ*sizeof(float)); 
    float* omeg = malloc(NZ*sizeof(float)); 
	float* mixr = malloc(NZ*sizeof(float)); 
	float* theta = malloc(NZ*sizeof(float)); 
	float* thetae = malloc(NZ*sizeof(float)); 
	float* buoy = malloc(NZ*sizeof(float));
	float* mse = malloc(NZ*sizeof(float));
	float* uwin = malloc(NZ*sizeof(float));
	float* vwin = malloc(NZ*sizeof(float));

    // load the sounding from disk, and compute the remaining
    // variables needed from the base data.
    start = clock();
    load_sounding(filename, pres, hght, tmpk, dwpk, mixr, vtmp, theta, thetae,
                  mse, wdir, wspd, uwin, vwin, omeg, NZ);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and preprocessing time: %fs\n", cpu_time_used);

    // print the contents of these arrays so that we
    // can see we've successfully loaded the data!
    for (int i = 0; i < NZ; ++i) {
        printf(
            "pres[%d] = %f hght[%d] = %f tmpk[%d] = %f dwpk[%d] = %f wdir[%d] "
            "= %f wspd[%d] = %f omeg[%d] = %f\n",
            i, pres[i], i, hght[i], i, tmpk[i], i, dwpk[i], i, wdir[i], i,
            wspd[i], i, omeg[i]);
    }

    // lets compute CAPE and CINH for a surface based, mixed layer,
    // and most unstable parcel. First, we need to create the parcels.
    sharp_Parcel_t* sb_pcl = sharp_Parcel_create();
    sharp_Parcel_t* mu_pcl = sharp_Parcel_create();
    sharp_Parcel_t* ml_pcl = sharp_Parcel_create();

    start = clock();
    // set parcel attributes to the Surface Based Parcel
	sharp_define_parcel(pres, tmpk, dwpk, mixr, theta, thetae, NZ, sb_pcl, 1);
    // Then the Most Unstable parcel
	sharp_define_parcel(pres, tmpk, dwpk, mixr, theta, thetae, NZ, mu_pcl, 3);
    // Finally the mixed layer parcel
	sharp_define_parcel(pres, tmpk, dwpk, mixr, theta, thetae, NZ, ml_pcl, 4);

	// NOTE: That's only if you want to use the convenience functions
	// that compute the 100 mb mixed layer, the most unstable layer, 
	// the forecast surface layer, and so on. If you have a known 
	// starting parcel attribute, you can do the following:
	//
	// float pcl_pres = 950.0;
	// float pcl_tmpk = 300.0;
	// float pcl_dwpk = 295.0;
	// sharp_define_custom_parcel(custom_pcl, pcl_pres, pcl_tmpk, pcl_dwpk);

    // Lift the parcel to compute buoyancy, followed by
    // integrating the CAPE and CINH. The buoyancy data is stored
    // a single array, so we need to do these operations separately
    // for each parcel.
    sharp_lift_parcel_wobf(pres, vtmp, buoy, NZ, sb_pcl);
	sharp_cape_cinh(pres, hght, buoy, NZ, sb_pcl);

    // Now the MU parcel
    sharp_lift_parcel_wobf(pres, vtmp, buoy, NZ, mu_pcl);
	sharp_cape_cinh(pres, hght, buoy, NZ, mu_pcl);

    // Finally do the ML parcel
    sharp_lift_parcel_wobf(pres, vtmp, buoy, NZ, ml_pcl);
	sharp_cape_cinh(pres, hght, buoy, NZ, ml_pcl);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    float sb_cape = sharp_Parcel_get_cape(sb_pcl);
    float sb_cinh = sharp_Parcel_get_cinh(sb_pcl);

    float mu_cape = sharp_Parcel_get_cape(mu_pcl);
    float mu_cinh = sharp_Parcel_get_cinh(mu_pcl);

    float ml_cape = sharp_Parcel_get_cape(ml_pcl);
    float ml_cinh = sharp_Parcel_get_cinh(ml_pcl);

    printf("SB CAPE: %f J/kg\tSB CINH: %f J/kg\n", sb_cape, sb_cinh);
    printf("MU CAPE: %f J/kg\tMU CINH: %f J/kg\n", mu_cape, mu_cinh);
    printf("ML CAPE: %f J/kg\tML CINH: %f J/kg\n", ml_cape, ml_cinh);
    printf("Lifting 3 parcels took %fs\n", cpu_time_used);

    // We MUST clean up any objects we create,
    // otherwise we will leak memory! DON'T LEAK MEMORY!
    sharp_Parcel_delete(sb_pcl);
    sharp_Parcel_delete(mu_pcl);
    sharp_Parcel_delete(ml_pcl);
	free(pres);
    free(hght);
    free(tmpk);
	free(vtmp);
    free(dwpk);
    free(wdir);
    free(wspd);
    free(omeg);
	free(mixr);
	free(theta);
	free(thetae); 
	free(buoy);
	free(mse);
	free(uwin);
	free(vwin);
    return 0;
}
