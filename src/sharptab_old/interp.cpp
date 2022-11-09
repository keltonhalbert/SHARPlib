#include "math.h"
#include "constants.cu"
#include "../utils/grid.h"
#include "profile.cu"

#pragma once
//using namespace std;

// 1D interpolation function for height and pressure coordinates
// requesting only a single value


#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _interp1d(float *x_data, float *y_data, float x_pt, bool is_logp, int *idx, int NX, int NY, int NZ) {
    float x0, x1, x;
    float y0, y1, y;
    int i = idx[0]; int j = idx[1]; int t = idx[3]; // idx[2] is reserved for functions where k is given
    int k, kbot, ktop;
    if (isnan(x_pt)) return nanf("");
    // We need to do a safety check here as 
    // there may be "below ground" levels
    // where data is missing, so find where
    // the data starts. This assumes that the
    // data is well ordered and that there are no
    // missing levels in the middle of the profile

    // find the first valid point and
    // last valid point to handle missing data
    kbot = 0; ktop = NZ-1;
    x0 = x_data[P3(i, j, kbot, t, NX, NY, NZ)];
    x1 = x_data[P3(i, j, ktop, t, NX, NY, NZ)];
    while (isnan(x0)) {
        x0 = x_data[P3(i, j, kbot+1, t, NX, NY, NZ)];
        kbot += 1;
    }
    while (isnan(x1)) {
        x1 = x_data[P3(i, j, ktop-1, t, NX, NY, NZ)];
        ktop -= 1;
    }

    if ((!is_logp) && (( x_pt < x_data[P3(i, j, kbot, t, NX, NY, NZ)] ) || (x_pt > x_data[P3(i, j, ktop, t, NX, NY, NZ)]))) {
        return nanf("");
    }

    if ((is_logp) && ((x_pt > x_data[P3(i, j, kbot, t, NX, NY, NZ)]) || (x_pt < x_data[P3(i, j, ktop, t, NX, NY, NZ)]))) {
        return nanf("");
    }
    if ((x_pt == x_data[P3(i, j, kbot, t, NX, NY, NZ)])) return y_data[P3(i, j, kbot, t, NX, NY, NZ)];
    if ((x_pt == x_data[P3(i, j, ktop, t, NX, NY, NZ)])) return y_data[P3(i, j, ktop, t, NX, NY, NZ)];

    if (is_logp) {
        x_pt = log10(x_pt);

        x = x_pt;
        // use the starting point of the valid data
        k = kbot;
        x1 = log10(x_data[P3(i, j, k, t, NX, NY, NZ)]);
        y1 = y_data[P3(i, j, k, t, NX, NY, NZ)];
        while ((log10(x_data[P3(i, j, k, t, NX, NY, NZ)]) > x) & (k < NZ - 1)) {
            x1 = log10(x_data[P3(i, j, k, t, NX, NY, NZ)]);
            y1 = y_data[P3(i, j, k, t, NX, NY, NZ)];
            k += 1;
        }

        k = ktop;
        x0 = log10(x_data[P3(i, j, k, t, NX, NY, NZ)]);
        y0 = y_data[P3(i, j, k, t, NX, NY, NZ)];
        while ((log10(x_data[P3(i, j, k, t, NX, NY, NZ)]) <= x) & (k > 0)) {
            x0 = log10(x_data[P3(i, j, k, t, NX, NY, NZ)]);
            y0 = y_data[P3(i, j, k, t, NX, NY, NZ)];
            k -= 1;
        }
    }

    else {
        x = x_pt;

        k = kbot;
        x0 = x_data[P3(i, j, k, t, NX, NY, NZ)];
        y0 = y_data[P3(i, j, k, t, NX, NY, NZ)];
        while ((x_data[P3(i, j, k, t, NX, NY, NZ)] <= x) & (k < NZ-1)) {
            x0 = x_data[P3(i, j, k, t, NX, NY, NZ)];
            y0 = y_data[P3(i, j, k, t, NX, NY, NZ)];
            k += 1;
        }

        k = ktop;
        x1 = x_data[P3(i, j, k, t, NX, NY, NZ)];
        y1 = y_data[P3(i, j, k, t, NX, NY, NZ)];
        while ((x_data[P3(i, j, k, t, NX, NY, NZ)] > x) & (k > 0)) {
            x1 = x_data[P3(i, j, k, t, NX, NY, NZ)];
            y1 = y_data[P3(i, j, k, t, NX, NY, NZ)];
            k -= 1;
        }
    }

    y = y0 * ( 1 - (x - x0) / (x1 - x0)) + y1 * ( (x - x0) / (x1 - x0));
    return y;
}

// 1D interpolation function for height and pressure coordinates
// requesting only an array of values
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float* _interp1d(float *x_data, float *y_data, float *x_pts, bool is_logp, int *idx, int NX, int NY, int NZ, int nPts) {
    float *y = new float[nPts];
    float x;

    // loop over the interpolation points
    // and do the interpolation at that point
    for (int ptidx = 0; ptidx < nPts; ptidx++) {

        x = x_pts[ptidx];
        y[ptidx] = _interp1d( x_data, y_data, x, is_logp, idx, NX, NY, NZ);
    }

    return y;
}

// interpolate the height in agl coordinates
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _interp_to_agl(float *pres_arr, float *hght_arr, float pval, int *idx, int NX, int NY, int NZ) {
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    // find the first valid point and
    // last valid point to handle missing data
    int kbot = 0;
    while ((isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) && (kbot < NZ) ) {
        kbot += 1;
    }
    int ktop = NZ-1;
    while ((isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) && (ktop > 0)) {
        ktop -= 1;
    }
    if ((ktop == 0) || (kbot == NZ) || (isnan(pval)) || (pval == MISSING)) return nanf("");
    return _interp1d(pres_arr, hght_arr, pval, true, idx, NX, NY, NZ) - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)];
}

// interpolate the height in msl coordinates
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _interp_to_msl(float *hght_arr, float hval, int *idx, int NX, int NY, int NZ) {
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    int kbot = 0;
    while ((isnan(hght_arr[P3(i, j, kbot, t, NX, NY, NZ)])) && (kbot < NZ) ) {
        kbot += 1;
    }
    int ktop = NZ-1;
    while ((isnan(hght_arr[P3(i, j, ktop, t, NX, NY, NZ)])) && (ktop > 0)) {
        ktop -= 1;
    }
    if ((ktop == 0) || (kbot == NZ) || (isnan(hval)) || (hval == MISSING)) return nanf("");
    return hval + hght_arr[P3(i, j, kbot, t, NX, NY, NZ)];
} 

