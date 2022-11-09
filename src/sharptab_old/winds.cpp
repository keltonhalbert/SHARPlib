#include <iostream>
#include "../utils/math.cu"
#include "constants.cu"
#include "profile.cu"
//using namespace std;

#pragma once

#ifdef __CUDACC__
__host__ __device__ 
#endif
void _mean_wind(profile *prof, float pbot, float ptop, float *out_arr, int *idx){
    float *pres_arr = prof->pres;
    float *u_arr = prof->uwin;
    float *v_arr = prof->vwin;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int dp = -1;
    // int stu = 0;
    // int stv = 0;
    float ps[1000];
    float u[1000];
    float v[1000];

    int num_p = ((pbot - (ptop + dp))/(-1*dp));

    for (int ip = 0; ip < num_p; ++ip) {
        ps[ip] = pbot + dp*ip;
        u[ip] = _interp1d(pres_arr, u_arr, ps[ip], true, idx, NX, NY, NZ);
        v[ip] = _interp1d(pres_arr, v_arr, ps[ip], true, idx, NX, NY, NZ);
    }

    float umean = weighted_average(u, ps, num_p);
    float vmean = weighted_average(v, ps, num_p);
    out_arr[0] = umean;
    out_arr[1] = vmean;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
void _mean_wind_npw(profile *prof, float pbot, float ptop, float *out_arr, int *idx){

    float *pres_arr = prof->pres;
    float *u_arr = prof->uwin;
    float *v_arr = prof->vwin;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int dp = -1;
    float ps[1100];
    float weights[1100];
    float u[1100];
    float v[1100];

    int num_p = ((pbot - (ptop + dp))/(-1*dp));

    for (int ip = 0; ip < num_p; ++ip) {
        ps[ip] = pbot + dp*ip;
        weights[ip] = 1.;
        u[ip] = _interp1d(pres_arr, u_arr, ps[ip], true, idx, NX, NY, NZ);
        v[ip] = _interp1d(pres_arr, v_arr, ps[ip], true, idx, NX, NY, NZ);
    }

    float umean = weighted_average(u, weights, num_p);
    float vmean = weighted_average(v, weights, num_p);
    out_arr[0] = umean;
    out_arr[1] = vmean;
}


#ifdef __CUDACC__
__host__ __device__ 
#endif
void _helicity(profile *prof, float lower, float upper, float stu, float stv, float *out_arr, int *idx){

    float *pres_arr = prof->pres;
    float *hght_arr = prof->hght;
    float *u_arr = prof->uwin;
    float *v_arr = prof->vwin;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int dp = -1;
    int ip_u;
    float plower, pupper;
    float layer;
    float phel = 0.f;
    float nhel = 0.f;

    float ps[1000];
    float u[1000];
    float v[1000];


    if (lower != upper){
        lower = _interp_to_msl(hght_arr, lower, idx, NX, NY, NZ);
        upper = _interp_to_msl(hght_arr, upper, idx, NX, NY, NZ);
        plower = _interp1d(hght_arr, pres_arr, lower, false, idx, NX, NY, NZ);
        pupper = _interp1d(hght_arr, pres_arr, upper, false, idx, NX, NY, NZ);

        if (isnan(plower) || isnan(pupper) || (plower == MISSING) || (pupper == MISSING) || \
           (stu == MISSING) || (stv == MISSING) || (isnan(stu)) || (isnan(stv))) {
            out_arr[0] = nanf(""); out_arr[1] = nanf(""); out_arr[2] = nanf("");
            return;
        }

        // calculate the values in the interpolation interval
        int num_p = ((plower - (pupper + dp))/(-1*dp));
        for (int ip = 0; ip < num_p; ++ip) {
            ps[ip] = plower + dp*ip;
            // go ahead and subtract out the storm motion during this step
            u[ip] = _interp1d(pres_arr, u_arr, ps[ip], true, idx, NX, NY, NZ) - stu;
            v[ip] = _interp1d(pres_arr, v_arr, ps[ip], true, idx, NX, NY, NZ) - stv;
        }

        // loop again to do the layer sum. While the loop is
        // kind of redundant, it prevents redundant calls to interp1d,
        // which means there should be less overhead
        for (int ip_l = 0; ip_l < num_p - 1; ++ip_l) {
            ip_u = ip_l + 1;
            layer = (u[ip_u] * v[ip_l]) - (u[ip_l] * v[ip_u]);
            if (layer > 0.f) phel += layer;
            else if (layer < 0.f) nhel += layer; 
        }

    }
    else{
         phel = nanf("");
         nhel = nanf("");
    }


    out_arr[0] = phel + nhel;
    out_arr[1] = phel;
    out_arr[2] = nhel;

}


#ifdef __CUDACC__
__host__ __device__ 
#endif
void _wind_shear(profile *prof, float pbot, float ptop, float *out_arr, int *idx){

    float *pres_arr = prof->pres;
    float *u_arr = prof->uwin;
    float *v_arr = prof->vwin;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    float ubot;
    float utop;
    float vbot;
    float vtop;
    float shu;
    float shv;
    
    if ( (pbot == MISSING) || (ptop == MISSING) || isnan(ptop) || isnan(pbot) ) {
        out_arr[0] = 0.f; out_arr[1] = 0.f;
        return;
    }


    ubot = _interp1d(pres_arr, u_arr, pbot, true, idx, NX, NY, NZ);
    utop = _interp1d(pres_arr, u_arr, ptop, true, idx, NX, NY, NZ);
    vbot = _interp1d(pres_arr, v_arr, pbot, true, idx, NX, NY, NZ);
    vtop = _interp1d(pres_arr, v_arr, ptop, true, idx, NX, NY, NZ);

    shu = utop - ubot;
    shv = vtop - vbot;

    out_arr[0] = shu; out_arr[1] = shv;


}

#ifdef __CUDACC__
__host__ __device__ 
#endif
void _non_parcel_bunkers_motion(profile *prof, float *out_arr, int *idx) {

    float *pres_arr = prof->pres;
    float *hght_arr = prof->hght;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float d = 7.5f; // deviation value empirically derived as 7.5 m/s

    int kbot = 0;
    while ((isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) && (kbot < NZ) ) {
        kbot += 1;
    }
    if (kbot == NZ) {
        out_arr[0] = nanf(""); out_arr[1] = nanf("");
        out_arr[2] = nanf(""); out_arr[3] = nanf("");
        return;
    }

    float msl500m = _interp_to_msl(hght_arr, 500.f, idx, NX, NY, NZ);
    float msl5500m = _interp_to_msl(hght_arr, 5500.f, idx, NX, NY, NZ);
    float msl6000m = _interp_to_msl(hght_arr, 6000.f, idx, NX, NY, NZ);
    float p500m = _interp1d(hght_arr, pres_arr, msl500m, false, idx, NX, NY, NZ);
    float p5500m = _interp1d(hght_arr, pres_arr, msl5500m, false, idx, NX, NY, NZ);
    float p6000m = _interp1d(hght_arr, pres_arr, msl6000m, false, idx, NX, NY, NZ);


    // SFC-6km mean wind
    float mn_wind_500m[2];
    float mn_wind_6000m[2];
    float mn_wind_sfc_6km[2];
    _mean_wind(prof, pres_arr[P3(i, j, kbot, t, NX, NY, NZ)], p500m, mn_wind_500m, idx);
    _mean_wind(prof, p5500m, p6000m, mn_wind_6000m, idx);
    _mean_wind_npw(prof, pres_arr[P3(i, j, kbot, t, NX, NY, NZ)], p6000m, mn_wind_sfc_6km, idx);

    float shr[2];
    shr[0] = mn_wind_6000m[0] - mn_wind_500m[0];
    shr[1] = mn_wind_6000m[1] - mn_wind_500m[1];

    // bunkers right and left motions
    float mag = sqrt(pow(shr[0], 2.f) + pow(shr[1], 2.f));
    float tmp = d / mag;
    float rstu = mn_wind_sfc_6km[0] + (tmp * shr[1]);
    float rstv = mn_wind_sfc_6km[1] - (tmp * shr[0]);
    float lstu = mn_wind_sfc_6km[0] - (tmp * shr[1]);
    float lstv = mn_wind_sfc_6km[1] + (tmp * shr[0]);

    out_arr[0] = rstu; out_arr[1] = rstv;
    out_arr[2] = lstu; out_arr[3] = lstv;
}

/* Calculates the critical angle (degrees) as specified by Esterheld and Giuliano (2008).
   If the critical angle is 90 degrees, this indicates that the lowest 500 meters of 
   the storm is experiencing pure streamwise vorticity. */
#ifdef __CUDACC__
__host__ __device__ 
#endif
float _critical_angle(profile *prof, float stu, float stv, int *idx ) {
    float *pres_arr = prof->pres;
    float *hght_arr = prof->hght;
    float *u_arr = prof->uwin;
    float *v_arr = prof->vwin;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float vec1_u, vec2_u;
    float vec1_v, vec2_v;
    float vec1_mag, vec2_mag;
    float dot, angle;
    int kbot = 0;
    // we need to find the valid first k level of
    // data since there can be below-ground missing data
    while ((isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) && (kbot < NZ) ) {
        kbot+=1;
    }
    if (kbot == NZ) {
        return nanf("");
    }
    // convert 500 meters AGL to it's corresponding
    // height above sea level, and then get the pressure
    // value at that height
    float lower = _interp_to_msl(hght_arr, 500., idx, NX, NY, NZ);
    float p500m = _interp1d(hght_arr, pres_arr, lower, false, idx, NX, NY, NZ);
 
    // Get the 500 meter u and v wind components
    float u500m = _interp1d(pres_arr, u_arr, p500m, true, idx, NX, NY, NZ);
    float v500m = _interp1d(pres_arr, v_arr, p500m, true, idx, NX, NY, NZ);

    // get the surface u and v wind components
    float usfc = u_arr[P3(i, j, kbot, t, NX, NY, NZ)]; float vsfc = v_arr[P3(i, j, kbot, t, NX, NY, NZ)];

    vec1_u = u500m - usfc; vec1_v = v500m - vsfc;
    vec2_u = stu - usfc; vec2_v = stv - vsfc;

    // compute the vector magnitudes
    vec1_mag = sqrt( powf(vec1_u, 2.) + powf(vec1_v, 2.) ); 
    vec2_mag = sqrt( powf(vec2_u, 2.) + powf(vec2_v, 2.) );

    // compute the dot product and get the angle 
    // between the two vectors, then convert
    // from radians to degrees and return it
    dot = vec1_u * vec2_u + vec1_v * vec2_v;
    angle =  (180. / PI ) * acos( dot / (vec1_mag * vec2_mag) );
    if ((i == 176) && (j == 112)) printf("%d %d %f %f\n",i, j, lower, p500m);
    return angle;
} 

#ifdef __CUDACC__
__host__ __device__ 
#endif
void comp2vec(float *uv, float *spd_dir) {
    const float PI = 3.14159265;
    const float TOL = 1e-10; 
    float u = uv[0]; float v = uv[1];
    if ((isnan(u)) || (isnan(v))) {
        spd_dir[0] = nanf(""); spd_dir[1] = nanf("");
    }
    float wdir = atan2(-u, -v) * 180 / PI;

    if (wdir < 0) wdir += 360.;
    if (fabs(wdir) < TOL) wdir = 0.;

    spd_dir[0] = sqrt(powf(u, 2.) + powf(v, 2.));
    spd_dir[1] = wdir;
    return;
}

