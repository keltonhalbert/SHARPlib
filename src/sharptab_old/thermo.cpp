#include <iostream>
#include "math.h"
#include "./constants.cu"
#include "../utils/grid.h"
#include "profile.cu"

#pragma once
//using namespace std;

// implementation of the Wobus Function for computing
// moist adiabats.
// Input: t, degrees Celsius
// Output: Correction to theta (C) for the calculation
// of saturated potential temperature
#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _wobf(float t) {
    float npol;
    float ppol;
    float tmp;
    float correction;
    
    tmp = t - 20.f;
    if (isnan(t)) return t;
    // check for missing data
    if (t == MISSING) return t;


    if (tmp <= 0.f) {
        npol = 1.f + tmp * (-8.841660499999999e-3f + tmp * ( 1.4714143e-4f + tmp * (-9.671989000000001e-7f + tmp * (-3.2607217e-8f + tmp * (-3.8598073e-10f)))));
        npol = 15.13f / (powf(npol, 4.f));
        correction =  npol;
     }

    else if (tmp > 0.f) {
        ppol = tmp * (4.9618922e-07f + tmp * (-6.1059365e-09f + tmp * (3.9401551e-11f + tmp * (-1.2588129e-13f + tmp * (1.6688280e-16f)))));
        ppol = 1.f + tmp * (3.6182989e-03f + tmp * (-1.3603273e-05f + ppol));
        ppol = (29.93f / powf(ppol,4.f)) + (0.96f * tmp) - 14.8f;
        correction = ppol;

    }

    return correction;
}

// compute the potential temperature of a parcel
// when p2=1000. hPa
#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _theta(float p, float t, float p2) {
    if ((p == MISSING) || (t == MISSING) || (p2 == MISSING)) return nanf("");
    return ((t + ZEROCNK) * powf((p2/p), ROCP)) - ZEROCNK;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _lcltemp(float t, float td) {
    float s; 
    float dlt;
    if ((t == MISSING) || (td == MISSING)) return nanf("");

    s = t - td;
    dlt = s * (1.2185f + 0.001278f * t + s * (-0.00219f + 1.173e-5f * s - 0.0000052f * t));
    return t - dlt;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _thalvl(float theta, float t) {
    if ((theta == MISSING) || (t == MISSING)) return nanf("");
    t = t + ZEROCNK;
    theta = theta + ZEROCNK;
    return 1000.f / ( powf( ( theta/t ), ( 1.f/ROCP ) ) );
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _vappres(float t) {
    float pol;
    if (t == MISSING) return nanf("");

    pol = t * (1.1112018e-17f + (t * -3.0994571e-20f));
    pol = t * (2.1874425e-13f + (t * (-1.789232e-15f + pol)));
    pol = t * (4.3884180e-09f + (t * (-2.988388e-11f + pol)));
    pol = t * (7.8736169e-05f + (t * (-6.111796e-07f + pol)));
    pol = 0.99999683f + (t * (-9.082695e-03f + pol));
    return 6.1078f / powf(pol, 8.f);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _mixratio(float p, float t) {
    if ((p ==  MISSING) || (t == MISSING)) return nanf("");

    float x;
    float wfw;
    float fwesw;

    x = 0.02f * (t - 12.5f + (7500.f / p));
    wfw = 1.f + (0.0000045f * p) + (0.0014f * x * x);
    fwesw = wfw * _vappres(t);
    return 621.97f * (fwesw / (p - fwesw));
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _virtemp(float p, float t, float td) {
    float tk;
    float w;
    float vt;
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) return nanf("");

    tk = t + ZEROCNK;
    w = 0.001f * _mixratio(p, td);
    vt = (tk * (1.f + w / eps) / (1.f + w)) - ZEROCNK;
    if (isnan(vt)) return t;
    else return vt;
}

#ifdef __CUDACC__
__global__ void virtemp(profile *prof, float *out_arr) {

    float *pres_arr = prof->pres;
    float *tmpc_arr = prof->tmpc;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i, j, k;
    int idx3D;
    float vtmp;
    
    i = (blockIdx.x * blockDim.x) + threadIdx.x;
    j = (blockIdx.y * blockDim.y) + threadIdx.y;
    k = (blockIdx.z * blockDim.z) + threadIdx.z;
    if ((j < NY) & (i < NX) & (k < NZ)) {

        idx3D = P3(i, j, k, NX, NY);
        vtmp = _virtemp(pres_arr[idx3D], tmpc_arr[idx3D], dwpc_arr[idx3D]);
        out_arr[idx3D] = vtmp;
    }

    return;
}
#endif 

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _temp_at_mixrat(float w, float p) {
    float x;
    if ((w == MISSING) || (p == MISSING)) return nanf("");
    
    x = log10(w * p / (622.f + w));
    x = ( powf(10.f, ((c1 * x) + c2 )) - c3 + (c4 * powf((powf(10.f, (c5 * x)) - c6), 2.f))) - ZEROCNK; 
    return x;
}

#ifdef __CUDACC__
__global__ void temp_at_mixrat(profile *prof, float *out_arr) {

    float *mixr_arr = prof->mixr;
    float *pres_arr = prof->pres;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i, j, k;
    int idx3D;
    float dwpc;
    
    i = (blockIdx.x * blockDim.x) + threadIdx.x;
    j = (blockIdx.y * blockDim.y) + threadIdx.y;
    k = (blockIdx.z * blockDim.z) + threadIdx.z;
    if ((j < NY) & (i < NX) & (k < NZ)) {

        idx3D = P3(i, j, k, NX, NY);
        dwpc = _temp_at_mixrat(mixr_arr[idx3D], pres_arr[idx3D]);
        out_arr[idx3D] = dwpc;
    }

    return;
}
#endif 


#ifdef __CUDACC__
__host__ __device__ 
#endif
inline void _drylift(float p, float t, float td, float *out) {
    float t2;
    float p2;
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) {
        out[0] = MISSING; out[1] = MISSING;
        return;
    }
    //float out[2];
    //float *out = new float[2];

    t2 = _lcltemp(t, td);
    p2 = _thalvl(_theta(p, t, 1000.f), t2);
    out[0] = t2; out[1] = p2;
    //return out;
}

// Calculate the temperature of a saturated parcel
// when lifted to a new pressure level.
// Input: pres, hPa | thetam, degrees Celcius
// Output: Temperature of the saturated parcel at
// the new pressure level. 
#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _satlift(float pres, float thetam) {
    if ((pres == MISSING) || (thetam == MISSING)) return nanf("");
    if ( ( fabs( pres - 1000.f ) - 0.001f ) <= 0.0f ) {
        return thetam;
    }

    float pwrp;
    float t1;
    float t2;
    float e1;
    float e2;
    float eor = 999.0f;
    float rate;
    while ( ( fabs(eor) - 0.01f ) > 0.f ) {
        if (eor == 999.0f) {
            pwrp = powf((pres/1000.f), ROCP);
            t1 = (thetam + ZEROCNK) * pwrp - ZEROCNK;
            e1 = _wobf(t1) - _wobf(thetam);
            rate = 1.f;
        }
        else {
            rate = (t2 - t1) / (e2 - e1);
            t1 = t2;
            e1 = e2;
        }
        t2 = t1 - (e1 * rate);
        e2 = (t2 + ZEROCNK) / pwrp - ZEROCNK;
        e2 += _wobf(t2) - _wobf(e2) - thetam;
        eor = e2 * rate;
    }
    return t2 - eor;
}

// lift the parcel moist adiabatically to it's new level
#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _wetlift(float p, float t, float p2) {
    float thta;
    float thetam;
    if ((p == MISSING) || (t == MISSING) || (p2 == MISSING)) return nanf("");

    thta = _theta(p, t, 1000.f);
    if (isnan(thta)) return nanf("");
    thetam = thta - _wobf(thta) + _wobf(t);
    return _satlift(p2, thetam);
} 

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _lifted(float p, float t, float td, float lev) {
    float dry[2];
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) return nanf("");
    _drylift(p, t, td, dry);
    return _wetlift(dry[1], dry[0], lev);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _wetbulb(float p, float t, float td) {
    float dry[2];
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) return nanf("");
    _drylift(p, t, td, dry);
    return _wetlift(dry[1], dry[0], p);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _thetaw(float p, float t, float td) {
    float dry[2];
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) return nanf("");
    _drylift(p, t, td, dry);
    return _wetlift(dry[1], dry[0], 1000.f);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _thetae(float p, float t, float td) {
    float dry[2];
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) return nanf("");
    _drylift(p, t, td, dry);
    return _theta(100.f, _wetlift(dry[1], dry[0], 100.f), 1000.f);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _relh(float p, float t, float td) {
    if ((p == MISSING) || (t == MISSING) || (td == MISSING)) return nanf("");
    return 100.f * _mixratio(p, td) / _mixratio(p, t);
}

#ifdef __CUDACC__
__host__ __device__ 
#endif
inline float _clausius_clapeyron_t(float p, float t, float w) {
    if ((p == MISSING) || (t == MISSING) || (w == MISSING)) return nanf("");

    float e0 = 6.1173f;
    float Rv = 461.50f;
    float Lv = 2.501f * pow(10.f, 6.f);
    float K1 = Lv / Rv;
    float K2 = 1.f/ZEROCNK;
    float K3 = 1.f/(t+ZEROCNK);

    float es = e0 * exp(K1*(K2-K3));
    float rh = w / (_mixratio(p, t));
    float e = rh * es;
    float K1_inv = Rv/Lv;
    float K4 = log(e/e0);
    float term1 = K2 - (K1_inv * K4);
    return (1.f/term1) - ZEROCNK;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
inline float _ctok(float t) {
    if (t == MISSING) return MISSING;
    return t + ZEROCNK;
}
