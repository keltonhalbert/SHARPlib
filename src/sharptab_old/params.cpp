#include "thermo.cu"
#include "winds.cu"
#include "interp.cu"
#include "constants.cu"
#include "../utils/math.cu"
#include "profile.cu"
#include <iostream>
#include <stdio.h>

#pragma once

// This will serve as out C++ Parcel object that
// retains information throughout the parcel lift.
// This also stores the parcel lifting level
// that gets set by defineParcel
struct Parcel {

    float pres = MISSING; // the beginning pressure 
    float tmpc = MISSING; // the beginning temperature
    float dwpc = MISSING; // the beginning dewpoint

    float lclpres = MISSING; // LCL pressure (mb)
    float lclhght = MISSING; // LCL height (m)
    float lfcpres = MISSING; // LFC pressure (mb)
    float lfchght = MISSING; // LFC height (m)
    float elpres = MISSING; // EL pressure (mb)
    float elhght = MISSING; // EL height (m)
    float mplpres = MISSING; // maximum parcel level (mb)
    float mplhght = MISSING; // maximum parcel level (m)
    float effpbot = MISSING; // effective inflow layer bottom (hPa)
    float effptop = MISSING; // effective inflow layer top (hPa)
    float bplus = MISSING; // CAPE
    float bminus = MISSING; // CINH
    float b3km = MISSING; // 0-3km CAPE
    float b6km = MISSING; // 0-6km CAPE
    float li5 = MISSING; // 500mb Lifted Index
    float li3 = MISSING; // 300mb Lifted Index
    float brnshear = MISSING; // Bulk Richardson Number Shear
    float brnu = MISSING; // Bulk Richardson Number U
    float brnv = MISSING; // Bulk Richardson Number V
    float brn = MISSING; // Bulk Richardson Number
    float limax = MISSING; // Maximum Lifted Index
    float limaxpres = MISSING; // Pressure Level of max LI (mb)
    float cap = MISSING; // Cap Strength (C)
    float cappres = MISSING; // Cap Pressrure (mb)
    float bmin = MISSING; // Buoyancy Minimum
    float bminpres = MISSING; // Pressure level of bmin

    float rstu = MISSING; // bunkers right motion u (m/s)
    float rstv = MISSING; // bunkers right motion v
    float lstu = MISSING;
    float lstv = MISSING;
};

#ifdef __CUDACC__
__host__ __device__
#endif
void clear_pcl(Parcel pcl) {
    pcl.pres = MISSING; // the beginning pressure 
    pcl.tmpc = MISSING; // the beginning temperature
    pcl.dwpc = MISSING; // the beginning dewpoint

    pcl.lclpres = MISSING; // LCL pressure (mb)
    pcl.lclhght = MISSING; // LCL height (m)
    pcl.lfcpres = MISSING; // LFC pressure (mb)
    pcl.lfchght = MISSING; // LFC height (m)
    pcl.elpres = MISSING; // EL pressure (mb)
    pcl.elhght = MISSING; // EL height (m)
    pcl.mplpres = MISSING; // maximum parcel level (mb)
    pcl.mplhght = MISSING; // maximum parcel level (m)
    pcl.effpbot = MISSING; // effective inflow layer bottom (hPa)
    pcl.effptop = MISSING; // effective inflow layer top (hPa)
    pcl.bplus = MISSING; // CAPE
    pcl.bminus = MISSING; // CINH
    pcl.b3km = MISSING; // 0-3km CAPE
    pcl.b6km = MISSING; // 0-6km CAPE
    pcl.li5 = MISSING; // 500mb Lifted Index
    pcl.li3 = MISSING; // 300mb Lifted Index
    pcl.brnshear = MISSING; // Bulk Richardson Number Shear
    pcl.brnu = MISSING; // Bulk Richardson Number U
    pcl.brnv = MISSING; // Bulk Richardson Number V
    pcl.brn = MISSING; // Bulk Richardson Number
    pcl.limax = MISSING; // Maximum Lifted Index
    pcl.limaxpres = MISSING; // Pressure Level of max LI (mb)
    pcl.cap = MISSING; // Cap Strength (C)
    pcl.cappres = MISSING; // Cap Pressrure (mb)
    pcl.bmin = MISSING; // Buoyancy Minimum
    pcl.bminpres = MISSING; // Pressure level of bmin

    pcl.rstu = MISSING; // bunkers right motion u (m/s)
    pcl.rstv = MISSING; // bunkers right motion v
    pcl.lstu = MISSING;
    pcl.lstv = MISSING;
}

// compute the most unstable parcel level by calculating the maximum thetae
// in the defined layer
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float most_unstable_level(profile *prof, float pbot, float ptop, int *idx) {

    // these are arrays that will hold data from 
    // pbot to ptop in dp increments - I figured 
    // being of size 1000 would be more than enough
    // since theoretically it's the size of the whole
    // column if the surface is at 1000mb and ptop is 0.
    float *pres_arr = prof->pres;
    float *tmpc_arr = prof->tmpc;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    int k;
    float p[1000];
    float tmp[1000];
    float d[1000];
    float mt;
    float max = -999.0f;
    int mi = -1;
    float dp = -1.f;

    // we need to check and make sure that the
    // ptop requested isn't outisde of the range
    // of data we actually posess. If it is outside
    // the range, just use the data we have.
    k = NZ-1;
    float datatop = pres_arr[P3(i, j, k, t, NX, NY, NZ)];
    while (datatop == MISSING) {
        datatop = pres_arr[P3(i, j, k-1, t, NX, NY, NZ)];
        k -= 1;
    }
    if (ptop < datatop) {
        ptop = datatop;
    }
    k = 0;
    float databot = pres_arr[P3(i, j, k, t, NX, NY, NZ)];
    while (databot == MISSING) {
        databot = pres_arr[P3(i, j, k+1, t, NX, NY, NZ)];
        k+=1;
    }
    if (pbot > databot) {
        pbot = databot;
    }

    // iterate from pbot to ptop in dp increments
    int num_p = (int)((pbot - (ptop + dp)) / (-1*dp));
    for (int i = 0; i < num_p; ++i) {
        p[i] = pbot +dp*i;
        // get the temperature and dewpoint at the level
        tmp[i] = _interp1d(pres_arr, tmpc_arr, p[i], true, idx, NX, NY, NZ);
        d[i] = _interp1d(pres_arr, dwpc_arr, p[i], true, idx, NX, NY, NZ);
        // essentially computing the thetae
        mt = _thetae(p[i], tmp[i], d[i]);

        // check if it's the max thetae,
        // store it's value and store
        // its level
        if ((mt > max) && !(isnan(mt))) {
            max = mt;
            mi = i;
        }
    }

    // return the pressure level
    // of the most unstable parcel
    return p[mi]; 
}

// calculate the pressure weighted mean theta between 
// pbot and ptop by interpolating in dp=-1 increments
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float mean_theta(profile *prof, float pbot, float ptop, int *idx) {

    float *pres_arr = prof->pres;
    float *tmpc_arr = prof->tmpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float p[1000];
    float tmp[1000];
    float thta[1000];
    float dp = -1.f;
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }
    if (pbot > pres_arr[P3(i, j, kbot, t, NX, NY, NZ)]) pbot = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    if (ptop < pres_arr[P3(i, j, ktop, t, NX, NY, NZ)]) ptop = pres_arr[P3(i, j, ktop, t, NX, NY, NZ)];

    int num_p = (int)((pbot - (ptop + dp))/(-1*dp));
    for (int ip = 0; ip < num_p; ++ip) {
        p[ip] = pbot + dp*ip;
        tmp[ip] = _interp1d(pres_arr, tmpc_arr, p[ip], true, idx, NX, NY, NZ);
        thta[ip] = _theta(p[ip], tmp[ip], 1000.f);
    }

    return weighted_average(thta, p, num_p);

}


/* Calculate the Precipitable Water Vapor in Inches */
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float pwat(profile *prof, float pbot, float ptop, int *idx) {
    // these are arrays that will hold data from 
    // pbot to ptop in dp increments - I figured 
    // being of size 1000 would be more than enough
    // since theoretically it's the size of the whole
    // column if the surface is at 1000mb and ptop is 0.
    float *pres_arr = prof->pres;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float p[1000];
    float w[1000];
    float d[1000];
    float dp = -1.f;
    float pw = 0.0;
    // we need to check and make sure that the
    // ptop requested isn't outisde of the range
    // of data we actually posess. If it is outside
    // the range, just use the data we have.
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }
    float datatop = pres_arr[P3(i, j, ktop, t, NX, NY, NZ)];
    float databot = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    if (ptop < datatop) {
        ptop = datatop;
    }
    if (pbot > databot) {
        pbot = databot;
    }

    // compute the column mixing ratio in 
    // interpolated dp intervals
    int num_p = (int)((pbot - (ptop + dp)) / (-1*dp));
    for (int i = 0; i < num_p; ++i) {
        p[i] = pbot + dp*i;
        // get the temperature and dewpoint at the level
        d[i] = _interp1d(pres_arr, dwpc_arr, p[i], true, idx, NX, NY, NZ);
        w[i] = _mixratio(p[i], d[i]);
    }

    // integrate mixing ratio to get
    // pwat
    for (int i = 1; i < num_p; ++i) {
        pw += ( ( ( w[i-1] + w[i] ) / 2. ) * ( p[i-1] - p[i] ) ) * 0.00040173;
    }
    return pw;
}

// compute the pressure weighted mean mixing ratio between
// pbot and ptop by interpolating in dp=-1 increments
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float mean_mixratio(profile *prof, float pbot, float ptop, int *idx) {

    float *pres_arr = prof->pres;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float p[1000];
    float d[1000];
    float m[1000];
    float dp = -1.f;
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }
    if (pbot > pres_arr[P3(i, j, kbot, t, NX, NY, NZ)]) pbot = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    if (ptop < pres_arr[P3(i, j, ktop, t, NX, NY, NZ)]) ptop = pres_arr[P3(i, j, ktop, t, NX, NY, NZ)];

    int num_p = (int)((pbot - (ptop + dp))/(-1*dp));
    for (int ip = 0; ip < num_p; ++ip) {
        p[ip] = pbot + dp*ip;
        d[ip] = _interp1d(pres_arr, dwpc_arr, p[ip], true, idx, NX, NY, NZ);
        m[ip] = _mixratio(p[ip], d[ip]);
    }

    float mmr = weighted_average(m, p, num_p);
    return mmr;

}



// sets a parcels lifting paramaters to the surface
#ifdef __CUDACC__
__host__ __device__ 
#endif 
void _sfc(profile *prof, Parcel *pcl, int *idx) {
    float *pres_arr = prof->pres;
    float *tmpc_arr = prof->tmpc;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }
    pcl->pres = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    pcl->tmpc = tmpc_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    pcl->dwpc = dwpc_arr[P3(i, j, kbot, t, NX, NY, NZ)];
}

// sets a parcels lifting paramaters to the most unstable parcel
#ifdef __CUDACC__
__host__ __device__ 
#endif 
void _mu(profile *prof, Parcel *pcl, float ptop, int *idx) {

    float *pres_arr = prof->pres;
    float *tmpc_arr = prof->tmpc;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }
    if (ptop < pres_arr[P3(i, j, ktop, t, NX, NY, NZ)]) ptop = pres_arr[P3(i, j, ktop, t, NX, NY, NZ)];
    float pbot = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    float mup = most_unstable_level(prof, pbot, pbot-ptop, idx);

    float mut = _interp1d(pres_arr, tmpc_arr, mup, true, idx, NX, NY, NZ);
    float mud = _interp1d(pres_arr, dwpc_arr, mup, true, idx, NX, NY, NZ);
    pcl->pres = mup; 
    pcl->tmpc = mut;
    pcl->dwpc = mud; 
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
void _ml(profile *prof, Parcel *pcl, float pdepth, int *idx) {

    float *pres_arr = prof->pres;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }
    float pbot = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    float ptop = pbot - pdepth;
    // make sure ptop is in the domain 
    if (ptop < pres_arr[P3(i, j, ktop, t, NX, NY, NZ)]) ptop = pres_arr[P3(i, j, ktop, t, NX, NY, NZ)];


    float mtheta = mean_theta(prof, pbot, ptop, idx);
    float mmr = mean_mixratio(prof, pbot, ptop, idx);

    pcl->pres = pbot;
    pcl->tmpc = _theta(1000.f, mtheta, pbot);
    pcl->dwpc = _temp_at_mixrat(mmr, pbot);

}

// definteParcel is used to set whether or not this is the surface, most unstable, or mixed layer parcel
#ifdef __CUDACC__
__host__ __device__ 
#endif 
void defineParcel(profile *prof, Parcel *pcl, int *idx, int flag) {

    if (flag == 1) {
        _sfc(prof, pcl, idx);
    }

    if (flag == 3) {
        _mu(prof, pcl, 300.f, idx);
    }

    if (flag == 4) {
        _ml(prof, pcl, 100.f, idx);
    }
    if (flag == 5) {
        // the parcel has already been defined - do nothing
    }

}


// device function for computing CAPE given 3D grid arrays
// for pressure, temperature, and dewpoint, with NZ being the
// number of vertical levels, i.e. len(pres), nY being the latitudinal
// dimension, and nX being the longitudinal dimension

#ifdef __CUDACC__
__host__ __device__ 
#endif 
void _cape(profile *prof, Parcel *pcl, int *idx, int flag) {
    float *pres_arr = prof->pres;
    float *hght_arr = prof->hght;
    float *tmpc_arr = prof->tmpc;
    float *vtmp_arr = prof->vtmp;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float pe1, pe2, pe3, pe4, pelast;
    float h1, h2, h3, h4;
    float tp1, tp2, tp3;
    float te1, te2, te3;
    float dry_arr[2];
    float blupper;
    float theta_parcel;
    float blmr;
    float tmp_env_theta;
    float tmp_env_dwpt;
    float tv_env;
    float tmp1;
    float tdef, tdef1, tdef2, tdef3;
    float lyre, lyrf, lyrlast;
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }

    // define the beginning pressure, temperature, and dewpoint
    // of the parcel based on the flag provided
    defineParcel(prof, pcl, idx, flag);

    int idx3D = P3(i, j, ktop, t, NX, NY, NZ);
    float pbot = pcl->pres;
    float ptop = pres_arr[idx3D];
    float tmpc = pcl->tmpc;
    float dwpc = pcl->dwpc;
    float pres = pbot;
    float cap_strength = MISSING;
    float cap_strengthpres = MISSING;
    float li_max = MISSING;
    float li_maxpres = MISSING;
    float totp = 0.f;
    float totn = 0.f;
    float tote = 0.f;
    float bplus = 0.f;
    float bminus = 0.f;
    float dp = -1.f;
    float mli;
    float mcap;






    // begin with the mixing layer
    pe1 = pbot;
    h1 = _interp1d(pres_arr, hght_arr, pbot, true, idx, NX, NY, NZ);
    tp1 = _virtemp(pres, tmpc, dwpc);

    // lift parcel and return LCL pres (hPa)
    // and LCL temp (C)
    _drylift(pres, tmpc, dwpc, dry_arr);
    pe2 = dry_arr[1];
    tp2 = dry_arr[0];

    // make sure the LCL pressure isn't below the surface
    if (pe2 > pres_arr[P3(i, j, kbot, t, NX, NY, NZ)]) {
        pcl->lclpres = pres_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    }
    else {
        pcl->lclpres = pe2;
    }
    // we want the lcl height above ground and not above sea level
    // so subract off the surface height
    pcl->lclhght = _interp1d(pres_arr, hght_arr, pcl->lclpres, true, idx, NX, NY, NZ) - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)];
    blupper = pe2;
    if (blupper < pres_arr[P3(i, j, ktop, t, NX, NY, NZ)]) {
        pcl->bplus = 0.f;
        pcl->bminus = 0.f;
        return;
    } 
    h2 = _interp1d(pres_arr, hght_arr, pe2, true, idx, NX, NY, NZ);
    te2 = _interp1d(pres_arr, vtmp_arr, pe2, true, idx, NX, NY, NZ);

    //Calculate lifted parcel theta for use in iterative CINH loop below
    //RECALL: lifted parcel theta is CONSTANT from LPL to LCL
    theta_parcel = _theta(pe2, tp2, 1000.f);

    //Environmental theta and mixing ratio at LPL
    blmr = _mixratio(pres, dwpc);

    // ACCUMULATED CINH IN THE MIXING LAYER BELOW THE LCL
    // this will be done in 'dp' increments and will use the virtual
    // temperature correction where possible
    int num_pp = (int)((pbot - (blupper + dp))/(-1 * dp));
    num_pp += 1;
    float hh[1000];
    float tdef_arr[1000];
    float pp_val;
    // initialize the array calle pp
    for (int ppidx = 0; ppidx < num_pp; ++ppidx) {
        // dp is -1
        pp_val = pbot + ppidx * dp;
        hh[ppidx] = _interp1d(pres_arr, hght_arr, pp_val, true, idx, NX, NY, NZ);

        tmp_env_theta = _theta(pp_val, _interp1d(pres_arr, tmpc_arr, pp_val, true, idx, NX, NY, NZ), 1000.);
        tmp_env_dwpt = _interp1d(pres_arr, dwpc_arr, pp_val, true, idx, NX, NY, NZ);
        tv_env = _virtemp(pp_val, tmp_env_theta, tmp_env_dwpt);
        tmp1 = _virtemp(pp_val, theta_parcel, _temp_at_mixrat(blmr, pp_val));
        tdef = (tmp1 - tv_env) / _ctok(tv_env);
        tdef_arr[ppidx] = tdef;
    }

    // calc the layer energy and accumulate the negative energy
    float lyre_loop = 0;
    int idx2;
    for (int idx1 = 0; idx1 < num_pp - 1; idx1++) {
        idx2 = idx1 + 1;
        lyre_loop = G * ( tdef_arr[idx1] + tdef_arr[idx2] ) / 2.f * ( hh[idx2] - hh[idx1] );
        if (lyre_loop < 0) {
            totn += lyre_loop;
        }
    }


    // move the bottom layer to the top of the boundary layer
    if (pbot > pe2) {
        pbot = pe2;
    }

    // check for the case where the LCL is above the upper
    // boundary of the data (e.g. a dropsonde)
    idx3D = P3(i, j, ktop, t, NX, NY, NZ);
    if (pbot < pres_arr[idx3D]) {
        pcl->bplus = bplus;
        pcl->bminus = bminus;
        return;
    }


    // find the lowest observation in layer
    int lptr = 999;
    int uptr = 0;
    for (int pidx = kbot; pidx < ktop+1; pidx++) {
        // starting at the surface, find the lowest observation 
        // where pbot > pres
        idx3D = P3(i, j, pidx, t, NX, NY, NZ);
        if ((pbot >= pres_arr[idx3D]) && (pidx < lptr)) lptr = pidx;
        // starting from the surface, find the highest observation
        // where pbot < pres
        if ((ptop <= pres_arr[idx3D]) && (pidx > uptr)) uptr = pidx;
    }

    // START WITH INTERPOLATED BOTTOM LAYER
    // begin moist ascent from lifted parcel LCL (pe2, tp2)
    pe1 = pbot;
    h1 = _interp1d(pres_arr, hght_arr, pe1, true, idx, NX, NY, NZ); 
    te1 = _interp1d(pres_arr, vtmp_arr, pe1, true, idx, NX, NY, NZ);
    tp1 = _wetlift(pe2, tp2, pe1);

    lyre = 0;


    for (int kz = lptr; kz < ktop+1; kz++) {
        idx3D = P3(i, j, kz, t, NX, NY, NZ);
        pe2 = pres_arr[idx3D];
        h2 = hght_arr[idx3D];
        te2 = vtmp_arr[idx3D];
        tp2 = _wetlift(pe1, tp1, pe2);
        tdef1 = (_virtemp(pe1, tp1, tp1) - te1) / _ctok(te1);
        tdef2 = (_virtemp(pe2, tp2, tp2) - te2) / _ctok(te2);
        lyrlast = lyre;
        lyre = G * (tdef1 + tdef2) / 2.f * (h2 - h1);

        if (lyre > 0) totp += lyre;
        else {
            if (pe2 > 500.f) totn += lyre;
        }
        /*
        // check for max li
        mli = _virtemp(pe2, tp2, tp2) - te2;
        if (mli > li_max) {
            li_max = mli;
            li_maxpres = pe2;
        }

        // check for max cap strength
        mcap = te2 - mli;
        if (mcap > cap_strength) {
            cap_strength = mcap;
            cap_strengthpres = pe2;
        }
        */

        tote += lyre;
        pelast = pe1;
        pe1 = pe2;
        te1 = te2;
        tp1 = tp2;

        // is this the top of the specified layer
        if (kz >= uptr) {
            pe3 = pe1;
            h3 = h1;
            te3 = te1;
            tp3 = tp1;
            lyrf = lyre;
            if (lyrf > 0.f) {
                bplus = totp - lyrf;
                bminus = totn;
            }
            else {
                bplus = totp;
                if (pe2 > 500.f) bminus = totn + lyrf;
                else bminus = totn;
            }
            pe2 = ptop;
            h2 = _interp1d(pres_arr, hght_arr, pe2, true, idx, NX, NY, NZ);
            te2 = _interp1d(pres_arr, vtmp_arr, pe2, true, idx, NX, NY, NZ);
            tp2 = _wetlift(pe3, tp3, pe2);
            tdef3 = (_virtemp(pe3, tp3, tp3) - te3) / _ctok(te3);
            tdef2 = (_virtemp(pe2, tp2, tp2) - te2) / _ctok(te2);
            lyrf = G * (tdef3 + tdef2) / 2.f * (h2 - h3);
            if (lyrf > 0.f) bplus += lyrf;
            else {
                if (pe2 > 500.f) bminus += lyrf;
            }
            if (bplus == 0.f) bminus = 0.f;
        }


        // is this the 3km level? If so, compute 0-3 CAPE
        if (pcl->lclhght < 3000.f) {
            if ((h1 - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)] <= 3000.f) && (h2 - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)] >= 3000.f)) {
                pe3 = pelast;
                h3 = _interp1d(pres_arr, hght_arr, pe3, true, idx, NX, NY, NZ);
                te3 = _interp1d(pres_arr, vtmp_arr, pe3, true, idx, NX, NY, NZ);
                tp3 = _wetlift(pe1, tp1, pe3);
                lyrf = lyre;
                if (lyrf > 0.f) {pcl->b3km = totp - lyrf;}
                else {pcl->b3km = totp;}
                h4 = _interp_to_msl(hght_arr, 3000.f, idx, NX, NY, NZ);
                pe4 = _interp1d(hght_arr, pres_arr, h4, false, idx, NX, NY, NZ);
                if ((pe2 != MISSING) && !(isnan(pe2))) {
                    te2 = _interp1d(pres_arr, vtmp_arr, pe4, true, idx, NX, NY, NZ);
                    tp2 = _wetlift(pe3, tp3, pe4);
                    tdef3 = (_virtemp(pe3, tp3, tp3) - te3) / (_ctok(te3));
                    tdef2 = (_virtemp(pe4, tp2, tp2) - te2) / (_ctok(te2));
                    lyrf = G * (tdef3 + tdef2) * 0.5f * (h4 - h3);
                    if (lyrf > 0.f) {pcl->b3km += lyrf;}
                }
            }
        }
        else pcl->b3km = 0.f;

        h1 = h2;

        
        // LFC possibility
        if ((lyre >= 0.f) && (lyrlast <=0.f)) {
            tp3 = tp1;
            pe2 = pe1;
            pe3 = pelast;

            if ( _interp1d(pres_arr, vtmp_arr, pe3, true, idx, NX, NY, NZ) < _virtemp(pe3, _wetlift(pe2, tp3, pe3), _wetlift(pe2, tp3, pe3))) {
                pcl->lfcpres = pe3;
                pcl->lfchght = _interp1d(pres_arr, hght_arr, pe3, true, idx, NX, NY, NZ) - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)];
            }
            else {
                while ((_interp1d(pres_arr, vtmp_arr, pe3, true, idx, NX, NY, NZ) > _virtemp(pe3, _wetlift(pe2, tp3, pe3), _wetlift(pe2, tp3, pe3))) && (pe3 > 0))
                {
                    pe3 -= 5.f;
                }

                if (pe3 > 0.f) {
                    pcl->lfcpres = pe3;
                    pcl->lfchght = _interp1d(pres_arr, hght_arr, pe3, true, idx, NX, NY, NZ) - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)];
                }
                else {
                    pcl->lfcpres = nanf("");
                    pcl->lfchght = nanf("");
                }
            }

            // hack to force the LFC to be at least the LCL
            if (pcl->lfcpres >= pcl->lclpres) {
                pcl->lfcpres = pcl->lclpres;
                pcl->lfchght = pcl->lclhght;
            }

        }

        // EL possibility
        if (lyre <= 0.f and lyrlast >= 0.f) {
            tp3 = tp1;
            pe2 = pe1;
            pe3 = pelast;

            while ( _interp1d(pres_arr, vtmp_arr, pe3, true, idx, NX, NY, NZ) < _virtemp(pe3, _wetlift(pe2, tp3, pe3), _wetlift(pe2, tp3, pe3)))
            {
                pe3 -= 5.f;
            }
            pcl->elpres = pe3;
            pcl->elhght = _interp1d(pres_arr, hght_arr, pe2, true, idx, NX, NY, NZ) - hght_arr[P3(i, j, kbot, t, NX, NY, NZ)];
            //pcl->limax = -1*li_max;
            //pcl->limaxpres = li_maxpres;

        }

        // 500 hPa lifted index
        if ((pres_arr[idx3D] <= 500.) && (pcl->li5 == MISSING)) {
            float a = _interp1d(pres_arr, vtmp_arr, 500., true, idx, NX, NY, NZ);
            float b = _wetlift(pe1, tp1, 500.);
            pcl->li5 = a - _virtemp(500., b, b);
        }

        // 300 hPa lifted index
        if ((pres_arr[idx3D] <= 300.) && (pcl->li3 != MISSING)) {
            float c = _interp1d(pres_arr, vtmp_arr, 300., true, idx, NX, NY, NX);
            float d = _wetlift(pe1, tp1, 300.);
            pcl->li3 = c - _virtemp(300., d, d);
        }

    }
    //out_bplus_bminus[0] = bplus;
    //out_bplus_bminus[1] = bminus;
    if (pcl->lfcpres == MISSING) pcl->lfcpres = nanf("");
    if (pcl->lfchght == MISSING) pcl->lfchght = nanf("");
    if (pcl->b3km == MISSING) pcl->b3km = 0.;
    pcl->bplus = bplus;
    pcl->bminus = bminus;
}


#ifdef __CUDACC__
__host__ __device__ 
#endif 
void _effective_inflow_layer(profile *prof, Parcel *mupcl, float ecape, float ecinh, int *idx) {
    float *pres_arr = prof->pres;
    float *tmpc_arr = prof->tmpc;
    float *dwpc_arr = prof->dwpc;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float mucape = mupcl->bplus; float mucinh = mupcl->bminus;
    float pbot = MISSING;
    float ptop = MISSING;
    int idx3D, bptr;
    Parcel loop_pcl;
    int kbot, ktop;
    // We need to find the valid data levels first instead
    // of assuming the surface is at k == 0
    kbot = 0;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }
    ktop = NZ-1;
    while (isnan(pres_arr[P3(i, j, ktop, t, NX, NY, NZ)])) {
        ktop -= 1;
    }

    if ((mucape != 0.f) && (mucape != MISSING)) {
        if ((mucape >= ecape) && (mucinh > ecinh)) {
            for (int iz = kbot; iz < ktop; ++iz) {
                idx3D = P3(i, j, iz, t, NX, NY, NZ);
                loop_pcl.pres = pres_arr[idx3D];
                loop_pcl.tmpc = tmpc_arr[idx3D];
                loop_pcl.dwpc = dwpc_arr[idx3D];
                loop_pcl.bplus = MISSING;
                loop_pcl.bminus = MISSING;
                _cape(prof, &loop_pcl, idx, 5);
                if ((loop_pcl.bplus >= ecape) && (loop_pcl.bminus > ecinh)) {
                    pbot = pres_arr[idx3D];
                    bptr = iz;
                    break;
                }
            }
            if (pbot == MISSING) {
                return;
            }
            for (int iz = bptr; iz < ktop; ++iz) {
                idx3D = P3(i, j, iz, t, NX, NY, NZ);
                loop_pcl.pres = pres_arr[idx3D];
                loop_pcl.tmpc = tmpc_arr[idx3D];
                loop_pcl.dwpc = dwpc_arr[idx3D];
                loop_pcl.bplus = MISSING;
                loop_pcl.bminus = MISSING;
                _cape(prof, &loop_pcl, idx, 5);
                if ((loop_pcl.bplus < ecape) || (loop_pcl.bminus < ecinh)) {
                    ptop = pres_arr[idx3D];
                    if(ptop > pbot) ptop = pbot;
                    break;
                }
            }

        }
        mupcl->effptop = ptop;
        mupcl->effpbot = pbot;
    }
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
void _bunkers_storm_motion(profile *prof, Parcel *mupcl, int *idx) {

    float *pres_arr = prof->pres;
    float *hght_arr = prof->hght;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    float d = 7.5f; // deviation value empirically derived at 7.5 m/s
    float mucape = mupcl->bplus;
    float muel = mupcl->elhght;
    float mn_wnd[2];
    float shr[2];
    float rstu, rstv;
    float lstu, lstv;


    float pbot = mupcl->effpbot;
    float base = _interp_to_agl(pres_arr, hght_arr, pbot, idx, NX, NY, NZ);
    if ((mucape > 100.f) && (mucape != MISSING) && (muel != MISSING) && (pbot != MISSING) && (!isnan(pbot))) {
        float depth = muel - base;
        float htop = base + (depth * 0.65f);
        float ptop = _interp1d(hght_arr, pres_arr, _interp_to_msl(hght_arr, htop, idx, NX, NY, NZ), false, idx, NX, NY, NZ);
        if ((ptop != MISSING) && (!isnan(ptop))) {
            _mean_wind(prof, pbot, ptop, mn_wnd, idx);
            _wind_shear(prof, pbot, ptop, shr, idx);
            float srmag = sqrt(pow(shr[0], 2.f) + pow(shr[1], 2.f));
            float uchg = d / srmag * shr[1];
            float vchg = d / srmag * shr[0];
            rstu = mn_wnd[0] + uchg;
            rstv = mn_wnd[1] - vchg;
            lstu = mn_wnd[0] - uchg;
            lstv = mn_wnd[1] + vchg;
        }
    }

    else {
        float npbm[4];
        _non_parcel_bunkers_motion(prof, npbm, idx);

        rstu = npbm[0];
        rstv = npbm[1];
        lstu = npbm[2];
        lstv = npbm[3];
    }

    mupcl->rstu = rstu;
    mupcl->rstv = rstv;
    mupcl->lstu = lstu;
    mupcl->lstv = lstv;

}

// lapse rate (c/km)
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _lapse_rate(profile *prof, float lower, float upper, bool isP, int *idx) {

    float *pres = prof->pres;
    float *hght = prof->hght;
    float *vtmp = prof->vtmp;
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    float p1, p2, z1, z2;

    if (isP) {
        p1 = lower;
        p2 = upper;
        z1 = _interp1d(pres, hght, lower, true, idx, NX, NY, NZ); 
        z2 = _interp1d(pres, hght, upper, true, idx, NX, NY, NZ); 
    }

    else {
        z1 = _interp_to_msl(hght, lower, idx, NX, NY, NZ);
        z2 = _interp_to_msl(hght, upper, idx, NX, NY, NZ);
        p1 = _interp1d(hght, pres, z1, false, idx, NX, NY, NZ);
        p2 = _interp1d(hght, pres, z2, false, idx, NX, NY, NZ);
    }

    float tv1 = _interp1d(pres, vtmp, p1, true, idx, NX, NY, NZ);
    float tv2 = _interp1d(pres, vtmp, p2, true, idx, NX, NY, NZ);

    return ((tv2 - tv1) / (z2 - z1)) * -1000.f;
    

}


#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _scp(float mucape, float srh, float ebwd) {

    if ((isnan(srh) )||(srh == MISSING)||(isnan(ebwd))) return 0.0;

    if (ebwd > 20.f) {
        ebwd = 20.f;
    }

    if (ebwd < 10.f) {
        ebwd = 10.f;
    }

    float muCAPE_term = mucape / 1000.f;
    float esrh_term = srh / 50.f;
    float ebwd_term = ebwd / 20.f;

    float scp = muCAPE_term * esrh_term * ebwd_term;
    return scp;

}


#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _stp_cin(float mlcape, float esrh, float ebwd, float mllcl, float mlcinh) {
    if ((isnan(esrh)) || (isnan(ebwd))) return 0.0;

    float cape_term = mlcape / 1500.f;
    float esrh_term = esrh / 150.f;
    float ebwd_term, lcl_term, cinh_term;

    if (ebwd < 12.5f) {
        ebwd_term = 0.f;
    }
    else if (ebwd > 30.f) {
        ebwd_term = 1.5f;
    }
    else {
        ebwd_term = ebwd / 20.f;
    }

    if (mllcl < 1000.f) {
        lcl_term = 1.0f;
    }
    else if (mllcl > 2000.f) {
        lcl_term = 0.f;
    }
    else {
        lcl_term = ((2000.f - mllcl) / 1000.f);
    }

    if (mlcinh > -50.f) {
        cinh_term = 1.0f;
    }

    else if (mlcinh < -200.f) {
        cinh_term = 0.f;
    }
    else {
        cinh_term = ((mlcinh + 200.f) / 150.f);
    }

    float stp_cin = cape_term * esrh_term * ebwd_term * lcl_term * cinh_term;
    if (stp_cin < 0.f) stp_cin = 0.f;
    return stp_cin;

}

/* Compute the Fixed Layer STP */
#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _stp_fixed(float sbcape, float sblcl, float srh01, float bwd6) {
    float stp_fixed;
    float lcl_term, bwd6_term, cape_term, srh_term;


    // compute the LCL term
    if (sblcl < 1000.) {
        lcl_term = 1.0;
    }
    else if (sblcl > 2000.) {
        lcl_term = 0.0;
    }
    else {
        lcl_term = ((2000.-sblcl)/1000.);
    }

    // compute the bulk wind difference term

    // greate than 30 m/s
    if (bwd6 > 30. ) { 
        bwd6 = 30.;
    }
    else if (bwd6 < 12.5) {
        bwd6 = 0.0;
    }

    bwd6_term = bwd6 / 20.;

    cape_term = sbcape / 1500.;
    srh_term = srh01 / 150.;

    stp_fixed = cape_term * lcl_term * srh_term * bwd6_term;
    return stp_fixed;
}


#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _vtp(float mlcape, float esrh, float ebwd, float mllcl, float mlcinh, float cape03, float lapse03, float effhbot) {

    if (effhbot > 0.f) return 0.f;
    if ((isnan(esrh))||(isnan(effhbot))) return 0.f;

    float cape_term = mlcape / 1500.f;
    float esrh_term = esrh / 150.f;
    float cape03_term = cape03 / 50.f;
    float lapse03_term = lapse03/6.5f;
    float ebwd_term, lcl_term, cinh_term;

    if (mllcl < 1000.f) {
        lcl_term = 1.0f;
    }
    else if (mllcl > 2000.f) {
        lcl_term = 0.f;
    }
    else {
        lcl_term = ((2000.f - mllcl) / 1000.f);
    }
    if (mlcinh > -50.f) {
        cinh_term = 1.0f;
    }
    else if (mlcinh < -200.f) {
        cinh_term = 0.f;
    }
    else {
        cinh_term = ((mlcinh + 200.f) / 150.f);
    }
    if (ebwd < 12.5f) {
        ebwd_term = 0.f;
    }
    else if (ebwd > 30.f) {
        ebwd_term = 1.5f;
    }
    else {
        ebwd_term = ebwd / 20.f;
    }
    float vtp = cape_term * esrh_term * ebwd_term * lcl_term * cinh_term * cape03_term * lapse03_term;
    return vtp;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _temp_lvl(profile *prof, float temp, int *idx, int NX, int NY, int NZ) {
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    float *tmpc_arr = prof->tmpc;
    float *pres_arr = prof->pres;
    float difft = 0.0;
    float difft_last = 0.0;
    bool found_t = false;
    float tval;

    int kbot = 0;
    int kidx;
    while (isnan(pres_arr[P3(i, j, kbot, t, NX, NY, NZ)])) {
        kbot += 1;
    }

    // find the first level (bottom-up) that we cross the given
    // temperature threshold
    for (kidx = kbot; kidx < NZ; ++kidx) {
        difft_last = difft; 
        tval = tmpc_arr[P3(i, j, kidx, t, NX, NY, NZ)];
        difft = temp - tval;
        if (difft == 0) {
            found_t = true;
            return pres_arr[P3(i, j, kidx, t, NX, NY, NZ)];
        }
        // have we reached an inflection point?
        else if (((difft < 0.0) && (difft_last > 0.0)) || ((difft > 0.0) && (difft_last < 0.0))) {
            found_t = true;
            break;
        }
    }

    if (found_t) {
        float interp_p[2];
        float interp_t[2];
        // the order of the x-coordinate (in this case, temperature) matters, as it needs
        // to be ascending. This should take care of that. 
        if (tmpc_arr[P3(i, j, kidx, t, NX, NY, NZ)] > tmpc_arr[P3(i, j, kidx-1, t, NX, NY, NZ)]) {
            interp_p[0] = log10(pres_arr[P3(i, j, kidx-1, t, NX, NY, NZ)]);
            interp_p[1] = log10(pres_arr[P3(i, j, kidx, t, NX, NY, NZ)]);
            interp_t[0] = tmpc_arr[P3(i, j, kidx-1, t, NX, NY, NZ)];
            interp_t[1] = tmpc_arr[P3(i, j, kidx, t, NX, NY, NZ)];

        }
        else {
            interp_p[0] = log10(pres_arr[P3(i, j, kidx, t, NX, NY, NZ)]);
            interp_p[1] = log10(pres_arr[P3(i, j, kidx-1, t, NX, NY, NZ)]);
            interp_t[0] = tmpc_arr[P3(i, j, kidx, t, NX, NY, NZ)];
            interp_t[1] = tmpc_arr[P3(i, j, kidx-1, t, NX, NY, NZ)];
        }
              
        // linearlly interpolate it to the best pressure value
        float pval = interp_p[0] + (temp - interp_t[0]) * ((interp_p[1] - interp_p[0])/(interp_t[1] - interp_t[0]));
        return powf(10., pval);
    }
    else {
        return nanf("");
    }
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
float _ship(profile *prof,Parcel *mupcl, float *shr06, int *idx, int NX, int NY, int NZ) {
    float mumr = _mixratio(mupcl->pres, mupcl->dwpc);
    float frz_lvl_prs = _temp_lvl(prof, 0.f, idx, NX, NY, NZ);
    float frz_lvl = _interp1d(prof->pres, prof->hght, frz_lvl_prs, true, idx, NX, NY, NZ);
    float h5_temp = _interp1d(prof->pres, prof->tmpc, 500., true, idx, NX, NY, NZ);
    float lr75 = _lapse_rate(prof, 700., 500., true, idx);
    float shr = sqrt(powf(shr06[0],2.) + powf(shr06[1], 2.)); 
    float ship = 0.f;

    if (shr > 27) shr = 27.;
    else if (shr < 7) shr = 7.;

    if (mumr > 13.6) mumr = 13.6;
    else if (mumr < 11.) mumr = 11.;

    if (h5_temp > -5.5) h5_temp = -5.5;

    ship = -1. * (mupcl->bplus * mumr * lr75 * h5_temp * shr) / 42000000.;
   
    if (mupcl->bplus < 1300) ship = ship*(mupcl->bplus/1300.);
    
    if (lr75 < 5.8) ship = ship*(lr75/5.8);

    if (frz_lvl < 2400) ship = ship * (frz_lvl/2400.);
    
    return ship;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
float c_totals(profile *prof, int *idx) {
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    float dwp_850 = _interp1d(prof->pres, prof->dwpc, 850., true, idx, NX, NY, NZ); 
    float tmp_500 = _interp1d(prof->pres, prof->tmpc, 500., true, idx, NX, NY, NZ);
    return dwp_850 - tmp_500;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
float v_totals(profile *prof, int *idx) {
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    float tmp_850 = _interp1d(prof->pres, prof->tmpc, 850., true, idx, NX, NY, NZ); 
    float tmp_500 = _interp1d(prof->pres, prof->tmpc, 500., true, idx, NX, NY, NZ);
    return tmp_850 - tmp_500;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
float t_totals(profile *prof, int *idx) {
    return c_totals(prof, idx) + v_totals(prof, idx); 
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
float k_index(profile *prof, int *idx) {
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;

    float tmp_850 = _interp1d(prof->pres, prof->tmpc, 850., true, idx, NX, NY, NZ); 
    float tmp_700 = _interp1d(prof->pres, prof->tmpc, 700., true, idx, NX, NY, NZ); 
    float tmp_500 = _interp1d(prof->pres, prof->tmpc, 500., true, idx, NX, NY, NZ); 
    float dwp_850 = _interp1d(prof->pres, prof->dwpc, 850., true, idx, NX, NY, NZ); 
    float dwp_700 = _interp1d(prof->pres, prof->dwpc, 700., true, idx, NX, NY, NZ); 

    return tmp_850 - tmp_500 + dwp_850 - (tmp_700 - dwp_700);
}


#ifdef __CUDACC__
__host__ __device__ 
#endif 
float sweat(profile *prof, int *idx) {
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    float vec_850[2];
    float vec_500[2];
    float uv_850[2];
    float uv_500[2];
    float term1, term2, term3, term4, term5;
    const float PI = 3.14159265;

    float dwp_850 = _interp1d(prof->pres, prof->dwpc, 850., true, idx, NX, NY, NZ); 
    uv_850[0] = _interp1d(prof->pres, prof->uwin, 850., true, idx, NX, NY, NZ) * 1.94384; // need to convert from m/s to kts
    uv_850[1] = _interp1d(prof->pres, prof->vwin, 850., true, idx, NX, NY, NZ) * 1.94384;
    uv_500[0] = _interp1d(prof->pres, prof->uwin, 500., true, idx, NX, NY, NZ) * 1.94384;
    uv_500[1] = _interp1d(prof->pres, prof->vwin, 500., true, idx, NX, NY, NZ) * 1.94384;
    comp2vec(uv_850, vec_850);
    comp2vec(uv_500, vec_500);
    float tt = t_totals(prof, idx);

    if (dwp_850 < 0.) term1 = 12. * dwp_850;
    else term1 = 0;


    if (tt < 49) term2 = 0;
    else term2 = 20. * (tt - 49);
    term3 = 2. * vec_850[0];
    term4 = vec_500[0];

    if ((130 <= vec_850[1]) && (250. >= vec_850[1]) && (210 <= vec_500[1]) && \
       (310 >= vec_500[1]) && (vec_500[1] > 0.) && (vec_850[0] >= 15.) && (vec_500[0] >= 15.)) {
        term5 = 125 * ( sin( (vec_500[1] - vec_850[1]) * (PI / 180.) ) + 0.2);
    }
    else term5 = 0.;
    free(uv_850);
    free(uv_500);
    free(vec_850);
    free(vec_500);

    return term1 + term2 + term3 + term4 + term5;
}

#ifdef __CUDACC__
__host__ __device__ 
#endif 
void bulk_rich(profile *prof, Parcel *pcl, int *idx) {
    int i = idx[0]; int j = idx[1]; int t = idx[3];
    int NX = prof->NX;
    int NY = prof->NY;
    int NZ = prof->NZ;
    int kbot = 0;
    while ((isnan(prof->pres[P3(i, j, kbot, t, NX, NY, NZ)])) && (kbot < NZ) ) {
        kbot += 1;
    }
    if (kbot == NZ) {
        pcl->brnu = nanf("");
        pcl->brnv = nanf("");
        pcl->brn = nanf("");
        pcl->brnshear = nanf(""); 
    }

    float pbot = prof->pres[P3(i, j, kbot, t, NX, NY, NZ)];
    float ptop = _interp1d(prof->hght, prof->pres, _interp_to_msl(prof->hght, 6000., idx, NX, NY, NZ), false, idx, NX, NY, NZ);
    float p500m = _interp1d(prof->hght, prof->pres,  _interp_to_msl(prof->hght, 500.f, idx, NX, NY, NZ), false, idx, NX, NY, NZ);


    // SFC-500m mean wind
    float mn_wind_500m[2];
    float mn_wind_6km[2];
    _mean_wind(prof, pbot, p500m, mn_wind_500m, idx);
    _mean_wind(prof, pbot, ptop, mn_wind_6km, idx);

    if ((isnan(pcl->bplus)) || (isnan(mn_wind_500m[0])) || (isnan(mn_wind_500m[1])) || (isnan(mn_wind_6km[0])) || (isnan(mn_wind_6km[1]))) {
        pcl->brnu = nanf("");
        pcl->brnv = nanf("");
        pcl->brn = nanf("");
        pcl->brnshear = nanf(""); 
    }

    float dx = mn_wind_6km[0] - mn_wind_500m[0];
    float dy = mn_wind_6km[1] - mn_wind_500m[1];
    pcl->brnu = dx;
    pcl->brnv = dy;
    pcl->brnshear = (powf(dx, 2.) + powf(dy, 2.)) / 2.;
    pcl->brn = pcl->bplus / pcl->brnshear;
}
