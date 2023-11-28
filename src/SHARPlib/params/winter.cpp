// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2023-11-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC.
#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params/winter.h>
#include <SHARPlib/thermo.h>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace sharp {

void posneg_temperature(float start, const float pres_arr[], const float height_arr[],
		const float temperature_arr[], const int N, float* pos, float* neg) noexcept {
        float upper, lower, pe1, h1, te1, tp1, totp, totn, pe2, h2, te2,
              tp2, tdef1, tdef2;
        float lyrlast, lyre, tote, pelast, ptop, pbot, lvl;
        short i, lptr, uptr, warmlayer=0, coldlayer=0, phase;

        // If there is no sounding or if there is no near-saturated layer
        // in the lowest 5 km, do not compute
        if (((interp_pressure(500,pres_arr,temperature_arr,N)<-273.15) && (interp_pressure(850,pres_arr,temperature_arr,N)<-273.15))||(start<-1)) return;

        /* ----- Get lowest and highest observations in layer ----- */
        lptr = 0;
        if (start=-1) {
        	upper = 500.0;
        }
        else {
          	upper = start;
        }

        i=N-1;
        while(pres_arr[i] < upper) {
          i--;
          if (i < 0) {
            // fprintf(stderr,"Warning: posneg_temp: Could not find a pressure greater than %.2f\n",upper);
            // fprintf(stderr, "Using %.2f as the upper level.\n",pres_arr[0]);
            i = 0;
            break;
          }
        }
        uptr = i;
        if (pres_arr[i] == upper) uptr--;

        /* ----- Start with top layer ----- */
        pe1 = upper;
        h1 =  interp_pressure(pe1,pres_arr,height_arr,N);
        te1 = interp_pressure(pe1,pres_arr,temperature_arr,N);
        tp1 = 0;

        totp = totn = tote = ptop = pbot = 0;

        for( i = uptr; i >= lptr; i--) {
		           if (temperature_arr[i]>-273.15) {
              /* ----- Calculate every level that reports a temp ----- */
              pe2 = pres_arr[i];
              h2  =  height_arr[i];
              te2 = temperature_arr[i];
              tp2 = 0;
              tdef1 = (0 - te1) / (te1 + 273.15);
              tdef2 = (0 - te2) / (te2 + 273.15);
              lyrlast = lyre;
              lyre = 4.9 * (tdef1 + tdef2) * (h2 - h1);

              /* Has a warm layer been found yet? */
              if (te2>0)
                if (warmlayer==0) {
                  warmlayer=1;
                  ptop=pe2;
                }

              /* Has a cold layer been found yet? */
              if (te2<0)
                if ((warmlayer==1) && (coldlayer==0)) {
                  coldlayer=1;
                  pbot=pe2;
                }

              if (warmlayer>0) {
                if (lyre>0)
                  totp += lyre;
                else
                  totn += lyre;

                tote += lyre;
                // printf("%4.0f - %4.0f E=%6.0f TOT=%6.0f Top=%6.0f Bot=%6.0f\n",
                //  pe1, pe2, lyre, tote, ptop, pbot);
              }

              pelast = pe1;
              pe1 = pe2;
              h1  = h2;
              te1 = te2;
              tp1 = tp2;
           }
        }
        if ((warmlayer==1) && (coldlayer==1)) {
          *pos = totp;
          *neg = totn;
          // printf("Tot= %.0f J/kg   Pos= %.0f J/kg   Neg= %.0f J/kg\n",
          //  tote, totp, totn);
          // printf("Top= %.0f        Bot= %.0f\n", ptop, pbot);
        }
        else {
          // printf ("Warm/Cold Layers not found.\n" );
          *pos = 0;
          *neg = 0;
        }
}

char *init_phase(const float pres_arr[], const float height_arr[], const float temperature_arr[],
		const float dewpoint_arr[], const float vvel_arr[], const int N, float *plevel, short *phase) noexcept {
        short       i, ok, avail;
        float       rh, p1, pbegin, pos, neg, p, w1;
        float       p_pres, p_phase, ptop, pbot;
        static char st[80];
        static char pt[80];

        *plevel = 0;
        *phase = -1;

        /* First, determine whether VVELS are available.  If they are,   */
        /* use them to determine level where precipitation will develop. */
        avail=0;
        for( i = 0; i < N; i++) {
            if (vvel_arr[i] > 0) avail++;
        }

        if (avail< 5) {
           /* No VVELS...must look for saturated level */

           /* ----- Find the highest near-saturated 50mb
 *               layer blo 5km agl ---- */
           for(i=N-1;i>0;i--) {
              ok = 0;
              pbegin = -999;
              if (height_arr[i]-height_arr[0] < 5000.0) {
                  rh = relh(pres_arr[i],temperature_arr[i],dewpoint_arr[i]); 
                  if (rh > 80) {
                    p1 = pres_arr[i]+5000;
                    if (relh(p1,interp_pressure(p1,pres_arr,temperature_arr,N),interp_pressure(p1,pres_arr,dewpoint_arr,N)) > 80) {
                                ok = 1;
                                pbegin = p1-2500.0;
                                break;
                    }
                }
              }
           }
        }
        else {
                /* ----- Find the highest near-saturated layer with UVV in 
 *                    the lowest 5km ----- */
                for(i=N-1;i>0;i--) {
                        ok=0;
                        pbegin=-999;
                        if ((height_arr[i]-height_arr[0] < 5000.0) && (vvel_arr[i] >= 0)) {
                  		rh = relh(pres_arr[i],temperature_arr[i],dewpoint_arr[i]); 
                                if (rh > 80) {
                                        p1 = pres_arr[i]+5000;
                    			if (relh(p1,interp_pressure(p1,pres_arr,temperature_arr,N),interp_pressure(p1,pres_arr,dewpoint_arr,N)) > 80) {
                                                ok = 1;
                                                pbegin = p1-2500;
                                                break;
                                	}
                                }
                        }
                }
        }

        if (!ok) {
           *plevel =  MISSING;
           *phase  = -1;
           strcpy(st, "N/A");
           return st;
        }

        p1 = interp_pressure(pbegin,pres_arr,temperature_arr,N);
        if(p1>0) {p_phase=0; strcpy(pt, "Rain"); }
        if((p1<=0) && (p1 > -5)) {p_phase=1; strcpy(pt, "Freezing Rain"); }
        if((p1<=-5) && (p1 > -9)) {p_phase=1; strcpy(pt, "ZR/S Mix" ); }
        if(p1 <= -9) {p_phase=3; strcpy(pt, "Snow"); }
        *plevel = pbegin;
        *phase = p_phase;
        return pt;
}

}  // end namespace sharp

