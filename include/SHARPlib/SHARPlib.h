/**
 * \file
 * \brief Wrapper header for the entire SHARP library
 *
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-12-01
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#ifndef __SHARPLIB_H__
#define __SHARPLIB_H__

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/layer.h>
#include <SHARPlib/params.h>
#include <SHARPlib/parcel.h>
#include <SHARPlib/profile.h>
#include <SHARPlib/thermo.h>
#include <SHARPlib/winds.h>

namespace sharp {

Profile *create_profile(const float pres[], const float hght[],
                        const float tmpk[], const float dwpk[],
                        const float wspd_or_u[], const float wdir_or_v[],
                        int NZ, Source sounding_type, bool windComponents) {
    Profile *prof = new Profile(NZ, sounding_type);

    for (int k = 0; k < NZ; k++) {
        float p = pres[k];
        float h = hght[k];
        float t = tmpk[k];
        float d = dwpk[k];

        float vtmp = virtual_temperature(p, t, d);
        float mixr = mixratio(p, d);
        float thta = theta(p, t, sharp::THETA_REF_PRESSURE);
        float thte = thetae(p, t, d);

        prof->pres[k] = p; 
        prof->hght[k] = h; 
        prof->tmpk[k] = t;
        prof->dwpk[k] = d; 

        prof->vtmp[k] = vtmp;
        prof->mixr[k] = mixr;
        prof->theta[k] = thta;
        prof->theta_e[k] = thte;
        prof->buoyancy[k] = 0.0f;

        float specific_humidity = (1.0f - mixr) * (mixr);
        if (mixr == MISSING) specific_humidity = MISSING;

        float height_agl = prof->hght[k] - prof->hght[0];
        prof->moist_static_energy[k] =
            moist_static_energy(height_agl, prof->tmpk[k], specific_humidity);

        if (windComponents) {
            // converting from knots to m/s
            prof->uwin[k] = wspd_or_u[k];
            prof->vwin[k] = wdir_or_v[k];
            if (prof->uwin[k] != MISSING) prof->uwin[k] *= 0.514444f;
            if (prof->vwin[k] != MISSING) prof->vwin[k] *= 0.514444f;

            WindVector vec = components_to_vector(prof->uwin[k], prof->vwin[k]);
            prof->wspd[k] = vec.speed;
            prof->wdir[k] = vec.direction;
        } else {
            prof->wspd[k] = wspd_or_u[k];
            if (prof->wspd[k] != MISSING) prof->wspd[k] *= 0.514444f;
            prof->wdir[k] = wdir_or_v[k];

            WindComponents cmp =
                vector_to_components(prof->wspd[k], prof->wdir[k]);
            prof->uwin[k] = cmp.u;
            prof->vwin[k] = cmp.v;
        }
    }

    return prof;
}

}  // end namespace sharp

#endif // __SHARPLIB_H__
