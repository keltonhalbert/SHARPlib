// Author: Kelton Halbert
// Email: kelton.halbert@noaa.gov
// License: Apache 2.0
// Date: 2022-10-13
//
// Written for the NWS Storm Predidiction Center
// Based on NSHARP routines originally written by
// John Hart and Rich Thompson at SPC. 
#pragma once

namespace sharp {
    

    struct profile {

        float *pres;
        float *hght;
        float *tmpc;
        float *dwpc;
        float *mixr;
        float *vtmp;
        float *uwin;
        float *vwin;

        int nlevs;
        int type; // used to determine what kind of profile (e.g. observed, model)

    };
}
