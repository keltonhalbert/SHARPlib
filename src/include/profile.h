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
