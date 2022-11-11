#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "constants.h"
#include "profile.h"
#include "parcel.h"

TEST_CASE("Testing parcel definitions") {

    sharp::Profile prof(10, sharp::Source::Observed);

    // need to fill with dummy data to prevent
    // segmentation fault
    for (int i = 0; i < 10; i++) {
        prof.pres[i] = 1000.0 - 100.0*i;
        prof.hght[i] = 0 + 1000.0*i;
        prof.tmpc[i] = i;
        prof.dwpc[i] = i;

    }

    sharp::Parcel pcl;

    sharp::lifter_wobus lifter;
    sharp::define_parcel(&prof, &pcl, sharp::LPL::SFC);

    sharp::lift_parcel<sharp::lifter_wobus>(lifter, &prof, &pcl);

}
