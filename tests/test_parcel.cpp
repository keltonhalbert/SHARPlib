#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "constants.h"
#include "profile.h"
#include "parcel.h"

TEST_CASE("Testing parcel definitions") {

    sharp::Profile prof(10, sharp::Source::Observed);
    sharp::Parcel pcl;

    sharp::lifter_wobus lifter;
    sharp::define_parcel(&prof, &pcl, sharp::LPL::SFC);

    sharp::lift_parcel<sharp::lifter_wobus>(lifter, &prof, &pcl);

}
