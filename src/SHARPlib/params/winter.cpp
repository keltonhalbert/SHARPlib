/**
 * \file
 * \brief Winter weather parameters
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2023-05-30
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#include <SHARPlib/constants.h>
#include <SHARPlib/interp.h>
#include <SHARPlib/params/winter.h>

#include <cstddef>

namespace sharp {

PressureLayer dendritic_layer(const float pressure[], const float temperature[],
                              const std::ptrdiff_t N) {
    sharp::PressureLayer dgz = {sharp::MISSING, sharp::MISSING};

    dgz.bottom =
        find_first_pressure(-12.0f + ZEROCNK, pressure, temperature, N);
    dgz.top = find_first_pressure(-17.0f, pressure, temperature, N);

    // do_stuff
    return dgz;
};

}  // namespace sharp
