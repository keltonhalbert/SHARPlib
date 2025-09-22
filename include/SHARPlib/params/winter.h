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

#ifndef SHARP_PARAMS_WINTER_H
#define SHARP_PARAMS_WINTER_H

#include <SHARPlib/constants.h>
#include <SHARPlib/layer.h>

#include <cstddef>

namespace sharp {
/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Get the layer encompassing the lowest dendritic growth zone (-12 to
 * -17 C)
 *
 * Search for and return the sharp::PressuerLayer of the lowest altitude
 * dendritic growth zone. If none is found, the top and bottom pressure levels
 * are set to sharp::MISSING.
 *
 * \param   pressure    (Pa)
 * \param   temperature (K)
 * \param   N           (length of arrays)
 *
 * \return   The top and bottom of the dendritic growth zone (Pa)
 */
[[nodiscard]] PressureLayer dendritic_layer(const float pressure[],
                                            const float temperature[],
                                            const std::ptrdiff_t N);
}  // namespace sharp

#endif  // SHARP_PARAMS_WINTER_H
