/**
 * \file
 * \brief Quality control: how SHARPlib recognizes missing data
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2026-10-01
 *
 * Written for the NWS Storm Predidiction Center \n
 */

#ifndef SHARP_QC_H
#define SHARP_QC_H

#include <SHARPlib/constants.h>

#include <cmath>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Whether a value is missing: sharp::MISSING or NaN
 *
 * NaN should never reach the library, but where it does, the quality-control
 * paths (builds without NO_QC) treat it the same as sharp::MISSING. Routines
 * that skip or bridge missing data use this check, so that every one of them
 * recognizes missing data the same way.
 *
 * \param   value   The value to check
 *
 * \return  Whether value is sharp::MISSING or NaN
 */
[[nodiscard]] inline bool is_missing(const float value) {
    return (value == MISSING) || std::isnan(value);
}

}  // namespace sharp

#endif  // SHARP_QC_H
