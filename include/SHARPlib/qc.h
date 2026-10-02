/**
 * \file
 * \brief Routines for detecting missing data
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 * \date   2026-10-01
 *
 * Written for the NWS Storm Prediction Center \n
 */

#ifndef SHARP_QC_H
#define SHARP_QC_H

#include <SHARPlib/constants.h>

#include <cmath>

namespace sharp {

/**
 * \author Kelton Halbert - NWS Storm Prediction Center
 *
 * \brief Check whether a value is sharp::MISSING or NaN
 *
 * Returns true if value is sharp::MISSING or NaN.
 *
 * \param   value   The value to check
 *
 * \return  Whether value is sharp::MISSING or NaN
 */
[[nodiscard]] inline bool is_missing(const float value) {
    return !std::islessgreater(value, MISSING);
}

}  // namespace sharp

#endif  // SHARP_QC_H
