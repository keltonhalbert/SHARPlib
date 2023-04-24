/**
 * \file
 * \brief Various algorithms used for searching, sorting, accumulating etc 
 * \author  
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2023-04-12
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC. 
 */

#ifndef __SHARP_ALGORITHMS
#define __SHARP_ALGORITHMS

#include <functional>

namespace sharp {


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Finds the index of the lower bound that does not satisfy element < value 
 *
 * Based on std::lower_bound, this iterates over an array using a binary search
 * to find the first element that does not satisfy the comparison condition. By
 * default, the comparison is std::less. 
 *
 * We use a custom implementation of sharp::lower_bound rather than 
 * std::lower_bound for a few reasons. First, we prefer raw pointer 
 * arrays over vectors for easy integration with SWIG bindings, 
 * C code, and the potential future in which this runs on CUDA 
 * architecture. Currently, the algorithms in the STL library are 
 * not supported by the CUDA STL, but the types in <functional> 
 * (i.e. std::less) are supported by the CUDA STL. 
 *
 * \param array     The array to search over 
 * \param N         The length of the array
 * \param value     The value for lower-bound comparison
 * \param cmp       The comparitor 
 */
template <typename T, typename C = std::less<>>
inline constexpr int lower_bound(const T* array, int N, T value, C cmp=C{}) noexcept {
    
    int idx = 0;
    int first = 0;
    int count = N-1;

    while (count > 0) {
        idx = first;
        int step = count / 2;
        idx += step;
        T element = array[idx];

        if (cmp(element, value)) {
            first = ++idx;
            count -= step + 1;
        }
        else {
            count = step;
        }
    }

    return first;
}


/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 * 
 * \brief Finds the first index that satisfies value < element 
 *
 * Based on std::upper_bound, this iterates over an array using a binary search
 * to find the first element that satisfies the comparison condition. By
 * default, the comparison is std::less. 
 *
 * We use a custom implementation of sharp::upper_bound rather than 
 * std::upper_bound for a few reasons. First, we prefer raw pointer 
 * arrays over vectors for easy integration with SWIG bindings, 
 * C code, and the potential future in which this runs on CUDA 
 * architecture. Currently, the algorithms in the STL library are 
 * not supported by the CUDA STL, but the types in <functional> 
 * (i.e. std::less) are supported by the CUDA STL. 
 *
 * \param array     The array to search over 
 * \param N         The length of the array
 * \param value     The value for upper-bound comparison
 * \param cmp       The comparitor 
 */
template <typename T, typename C = std::less<>>
inline constexpr int upper_bound(const T* array, int N, T value, C cmp=C{}) noexcept {

    int idx = 0;
    int first = 0;
    int count = N-1;

    while (count > 0) {
        idx = first;
        int step = count / 2;
        idx += step;
        T element = array[idx];

        if (!cmp(value, element)) {
            first = ++idx;
            count -= step + 1;
        }
        else {
            count = step;
        }
    }

    return first;
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief A kernel that computes the area under a curve using the trapezoidal rule.
 *
 * This is a function kernel that integrates an area using the trapezoidal rule. 
 * This function does not do the actual integration over an array, but is meant
 * to be wrapped by any function that does so. 
 *
 * \param var_top      The top value of the variable to integrate
 * \param var_bottom   The bottom value of the variable to integrate
 * \param coord_top    The top value of the coordinate to integrate
 * \param coord_bottom The bottom value of the coordinate to integrate
 * \param weights      The weights to accumulate
 * \param weighted     Whether or not to accumulate weights
 * \return             The area under the curve
 */
template <typename _T>
inline constexpr _T __integ_trapz(_T var_top, _T var_bottom, 
		                       _T coord_top, _T coord_bottom, 
							   _T& weights, bool weighted=false) noexcept {
	if (weighted) {
		weights += coord_top - coord_bottom;
	}
	return ((var_top + var_bottom) / 2.0) * (coord_top - coord_bottom);
}


} // end namespace sharp

#endif // __SHARP_ALGORITHMS
