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
#include <iostream>

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
int lower_bound(const T* array, int N, T value, C cmp=C{}) {
    
    int idx = 0;
    int first = 0;
    int count = N;

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
int upper_bound(const T* array, int N, T value, C cmp=C{}) {

    int idx = 0;
    int first = 0;
    int count = N;

    // TO-DO:
    // I have an overflow here somewhere... need to 
    // figure that out -_-
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



} // end namespace sharp

#endif // __SHARP_ALGORITHMS
