/**
 * \file
 * \brief Routines for linear interpolation of vertical atmospheric profiles
 * \author
 *   Kelton Halbert                  \n
 *   Email: kelton.halbert@noaa.gov  \n
 *   License: Apache 2.0             \n
 * \date   2022-10-13
 *
 * Written for the NWS Storm Predidiction Center \n
 * Based on NSHARP routines originally written by
 * John Hart and Rich Thompson at SPC.
 */

#ifndef __SHARP_INTERP_H__
#define __SHARP_INTERP_H__

#include <SHARPlib/constants.h>

#include <cmath>

namespace sharp {

/**
 * \author C++20 Standard Template Library <cmath>
 *
 * \brief Linearly interpolates between a and b with distance t.
 *
 * This routine was copied from the C++20 <cmath> standard template library.
 * It is used to linearly interpolate between A and B over a normalized
 * distance T, where 0 <= T <= 1. Full documentation and details can be found
 * in the paper located here: 
 * https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html
 *
 * This linear interpolation guarantees the following:
 *
 * 1. Exactness: lerp(a,b,0)==a & lerp(a,b,1)==b
 * 2. Monotonicity: cmp(lerp(a,b,t2), lerp(a,b,t1)) * cmp(t2, t1) <!-- 
 * -->* cmp(b,a) >= 0, where cmp is an arithmetic three-way comparison function
 * 3. Determinacy: result of NaN only for lerp(a,a,INFINITY)
 * 4. Boundedness: t<0 || t>1 || isfinite(lerp(a,b,t))
 * 5. Consistency: lerp(a,a,t) == a
 *
 * \param   __a     left value
 * \param   __b     right value
 * \param   __t     where 0 <= __t <= 1.0
 *
 * \return  The value between __a and __b at distance __t between them.
 */
template <typename _Fp>
[[nodiscard]] constexpr _Fp __lerp(_Fp __a, _Fp __b, _Fp __t)  {
    if ((__a <= 0 && __b >= 0) || (__a >= 0 && __b <= 0))
        return __t * __b + (1 - __t) * __a;

    if (__t == 1) return __b;  // exact

    // Exact at __t=0, monotonic except near __t=1,
    // bounded, determinate, and consistent:
    const _Fp __x = __a + __t * (__b - __a);
    return (__t > 1) == (__b > __a)
               ? (__b < __x ? __x : __b)
               : (__b > __x ? __x : __b);  // monotonic near __t=1
}

/**
 * \author C++20 Standard Template Library <cmath>
 *
 * \brief Linearly interpolates between a and b with distance t.
 *
 * This routine was copied from the C++20 <cmath> standard template library.
 * It is used to linearly interpolate between A and B over a normalized
 * distance T, where 0 <= T <= 1. Full documentation and details can be found
 * in the paper located here:
 https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html

 *
 * This linear interpolation guarantees the following:
 *
 * 1. Exactness: lerp(a,b,0)==a & lerp(a,b,1)==b
 * 2. Monotonicity: cmp(lerp(a,b,t2), lerp(a,b,t1)) * cmp(t2, t1) <!-- 
 * -->* cmp(b,a) >= 0, where cmp is an arithmetic three-way comparison function
 * 3. Determinacy: result of NaN only for lerp(a,a,INFINITY)
 * 4. Boundedness: t<0 || t>1 || isfinite(lerp(a,b,t))
 * 5. Consistency: lerp(a,a,t) == a
 *
 * \param   __a     left value
 * \param   __b     right value
 * \param   __t     where 0 <= __t <= 1.0
 *
 * \return  The value between __a and __b at distance __t between them.
 */
[[nodiscard]] constexpr float lerp(float __a, float __b, float __t)  {
    return __lerp(__a, __b, __t);
}

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Linearly interpolate a data field to a requested height level.
 *
 * Interpolates data in height coordinates (meters) when given the
 * array of heights, the array of data to interpolate, a height
 * level to interpolate to, and the number of levels in the arrays.
 *
 * The height array must be sorted in ascending order and monotonic.
 * For performance reasons, this routine assumes your data is well
 * ordered. If weird stuff happens, check the ordering of your array
 * values first.Duplicate height values or decreasing height values 
 * may produce unexpected results.
 *
 * \param   height_val  The height value to interpolate data to
 * \param   height_arr  The array of height values in the profile
 * \param   data_arr    The array of data values you want interpolated
 * \param   N           The the array length (number of elements)
 *
 * \return  The value of data_arr at the requested height_val.
 */
[[nodiscard]] float interp_height(float height_val, const float height_arr[],
                                  const float data_arr[], const int N) ;

/**
 * \author Kelton Halbert - NWS Storm Prediction Center/OU-CIWRO
 *
 * \brief Linearly interpolate a data field to a requested pressure level.
 *
 * Interpolates data in pressure coordinates (mb) when given the
 * array of pressure values, the array of data to interpolate, a
 * pressure level to interpolate to, and the number of levels in
 * the arrays.
 *
 * The pressure array must be sorted in descending order and monotonic.
 * Duplicate pressure values or increasing pressure values may produce
 * unexpected results.
 *
 * \param   pressure_val    The pressure value to interpolate data to
 * \param   pressure_arr    The array of pressure values in the profile
 * \param   data_arr        The array of data values you want interpolated
 * \param   N               The the array length (number of elements)
 *
 * \return  The value of data_arr at the requested pressure_val.
 */
[[nodiscard]] float interp_pressure(float pressure_val,
                                    const float pressure_arr[],
                                    const float data_arr[],
                                    const int N);
/**
 *  \author Nathan Dahl - NWS Storm Prediction Center/OU-CIWRO (mutated from Kelton's code)
 * 
 *  \brief Use linear interpolation to obtain the lowest height
 *  above the surface where a given variable equals a desired value
 *  (e.g., the height of the 0 C isotherm).
 *  
 *  The height array must be in ascending order and monotonic. Note that the
 *  given variable is assumed to generally decrease with height; if the surface value
 *  is less than te desired value, the code will return a "missing" value for the height.
 *  
 *  param   data_val        The variable value corresponding to the desired height
 *  param   data_arr        Array containing the vertical profile of the variable
 *  param   height_arr      Array containing the height levels of the profile
 *  param   N               The profile array length (number of elements)
 *  
 *  \return  The value of height_arr at the requested data_val.
*/
[[nodiscard]] float interp_hghtlevel(float data_val,
                                    const float data_arr[],
                                    const float height_arr[],
                                    const int N) noexcept;

/**
 *  \author Nathan Dahl - NWS Storm Prediction Center/OU-CIWRO
 *
 *  \brief Search for lowest altitude at which wetbulb temperature is zero, using
 *  linear interpolation between levels
 *
 *  param   p_arr           Array containing the pressure levels of the profile
 *  param   t_arr	    Array containing the vertical temperature profile
 *  param   td_arr          Array containing the vertical dewpoint profile
 *  param   height_arr      Array containing the height levels of the profile
 *  param   N               The array length (number of elements)
 *
 *  \return The lowest altitude at which the wetbulb temperature is zero
*/
[[nodiscard]] float interp_wbzh(const float p_arr[], const float t_arr[], const float td_arr[],
                          const float height_arr[], const int N) noexcept;

}  // end namespace sharp

#endif // __SHARP_INTERP_H__
