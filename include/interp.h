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
#pragma once 
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
 * in the paper located here: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html 

 *
 * This linear interpolation guarantees the following:
 *
 * 1. Exactness: lerp(a,b,0)==a & lerp(a,b,1)==b
 * 2. Monotonicity: cmp(lerp(a,b,t2), lerp(a,b,t1)) * cmp(t2, t1) * cmp(b,a) >= 0, where cmp is an arithmetic three-way comparison function
 * 3. Determinacy: result of NaN only for lerp(a,a,INFINITY)
 * 4. Boundedness: t<0 || t>1 || isfinite(lerp(a,b,t))
 * 5. Consistency: lerp(a,a,t) == a
 *
 * \param __a lef value
 * \param __b right value
 * \param __t where 0 <= __t <= 1.0
 *
 * \return The value between __a and __b at distance __t between them. 
 */
template<typename _Fp>
 constexpr _Fp
 __lerp(_Fp __a, _Fp __b, _Fp __t) noexcept
 {
    if ((__a <= 0 && __b >= 0) || (__a >= 0 && __b <= 0))
        return __t * __b + (1 - __t) * __a;

    if (__t == 1)
        return __b;                        // exact

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
 * in the paper located here: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html 

 *
 * This linear interpolation guarantees the following:
 *
 * 1. Exactness: lerp(a,b,0)==a & lerp(a,b,1)==b
 * 2. Monotonicity: cmp(lerp(a,b,t2), lerp(a,b,t1)) * cmp(t2, t1) * cmp(b,a) >= 0, where cmp is an arithmetic three-way comparison function
 * 3. Determinacy: result of NaN only for lerp(a,a,INFINITY)
 * 4. Boundedness: t<0 || t>1 || isfinite(lerp(a,b,t))
 * 5. Consistency: lerp(a,a,t) == a
 *
 * \param __a lef value
 * \param __b right value
 * \param __t where 0 <= __t <= 1.0
 *
 * \return The value between __a and __b at distance __t between them. 
 */
constexpr float
lerp(float __a, float __b, float __t) noexcept
{ return __lerp(__a, __b, __t); }

/**
 * \author C++20 Standard Template Library <cmath>
 *
 * \brief Linearly interpolates between a and b with distance t.
 *
 * This routine was copied from the C++20 <cmath> standard template library.
 * It is used to linearly interpolate between A and B over a normalized 
 * distance T, where 0 <= T <= 1. Full documentation and details can be found
 * in the paper located here: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html 

 *
 * This linear interpolation guarantees the following:
 *
 * 1. Exactness: lerp(a,b,0)==a & lerp(a,b,1)==b
 * 2. Monotonicity: cmp(lerp(a,b,t2), lerp(a,b,t1)) * cmp(t2, t1) * cmp(b,a) >= 0, where cmp is an arithmetic three-way comparison function
 * 3. Determinacy: result of NaN only for lerp(a,a,INFINITY)
 * 4. Boundedness: t<0 || t>1 || isfinite(lerp(a,b,t))
 * 5. Consistency: lerp(a,a,t) == a
 *
 * \param __a lef value
 * \param __b right value
 * \param __t where 0 <= __t <= 1.0
 *
 * \return The value between __a and __b at distance __t between them. 
 */
constexpr double
lerp(double __a, double __b, double __t) noexcept
{ return __lerp(__a, __b, __t); }

/**
 * \author C++20 Standard Template Library <cmath>
 *
 * \brief Linearly interpolates between a and b with distance t.
 *
 * This routine was copied from the C++20 <cmath> standard template library.
 * It is used to linearly interpolate between A and B over a normalized 
 * distance T, where 0 <= T <= 1. Full documentation and details can be found
 * in the paper located here: https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0811r3.html 

 *
 * This linear interpolation guarantees the following:
 *
 * 1. Exactness: lerp(a,b,0)==a & lerp(a,b,1)==b
 * 2. Monotonicity: cmp(lerp(a,b,t2), lerp(a,b,t1)) * cmp(t2, t1) * cmp(b,a) >= 0, where cmp is an arithmetic three-way comparison function
 * 3. Determinacy: result of NaN only for lerp(a,a,INFINITY)
 * 4. Boundedness: t<0 || t>1 || isfinite(lerp(a,b,t))
 * 5. Consistency: lerp(a,a,t) == a
 *
 * \param __a lef value
 * \param __b right value
 * \param __t where 0 <= __t <= 1.0
 *
 * \return The value between __a and __b at distance __t between them. 
 */
constexpr long double
lerp(long double __a, long double __b, long double __t) noexcept
{ return __lerp(__a, __b, __t); }

}
