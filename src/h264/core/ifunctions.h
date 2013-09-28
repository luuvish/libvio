#ifndef _IFUNCTIONS_H_
#define _IFUNCTIONS_H_


#include <algorithm>

using std::max;
using std::min;

template<typename T>
inline T clip3(T low, T high, T x)
{
    return min(max(low, x), high);
}

template<typename T>
inline T clip1(T high, T x)
{
    return clip3(0, high, x);
}

#include <cstdlib>

using std::abs;

#include <cmath>

using std::round;
using std::ceil;
using std::log2;
using std::log10;
using std::sqrt;


template<typename T>
inline T sign(T x)
{
    return ((x >= 0) - (x < 0));
}

template<typename T>
inline T median(T x, T y, T z)
{
	return x + y + z - min(x, min(y, z)) - max(x, max(y, z));
}


#endif
