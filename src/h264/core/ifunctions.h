#ifndef _IFUNCTIONS_H_
#define _IFUNCTIONS_H_


#include <algorithm>

using std::max;
using std::min;

template<typename T>
T clip1(T high, T x)
{
    return min(max(0, x), high);
}

template<typename T>
T clip3(T low, T high, T x)
{
    return min(max(low, x), high);
}

#include <cstdlib>

using std::abs;

#include <cmath>

using std::round;
using std::ceil;
using std::log2;


#endif
