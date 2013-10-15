#ifndef _DEFINES_H_
#define _DEFINES_H_


#define DISABLE_ERC               0    //!< Disable any error concealment processes

#define MVC_EXTENSION_ENABLE      1    //!< enable support for the Multiview High Profile
#define MAX_VIEW_NUM              1024   

typedef uint8_t  byte;
typedef uint16_t px_t;

#define MAX_CODED_FRAME_SIZE   8000000         //!< bytes for one frame
#define MCBUF_LUMA_PAD_X        32
#define MCBUF_LUMA_PAD_Y        12
#define MCBUF_CHROMA_PAD_X      16
#define MCBUF_CHROMA_PAD_Y      8


#include <algorithm>

using std::max;
using std::min;

template<typename T>
inline T max(T x, T y, T z)
{
	return max(max(x, y), z);
}

template<typename T>
inline T min(T x, T y, T z)
{
	return min(min(x, y), z);
}

template<typename T>
inline T min_positive(T x, T y)
{
	return (x >=0 && y >= 0) ? min(x, y) : max(x, y);
}

template<typename T>
inline T min_positive(T x, T y, T z)
{
	return min_positive(min_positive(x, y), z);
}

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
