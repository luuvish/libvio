#ifndef _IFUNCTIONS_H_
#define _IFUNCTIONS_H_


#include <math.h>
#include <limits.h>


static inline int imin(int a, int b)
{
  return ((a) < (b)) ? (a) : (b);
}

static inline int imax(int a, int b)
{
  return ((a) > (b)) ? (a) : (b);
}

static inline int iabs(int x)
{
  static const int INT_BITS = (sizeof(int) * CHAR_BIT) - 1;
  int y = x >> INT_BITS;
  return (x ^ y) - y;
}



static inline int iClip3(int low, int high, int x)
{
  x = imax(x, low);
  x = imin(x, high);

  return x;
}


/*!
 ************************************************************************
 * \brief
 *    calculate Ceil(Log2(uiVal))
 ************************************************************************
 */
static inline unsigned CeilLog2( unsigned uiVal)
{
  unsigned uiTmp = uiVal-1;
  unsigned uiRet = 0;

  while( uiTmp != 0 )
  {
    uiTmp >>= 1;
    uiRet++;
  }
  return uiRet;
}


#endif
