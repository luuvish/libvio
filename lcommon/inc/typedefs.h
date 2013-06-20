/*!
 *************************************************************************************
 * \file typedefs.h
 *
 * \brief
 *    Common type definitions
 *    Currently only supports Windows and Linux operating systems. 
 *    Need to add support for other "older systems such as VAX, DECC, Unix Alpha etc
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>
 *************************************************************************************
 */

#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "win32.h"

typedef unsigned char  byte;     //!< byte type definition
typedef unsigned char  uint8;    //!< type definition for unsigned char (same as byte, 8 bits)
typedef unsigned short uint16;   //!< type definition for unsigned short (16 bits)
typedef unsigned int   uint32;   //!< type definition for unsigned int (32 bits)

typedef          char  int8;
typedef          short int16;
typedef          int   int32;

typedef uint16 imgpel;

enum {
	FALSE,
	TRUE
};
#define Boolean int
/*
//! Boolean Type
#ifdef FALSE
#  define Boolean int
#else
typedef enum {
  FALSE,
  TRUE
} Boolean;
#endif
*/



#ifdef __cplusplus
}
#endif

#endif

