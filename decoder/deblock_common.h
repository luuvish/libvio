/*!
 ***************************************************************************
 * \file
 *    loop_filter.h
 *
 * \date
 *    25 April 2010
 *
 * \brief
 *    Headerfile for loopfilter processing
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis         <alexismt@ieee.org> 

 **************************************************************************
 */

#ifndef _LOOP_FILTER_H_
#define _LOOP_FILTER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"

void get_strength_ver    (Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p);
void get_strength_hor    (Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p);
void get_strength_ver_MBAff(byte *Strength, Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p);
void get_strength_hor_MBAff(byte *Strength, Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p);

void set_loop_filter_functions_mbaff(VideoParameters *p_Vid);
void set_loop_filter_functions_normal(VideoParameters *p_Vid);


#ifdef __cplusplus
}
#endif

#endif
