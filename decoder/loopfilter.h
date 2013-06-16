/*!
 ************************************************************************
 *  \file
 *     loopfilter.h
 *  \brief
 *     external deblocking filter interface
 ************************************************************************
 */

#ifndef _LOOPFILTER_H_
#define _LOOPFILTER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "mbuffer.h"

// For 4:4:4 independent mode
extern void change_plane_JV      ( VideoParameters *p_Vid, int nplane, Slice *pSlice);
extern void make_frame_picture_JV( VideoParameters *p_Vid );

extern void DeblockPicture(VideoParameters *p_Vid, StorablePicture *p) ;

void  init_Deblock(VideoParameters *p_Vid, int mb_aff_frame_flag);

#ifdef __cplusplus
}
#endif

#endif //_LOOPFILTER_H_
