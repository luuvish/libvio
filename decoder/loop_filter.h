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
#include "dpb.h"

struct slice_t;

// For 4:4:4 independent mode
void change_plane_JV      ( VideoParameters *p_Vid, int nplane, struct slice_t *pSlice);
void make_frame_picture_JV( VideoParameters *p_Vid );

void init_Deblock(VideoParameters *p_Vid, int mb_aff_frame_flag);
void DeblockPicture(VideoParameters *p_Vid, StorablePicture *p) ;

#ifdef __cplusplus
}
#endif

#endif //_LOOPFILTER_H_
