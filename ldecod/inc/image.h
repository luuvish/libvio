
/*!
 ************************************************************************
 * \file image.h
 *
 * \brief
 *    prototypes for image.c
 *
 ************************************************************************
 */

#ifndef _IMAGE_H_
#define _IMAGE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "mbuffer.h"

struct slice_t;

extern void init_picture(VideoParameters *p_Vid, struct slice_t *currSlice, InputParameters *p_Inp);

extern void calculate_frame_no(VideoParameters *p_Vid, StorablePicture *p);
extern void find_snr          (VideoParameters *p_Vid, StorablePicture *p, int *p_ref);

extern void decode_one_slice  (struct slice_t *currSlice);
extern int  read_new_slice    (struct slice_t *currSlice);
extern void exit_picture      (VideoParameters *p_Vid, StorablePicture **dec_picture);
extern int  decode_one_frame  (DecoderParams *pDecoder);

extern void init_old_slice(OldSliceParams *p_old_slice);
// For 4:4:4 independent mode
extern void copy_dec_picture_JV (VideoParameters *p_Vid, StorablePicture *dst, StorablePicture *src );

extern void frame_postprocessing(VideoParameters *p_Vid);
extern void field_postprocessing(VideoParameters *p_Vid);

#if (MVC_EXTENSION_ENABLE)
extern int GetViewIdx(VideoParameters *p_Vid, int iVOIdx);
extern int GetVOIdx(VideoParameters *p_Vid, int iViewId);
extern int get_maxViewIdx(VideoParameters *p_Vid, int view_id, int anchor_pic_flag, int listidx);
#endif

extern void init_slice(VideoParameters *p_Vid, struct slice_t *currSlice);
extern void decode_slice(struct slice_t *currSlice, int current_header);

#ifdef __cplusplus
}
#endif

#endif

