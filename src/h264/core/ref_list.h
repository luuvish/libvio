#ifndef _REF_LIST_H_
#define _REF_LIST_H_

#include "global.h"
#include "bitstream_cabac.h"


struct slice_t;

extern void update_pic_num   (slice_t *currSlice);

extern void gen_pic_list_from_frame_list(bool bottom_field_flag, frame_store **fs_list, int list_idx, storable_picture **list, char *list_size, int long_term);
extern void init_lists         (slice_t *currSlice);

extern void reorder_lists   (slice_t *currSlice);
extern void init_mbaff_lists(VideoParameters *p_Vid, slice_t *currSlice);


extern void insert_picture_in_dpb(VideoParameters *p_Vid, frame_store* fs, storable_picture* p);

extern void fill_frame_num_gap(VideoParameters *p_Vid, slice_t *pSlice);

extern void pad_dec_picture(VideoParameters *p_Vid, storable_picture *dec_picture);
extern void pad_buf(imgpel *pImgBuf, int iWidth, int iHeight, int iStride, int iPadX, int iPadY);


#if (MVC_EXTENSION_ENABLE)
extern void reorder_lists_mvc     (slice_t *currSlice, int currPOC);
extern void init_lists_p_slice_mvc(slice_t *currSlice);
extern void init_lists_b_slice_mvc(slice_t *currSlice);
extern void init_lists_i_slice_mvc(slice_t *currSlice);

extern void reorder_ref_pic_list_mvc(slice_t *currSlice, int cur_list, int **anchor_ref, int **non_anchor_ref,
                                     int view_id, int anchor_pic_flag, int currPOC, int listidx);

extern void reorder_short_term(slice_t *currSlice, int cur_list, int num_ref_idx_lX_active_minus1, int picNumLX, int *refIdxLX, int currViewID);
extern void reorder_long_term(slice_t *currSlice, storable_picture **RefPicListX, int num_ref_idx_lX_active_minus1, int LongTermPicNum, int *refIdxLX, int currViewID);
#endif


#endif // _REF_LIST_H_
