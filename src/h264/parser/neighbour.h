/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : neighbour.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _VIO_H264_NEIGHBOUR_H_
#define _VIO_H264_NEIGHBOUR_H_


namespace vio  {
namespace h264 {


struct position_t {
    int x;
    int y;
};

using pos_t = position_t;

struct location_t {
    int x;
    int y;
};

using loc_t = location_t;

loc_t operator + (const loc_t& l, const loc_t& r);
loc_t operator - (const loc_t& l, const loc_t& r);

struct neighbour_t {
    mb_t* mb;
    int   x;
    int   y;
};

using nb_t = neighbour_t;


class Neighbour {
public:
    mb_t*       mb_data;

    loc_t get_location (slice_t* slice, bool chroma, int mbAddr, const pos_t& offset={0,0});
    mb_t* get_mb       (slice_t* slice, bool chroma, int mbAddr, const pos_t& offset={0,0});
    nb_t  get_neighbour(slice_t* slice, bool chroma, int mbAddr, const pos_t& offset={0,0});
    mb_t* get_mb       (slice_t* slice, bool chroma, const loc_t& loc);
    nb_t  get_neighbour(slice_t* slice, bool chroma, const loc_t& loc);

    pos_t get_position(slice_t* slice, int mbAddr, int blkIdx);

    int   predict_nnz(mb_t* mb, int pl, int i, int j);
};


int mb_skip_flag_ctxIdxInc(mb_t& mb);
int mb_field_decoding_flag_ctxIdxInc(mb_t& mb);
int mb_type_si_slice_ctxIdxInc(mb_t& mb);
int mb_type_i_slice_ctxIdxInc(mb_t& mb);
int mb_type_b_slice_ctxIdxInc(mb_t& mb);
int transform_size_8x8_flag_ctxIdxInc(mb_t& mb);
int intra_chroma_pred_mode_ctxIdxInc(mb_t& mb);
int ref_idx_ctxIdxInc(mb_t& mb, uint8_t list, uint8_t x0, uint8_t y0);
int mvd_ctxIdxInc(mb_t& mb, uint8_t list, uint8_t x0, uint8_t y0, bool comp);
int cbp_luma_ctxIdxInc(mb_t& mb, uint8_t x0, uint8_t y0, uint8_t coded_block_pattern);
int cbp_chroma_ctxIdxInc(mb_t& mb);

int coded_block_flag_ctxIdxInc(mb_t& mb, int pl, bool chroma, bool ac, int blkIdx);
void update_coded_block_flag(mb_t* mb, int pl, bool chroma, bool ac, int blkIdx);


}
}


#endif // _VIO_H264_NEIGHBOUR_H_
