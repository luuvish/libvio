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
 *  File      : inter_prediction.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "inter_prediction.h"


namespace vio  {
namespace h264 {


static inline int RSHIFT_RND(int x, int a)
{
    return (a > 0) ? ((x + (1 << (a-1) )) >> a) : (x << (-a));
}
static inline int RSHIFT_RND_SF(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}


InterPrediction::InterPrediction()
{
}

InterPrediction::~InterPrediction()
{
}


void InterPrediction::init(slice_t& slice)
{
}

void InterPrediction::mc_prediction(
    imgpel* mb_pred, imgpel block[16][16], int block_size_y, int block_size_x,
    mb_t* mb, ColorPlane pl, short l0_refframe, int pred_dir)
{
    slice_t& slice = *mb->p_Slice;
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;
    shr_t& shr = slice.header;

    bool weighted_pred_flag = (pps.weighted_pred_flag && (shr.slice_type == P_slice || shr.slice_type == SP_slice)) ||
                              (pps.weighted_bipred_idc == 1 && (shr.slice_type == B_slice));
    int weight, offset, denom, color_clip;
    if (weighted_pred_flag) {
        int ref_idx = shr.MbaffFrameFlag && mb->mb_field_decoding_flag ? l0_refframe >> 1 : l0_refframe;

        weight = pl == 0 ? (pred_dir == 0 ? shr.luma_weight_l0 : shr.luma_weight_l1)[ref_idx] :
                           (pred_dir == 0 ? shr.chroma_weight_l0 : shr.chroma_weight_l1)[ref_idx][pl - 1];
        offset = pl == 0 ? (pred_dir == 0 ? shr.luma_offset_l0 : shr.luma_offset_l1)[ref_idx] :
                           (pred_dir == 0 ? shr.chroma_offset_l0 : shr.chroma_offset_l1)[ref_idx][pl - 1];
        offset <<= pl == 0 ? sps.bit_depth_luma_minus8 : sps.bit_depth_chroma_minus8;

        denom  = pl > 0 ? shr.chroma_log2_weight_denom : shr.luma_log2_weight_denom;
        color_clip = (1 << (pl > 0 ? sps.BitDepthC : sps.BitDepthY)) - 1;
    }

    for (int j = 0; j < block_size_y; ++j) {
        for (int i = 0; i < block_size_x; ++i) {
            if (weighted_pred_flag) {
                int result = RSHIFT_RND((weight * block[j][i]), denom);
                mb_pred[i] = (imgpel) clip1(color_clip, result + offset);
            } else
                mb_pred[i] = block[j][i];
        }
        mb_pred += 16;
    }
}

void InterPrediction::bi_prediction(
    imgpel* mb_pred, imgpel block_l0[16][16], imgpel block_l1[16][16],
    int block_size_y, int block_size_x, mb_t* mb, ColorPlane pl, short l0_refframe, short l1_refframe)
{
    slice_t& slice = *mb->p_Slice;
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;
    shr_t& shr = slice.header;

    int weight0, weight1, offset0, offset1;
    int offset, denom, color_clip;
    if (pps.weighted_bipred_idc) {
        int ref_idx0 = shr.MbaffFrameFlag && mb->mb_field_decoding_flag ? l0_refframe >> 1 : l0_refframe;
        int ref_idx1 = shr.MbaffFrameFlag && mb->mb_field_decoding_flag ? l1_refframe >> 1 : l1_refframe;

        if (pps.weighted_bipred_idc == 1) {
            weight0 = pl == 0 ? shr.luma_weight_l0[ref_idx0] : shr.chroma_weight_l0[ref_idx0][pl - 1];
            weight1 = pl == 0 ? shr.luma_weight_l1[ref_idx1] : shr.chroma_weight_l1[ref_idx1][pl - 1];
            offset0 = pl == 0 ? shr.luma_offset_l0[ref_idx0] : shr.chroma_offset_l0[ref_idx0][pl - 1];
            offset1 = pl == 0 ? shr.luma_offset_l1[ref_idx1] : shr.chroma_offset_l1[ref_idx1][pl - 1];
            offset0 <<= pl == 0 ? sps.bit_depth_luma_minus8 : sps.bit_depth_chroma_minus8;
            offset1 <<= pl == 0 ? sps.bit_depth_luma_minus8 : sps.bit_depth_chroma_minus8;
        }

        if (pps.weighted_bipred_idc == 2) {
            storable_picture* ref_pic0 = slice.RefPicList[LIST_0][ref_idx0];
            storable_picture* ref_pic1 = slice.RefPicList[LIST_1][ref_idx1];
            if (shr.MbaffFrameFlag && mb->mb_field_decoding_flag) {
                ref_pic0 = (mb->mbAddrX % 2 == l0_refframe % 2) ? ref_pic0->top_field : ref_pic0->bottom_field;
                ref_pic1 = (mb->mbAddrX % 2 == l1_refframe % 2) ? ref_pic1->top_field : ref_pic1->bottom_field;
            }

            int td = clip3(-128, 127, ref_pic1->poc - ref_pic0->poc);
            if (td == 0 || ref_pic1->is_long_term || ref_pic0->is_long_term) {
                weight0 = 32;
                weight1 = 32;
            } else {
                int poc = !shr.MbaffFrameFlag || !mb->mb_field_decoding_flag ? shr.PicOrderCnt :
                          mb->mbAddrX % 2 == 0 ? shr.TopFieldOrderCnt : shr.BottomFieldOrderCnt;
                int tb = clip3(-128, 127, poc - ref_pic0->poc);
                int tx = (16384 + abs(td / 2)) / td;
                int DistScaleFactor = clip3(-1024, 1023, (tx * tb + 32) >> 6);

                weight1 = DistScaleFactor >> 2;
                weight0 = 64 - weight1;
                if (weight1 < -64 || weight1 > 128) {
                    weight0 = 32;
                    weight1 = 32;
                }
            }
            offset0 = offset1 = 0;
        }

        offset = (offset0 + offset1 + 1) >> 1;
        denom  = pl > 0 ? shr.chroma_log2_weight_denom + 1 : shr.luma_log2_weight_denom + 1;
        color_clip = (1 << (pl > 0 ? sps.BitDepthC : sps.BitDepthY)) - 1;
    }

    for (int j = 0; j < block_size_y; ++j) {
        for (int i = 0; i < block_size_x; ++i) {
            if (pps.weighted_bipred_idc) {
                int result = RSHIFT_RND((weight0 * block_l0[j][i] + weight1 * block_l1[j][i]), denom);
                mb_pred[i] = (imgpel) clip1(color_clip, result + offset);
            } else
                mb_pred[i] = (imgpel) (((block_l0[j][i] + block_l1[j][i]) + 1) >> 1);
        }
        mb_pred += 16;
    }
}

void InterPrediction::get_block_luma(
    storable_picture* curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y,
    imgpel block[16][16], int shift_x, int maxold_x, int maxold_y, ColorPlane pl, mb_t* currMB)
{
    slice_t& slice = *currMB->p_Slice;
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;
    int max_imgpel_value = (1 << (pl > 0 ? sps.BitDepthC : sps.BitDepthY)) - 1;

    if (curr_ref->no_ref) {
        memset(block[0], max_imgpel_value, block_size_y * block_size_x * sizeof(imgpel));
        return;
    }

    imgpel **cur_imgY = (sps.separate_colour_plane_flag && shr.colour_plane_id > PLANE_Y) ?
                        curr_ref->imgUV[shr.colour_plane_id-1] : 
                        pl ? curr_ref->imgUV[pl - 1] : curr_ref->imgY;
    int dx = (x_pos & 3);
    int dy = (y_pos & 3);
    x_pos >>= 2;
    y_pos >>= 2;
    x_pos = clip3(-18, maxold_x + 2, x_pos);
    y_pos = clip3(-10, maxold_y + 2, y_pos);

    int tmp_res[16+5][16+5];

    if (dx == 0 && dy == 0) {
        for (int j = 0; j < block_size_y; ++j)
            memcpy(&block[j][0], &cur_imgY[y_pos + j][x_pos], 16 * sizeof(imgpel));
    } else if (dx == 0 || dy == 0) {
        if (dy == 0)
            shift_x = 1;

        int img_pos_x = x_pos - (dy == 0 ? 2 : 0);
        int img_pos_y = y_pos - (dx == 0 ? 2 : 0);
        int blk_pos_x = x_pos + (dx == 3 && dy == 0 ? 1 : 0);
        int blk_pos_y = y_pos + (dx == 0 && dy == 3 ? 1 : 0);

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                imgpel p0 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 0];
                imgpel p1 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 1];
                imgpel p2 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 2];
                imgpel p3 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 3];
                imgpel p4 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 4];
                imgpel p5 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                block[j][i] = (imgpel) clip1(max_imgpel_value, (result + 16) >> 5);

                if (dx == 1 || dx == 3 || dy == 1 || dy == 3) {
                    imgpel q0 = cur_imgY[blk_pos_y + j][blk_pos_x + i];
                    block[j][i] = (imgpel) ((block[j][i] + q0 + 1) >> 1);
                }
            }
        }
    } else if (dx == 2 || dy == 2) {
        if (dx == 2)
            shift_x = 1;
        int step_x = dx != 2 && dy == 2 ? 1 : 0;
        int step_y = dx == 2 || dy != 2 ? 1 : 0;

        int img_pos_x = x_pos - 2;
        int img_pos_y = y_pos - 2;
        int blk_pos_x = dy == 2 ? (dx == 3 ? 3 : 2) : 0;
        int blk_pos_y = dx == 2 ? (dy == 3 ? 3 : 2) : 0;

        for (int j = 0; j < block_size_y + (dx == 2 ? 5 : 0); ++j) {
            for (int i = 0; i < block_size_x + (dy == 2 ? 5 : 0); ++i) {
                imgpel p0 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 0];
                imgpel p1 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 1];
                imgpel p2 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 2];
                imgpel p3 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 3];
                imgpel p4 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 4];
                imgpel p5 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);

                tmp_res[j][i] = result;
            }
        }

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                int p0 = tmp_res[j + step_y * 0][i + step_x * 0];
                int p1 = tmp_res[j + step_y * 1][i + step_x * 1];
                int p2 = tmp_res[j + step_y * 2][i + step_x * 2];
                int p3 = tmp_res[j + step_y * 3][i + step_x * 3];
                int p4 = tmp_res[j + step_y * 4][i + step_x * 4];
                int p5 = tmp_res[j + step_y * 5][i + step_x * 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                block[j][i] = (imgpel) clip1(max_imgpel_value, (result + 512) >> 10);

                if (dx == 1 || dx == 3 || dy == 1 || dy == 3) {
                    imgpel q0 = clip1(max_imgpel_value, (tmp_res[blk_pos_y + j][blk_pos_x + i] + 16) >> 5);
                    block[j][i] = (imgpel) ((block[j][i] + q0 + 1) >> 1);
                }
            }
        }
    } else {
        int img_pos_x = x_pos - 2;
        int img_pos_y = y_pos + (dy == 3 ? 1 : 0);
        int blk_pos_x = x_pos + (dx == 3 ? 1 : 0);
        int blk_pos_y = y_pos - 2;

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                imgpel p0 = cur_imgY[img_pos_y + j][img_pos_x + i + 0];
                imgpel p1 = cur_imgY[img_pos_y + j][img_pos_x + i + 1];
                imgpel p2 = cur_imgY[img_pos_y + j][img_pos_x + i + 2];
                imgpel p3 = cur_imgY[img_pos_y + j][img_pos_x + i + 3];
                imgpel p4 = cur_imgY[img_pos_y + j][img_pos_x + i + 4];
                imgpel p5 = cur_imgY[img_pos_y + j][img_pos_x + i + 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                block[j][i] = (imgpel) clip1(max_imgpel_value, (result + 16) >> 5);
            }
        }

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                imgpel p0 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 0];
                imgpel p1 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 1];
                imgpel p2 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 2];
                imgpel p3 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 3];
                imgpel p4 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 4];
                imgpel p5 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                imgpel q0 = (imgpel) clip1(max_imgpel_value, (result + 16) >> 5);

                block[j][i] = (imgpel) ((block[j][i] + q0 + 1) >> 1);
            }
        }
    }
}

void InterPrediction::get_block_chroma(
    storable_picture* curr_ref, int x_pos, int y_pos,
    int maxold_x, int maxold_y, int block_size_x, int block_size_y,
    imgpel block1[16][16], imgpel block2[16][16], mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    imgpel no_ref_value = (imgpel)(1 << (sps->BitDepthC - 1));

    int shiftpel_x = sps->chroma_format_idc == YUV400 ? 0 :
                     sps->chroma_format_idc == YUV444 ? 2 : 3;
    int shiftpel_y = sps->chroma_format_idc == YUV400 ? 0 :
                     sps->chroma_format_idc == YUV420 ? 3 : 2;
    int total_scale = shiftpel_x + shiftpel_y;

    int subpel_x = sps->chroma_format_idc == YUV400 ? 0 :
                   sps->chroma_format_idc == YUV444 ? 3 : 7;
    int subpel_y = sps->chroma_format_idc == YUV400 ? 0 :
                   sps->chroma_format_idc == YUV420 ? 7 : 3;

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (sps->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (sps->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

    if (curr_ref->no_ref) {
        memset(block1, no_ref_value, block_size_y * block_size_x * sizeof(imgpel));
        memset(block2, no_ref_value, block_size_y * block_size_x * sizeof(imgpel));
        return;
    }

    short dx = (short) (x_pos & subpel_x);
    short dy = (short) (y_pos & subpel_y);
    x_pos = x_pos >> shiftpel_x;
    y_pos = y_pos >> shiftpel_y;
    //clip MV;
    assert(block_size_y <= iChromaPadY && block_size_x <= iChromaPadX);
    x_pos = clip3(-iChromaPadX, maxold_x, x_pos); //16
    y_pos = clip3(-iChromaPadY, maxold_y, y_pos); //8

    if (dx == 0 && dy == 0) {
        for (int j = 0; j < block_size_y; ++j) {
            memcpy(&block1[j][0], &curr_ref->imgUV[0][y_pos + j][x_pos], 16 * sizeof(imgpel));
            memcpy(&block2[j][0], &curr_ref->imgUV[1][y_pos + j][x_pos], 16 * sizeof(imgpel));
        }
    } else {
        short dxcur = (short) (subpel_x + 1 - dx);
        short dycur = (short) (subpel_y + 1 - dy);
        short w00 = dxcur * dycur;
        short w01 = dxcur * dy;
        short w10 = dx * dycur;
        short w11 = dx * dy;

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                imgpel p0 = curr_ref->imgUV[0][y_pos + j + 0][x_pos + i + 0];
                imgpel p1 = curr_ref->imgUV[0][y_pos + j + 1][x_pos + i + 0];
                imgpel p2 = curr_ref->imgUV[0][y_pos + j + 0][x_pos + i + 1];
                imgpel p3 = curr_ref->imgUV[0][y_pos + j + 1][x_pos + i + 1];
                int result = (w00 * p0 + w01 * p1 + w10 * p2 + w11 * p3);
                block1[j][i] = (imgpel) RSHIFT_RND_SF(result, total_scale);
            }
        }
        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                imgpel p0 = curr_ref->imgUV[1][y_pos + j + 0][x_pos + i + 0];
                imgpel p1 = curr_ref->imgUV[1][y_pos + j + 1][x_pos + i + 0];
                imgpel p2 = curr_ref->imgUV[1][y_pos + j + 0][x_pos + i + 1];
                imgpel p3 = curr_ref->imgUV[1][y_pos + j + 1][x_pos + i + 1];
                int result = (w00 * p0 + w01 * p1 + w10 * p2 + w11 * p3);
                block2[j][i] = (imgpel) RSHIFT_RND_SF(result, total_scale);
            }
        }
    }
}


void InterPrediction::check_motion_vector_range(mb_t& mb, const mv_t *mv, slice_t *pSlice)
{
    int max_vmv_r;
    const sps_t& sps = *pSlice->active_sps;
    const shr_t& shr = pSlice->header;
    if (sps.level_idc <= 10)
        max_vmv_r = 64 * 4;
    else if (sps.level_idc <= 20)
        max_vmv_r = 128 * 4;
    else if (sps.level_idc <= 30)
        max_vmv_r = 256 * 4;
    else
        max_vmv_r = 512 * 4; // 512 pixels in quarter pixels
    max_vmv_r = shr.field_pic_flag || mb.mb_field_decoding_flag ? max_vmv_r >> 1 : max_vmv_r;

    if (mv->mv_x > 8191 || mv->mv_x < -8192)
        fprintf(stderr, "WARNING! Horizontal motion vector %d is out of allowed range {-8192, 8191} in picture %d, macroblock %d\n",
                mv->mv_x, pSlice->p_Vid->number, mb.mbAddrX);
    if ((mv->mv_y > max_vmv_r - 1) || (mv->mv_y < -max_vmv_r))
        fprintf(stderr, "WARNING! Vertical motion vector %d is out of allowed range {%d, %d} in picture %d, macroblock %d\n",
                mv->mv_y, -max_vmv_r, max_vmv_r - 1,
                pSlice->p_Vid->number, mb.mbAddrX);
}

int InterPrediction::CheckVertMV(mb_t *currMB, int vec_y, int block_size_y)
{
    storable_picture *dec_picture = currMB->p_Slice->dec_picture;
    int y_pos = vec_y >> 2;
    int maxold_y = (currMB->mb_field_decoding_flag) ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y - 1;

    if (block_size_y <= (MCBUF_LUMA_PAD_Y-4))
        return 0;

    if (y_pos < (2 - MCBUF_LUMA_PAD_Y) ||
        y_pos > (maxold_y + MCBUF_LUMA_PAD_Y - block_size_y - 2))
        return 1;
    else
        return 0;
}

static int chroma_vector_adjustment(mb_t& mb, storable_picture* RefPic)
{
    slice_t& slice = *mb.p_Slice;
    shr_t& shr = slice.header;

    if (!shr.MbaffFrameFlag && shr.field_pic_flag) {
        if (shr.structure != RefPic->slice.structure)
            return !shr.bottom_field_flag ? -2 : 2;
    }
    if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag) {
        if (mb.mbAddrX % 2 == 0 && RefPic->slice.structure == BOTTOM_FIELD)
            return -2;
        if (mb.mbAddrX % 2 == 1 && RefPic->slice.structure == TOP_FIELD)
            return 2;
    }

    return 0;
}

void InterPrediction::perform_mc(mb_t* mb, ColorPlane pl, int pred_dir, int i, int j, int block_size_x, int block_size_y)
{
    slice_t& slice = *mb->p_Slice;
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;

    storable_picture* dec_picture = slice.dec_picture;

    mv_t *l0_mv_array, *l1_mv_array;
    short l0_refframe, l1_refframe;
    storable_picture* RefPicList0, *RefPicList1;

    auto mv_info = &dec_picture->mv_info[mb->mb.y * 4 + j][mb->mb.x * 4 + i];
    if (pred_dir != 2) {
        l0_mv_array = &mv_info->mv[pred_dir];
        l0_refframe = mv_info->ref_idx[pred_dir];
        RefPicList0 = get_ref_pic(*mb, slice.RefPicList[pred_dir], l0_refframe);
    } else {
        l0_mv_array = &mv_info->mv[LIST_0];
        l1_mv_array = &mv_info->mv[LIST_1];
        l0_refframe = mv_info->ref_idx[LIST_0];
        l1_refframe = mv_info->ref_idx[LIST_1];
        RefPicList0 = get_ref_pic(*mb, slice.RefPicList[LIST_0], l0_refframe);
        RefPicList1 = get_ref_pic(*mb, slice.RefPicList[LIST_1], l1_refframe);
    }

    int block_y_aff;
    if (shr.MbaffFrameFlag && mb->mb_field_decoding_flag)
        block_y_aff = (mb->mb.y * 4 - 4 * (mb->mbAddrX % 2)) / 2;
    else
        block_y_aff = mb->mb.y * 4;

    check_motion_vector_range(*mb, l0_mv_array, &slice);
    int vec1_x = (mb->mb.x * 4 + i) * 16 + l0_mv_array->mv_x;
    int vec1_y = (block_y_aff + j) * 16 + l0_mv_array->mv_y;
    int vec2_x, vec2_y;
    if (pred_dir == 2) {
        check_motion_vector_range(*mb, l1_mv_array, &slice);
        vec2_x = (mb->mb.x * 4 + i) * 16 + l1_mv_array->mv_x;
        vec2_y = (block_y_aff + j) * 16 + l1_mv_array->mv_y;
    }

    // vars for get_block_luma
    int shift_x  = dec_picture->iLumaStride;
    int maxold_x = dec_picture->size_x - 1;
    int maxold_y = (mb->mb_field_decoding_flag) ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y - 1;
    imgpel tmp_block_l0[16][16];
    imgpel tmp_block_l1[16][16];
    imgpel tmp_block_l2[16][16];
    imgpel tmp_block_l3[16][16];

    if (CheckVertMV(mb, vec1_y, block_size_y)) {
        get_block_luma(RefPicList0, vec1_x, vec1_y,
                       block_size_x, 8, tmp_block_l0, shift_x,
                       maxold_x, maxold_y, pl, mb);
        get_block_luma(RefPicList0, vec1_x, vec1_y+32,
                       block_size_x, block_size_y-8, tmp_block_l0 + 8, shift_x,
                       maxold_x, maxold_y, pl, mb);
    } else
        get_block_luma(RefPicList0, vec1_x, vec1_y,
                       block_size_x, block_size_y, tmp_block_l0, shift_x,
                       maxold_x, maxold_y, pl, mb);
    if (pred_dir == 2) {
        if (CheckVertMV(mb, vec2_y, block_size_y)) {
            get_block_luma(RefPicList1, vec2_x, vec2_y,
                           block_size_x, 8, tmp_block_l1, shift_x,
                           maxold_x, maxold_y, pl, mb);
            get_block_luma(RefPicList1, vec2_x, vec2_y+32,
                           block_size_x, block_size_y-8, tmp_block_l1 + 8, shift_x,
                           maxold_x, maxold_y, pl, mb);
        } else
            get_block_luma(RefPicList1, vec2_x, vec2_y,
                           block_size_x, block_size_y, tmp_block_l1, shift_x,
                           maxold_x, maxold_y, pl, mb);
    }

    if (pred_dir != 2)
        mc_prediction(&slice.mb_pred[pl][j * 4][i * 4],
                      tmp_block_l0, block_size_y, block_size_x,
                      mb, pl, l0_refframe, pred_dir); 
    else
        bi_prediction(&slice.mb_pred[pl][j * 4][i * 4],
                      tmp_block_l0, tmp_block_l1, block_size_y, block_size_x,
                      mb, pl, l0_refframe, l1_refframe);

    if (sps.chroma_format_idc != YUV400 && sps.chroma_format_idc != YUV444) {
        int maxold_x = dec_picture->size_x_cr - 1;
        int maxold_y = mb->mb_field_decoding_flag ? (dec_picture->size_y_cr >> 1) - 1 : dec_picture->size_y_cr - 1;

        int ioff_cr = sps.MbWidthC  == 16 ? i * 4 : i * 2;
        int joff_cr = sps.MbHeightC == 16 ? j * 4 : j * 2;
        int block_size_x_cr = sps.MbWidthC  == 16 ? block_size_x : block_size_x >> 1;
        int block_size_y_cr = sps.MbHeightC == 16 ? block_size_y : block_size_y >> 1;

        int vec2_y_cr, vec1_y_cr;
        if (pred_dir != 2) {
            if (sps.chroma_format_idc == 1)
                vec1_y_cr = vec1_y + vio::h264::chroma_vector_adjustment(*mb, RefPicList0);
            else
                vec1_y_cr = vec1_y;
        } else {
            if (sps.chroma_format_idc == 1) {
                vec1_y_cr = vec1_y + vio::h264::chroma_vector_adjustment(*mb, RefPicList0);
                vec2_y_cr = vec2_y + vio::h264::chroma_vector_adjustment(*mb, RefPicList1);
            } else {
                vec1_y_cr = vec1_y;
                vec2_y_cr = vec2_y;
            }
        }

        get_block_chroma(RefPicList0, vec1_x, vec1_y_cr,
                         maxold_x, maxold_y, block_size_x_cr, block_size_y_cr,
                         tmp_block_l0, tmp_block_l2, mb);
        if (pred_dir == 2)
            get_block_chroma(RefPicList1, vec2_x, vec2_y_cr,
                             maxold_x, maxold_y, block_size_x_cr, block_size_y_cr,
                             tmp_block_l1, tmp_block_l3, mb);

        if (pred_dir != 2) {
            mc_prediction(&slice.mb_pred[1][joff_cr][ioff_cr],
                          tmp_block_l0, block_size_y_cr, block_size_x_cr,
                          mb, PLANE_U, l0_refframe, pred_dir);
            mc_prediction(&slice.mb_pred[2][joff_cr][ioff_cr],
                          tmp_block_l2, block_size_y_cr, block_size_x_cr,
                          mb, PLANE_V, l0_refframe, pred_dir);
        } else {
            bi_prediction(&slice.mb_pred[1][joff_cr][ioff_cr],
                          tmp_block_l0, tmp_block_l1, block_size_y_cr, block_size_x_cr,
                          mb, PLANE_U, l0_refframe, l1_refframe);
            bi_prediction(&slice.mb_pred[2][joff_cr][ioff_cr],
                          tmp_block_l2, tmp_block_l3, block_size_y_cr, block_size_x_cr,
                          mb, PLANE_V, l0_refframe, l1_refframe);
        }
    }
}


void InterPrediction::motion_compensation(mb_t *mb)
{
}


}
}
