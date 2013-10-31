/*
 * =============================================================================
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
 * =============================================================================
 *
 *  File      : inter_prediction.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * =============================================================================
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "decoder.h"


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


void InterPrediction::init(slice_t& slice)
{
    this->sets.pic   = slice.dec_picture;
    this->sets.sps   = slice.active_sps;
    this->sets.pps   = slice.active_pps;
    this->sets.slice = &slice;
}

void InterPrediction::mc_prediction(
    px_t* mb_pred, px_t block[16][16], int block_size_y, int block_size_x,
    mb_t& mb, int pl, short l0_refframe, int pred_dir)
{
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;
    shr_t& shr = slice.header;

    bool weighted_pred_flag = (pps.weighted_pred_flag && (shr.slice_type == P_slice || shr.slice_type == SP_slice)) ||
                              (pps.weighted_bipred_idc == 1 && (shr.slice_type == B_slice));
    int weight, offset, denom, color_clip;
    if (weighted_pred_flag) {
        int ref_idx = shr.MbaffFrameFlag && mb.mb_field_decoding_flag ? l0_refframe >> 1 : l0_refframe;

        weight = shr.pred_weight_l[pred_dir][pl][ref_idx].weight;
        offset = shr.pred_weight_l[pred_dir][pl][ref_idx].offset;
        offset <<= pl == 0 ? sps.bit_depth_luma_minus8 : sps.bit_depth_chroma_minus8;

        denom  = pl > 0 ? shr.chroma_log2_weight_denom : shr.luma_log2_weight_denom;
        color_clip = (1 << (pl > 0 ? sps.BitDepthC : sps.BitDepthY)) - 1;
    }

    for (int j = 0; j < block_size_y; ++j) {
        for (int i = 0; i < block_size_x; ++i) {
            if (weighted_pred_flag) {
                int result = RSHIFT_RND((weight * block[j][i]), denom);
                mb_pred[i] = (px_t) clip1(color_clip, result + offset);
            } else
                mb_pred[i] = block[j][i];
        }
        mb_pred += 16;
    }
}

void InterPrediction::bi_prediction(
    px_t* mb_pred, px_t block_l0[16][16], px_t block_l1[16][16],
    int block_size_y, int block_size_x, mb_t& mb, int pl, short l0_refframe, short l1_refframe)
{
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;
    shr_t& shr = slice.header;

    int weight0, weight1, offset0, offset1;
    int offset, denom, color_clip;
    if (pps.weighted_bipred_idc) {
        int ref_idx0 = shr.MbaffFrameFlag && mb.mb_field_decoding_flag ? l0_refframe >> 1 : l0_refframe;
        int ref_idx1 = shr.MbaffFrameFlag && mb.mb_field_decoding_flag ? l1_refframe >> 1 : l1_refframe;

        if (pps.weighted_bipred_idc == 1) {
            weight0 = shr.pred_weight_l[0][pl][ref_idx0].weight;
            weight1 = shr.pred_weight_l[1][pl][ref_idx1].weight;
            offset0 = shr.pred_weight_l[0][pl][ref_idx0].offset;
            offset1 = shr.pred_weight_l[1][pl][ref_idx1].offset;
            offset0 <<= pl == 0 ? sps.bit_depth_luma_minus8 : sps.bit_depth_chroma_minus8;
            offset1 <<= pl == 0 ? sps.bit_depth_luma_minus8 : sps.bit_depth_chroma_minus8;
        }

        if (pps.weighted_bipred_idc == 2) {
            storable_picture* ref_pic0 = slice.RefPicList[LIST_0][ref_idx0];
            storable_picture* ref_pic1 = slice.RefPicList[LIST_1][ref_idx1];
            if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag) {
                ref_pic0 = (mb.mbAddrX % 2 == l0_refframe % 2) ? ref_pic0->top_field : ref_pic0->bottom_field;
                ref_pic1 = (mb.mbAddrX % 2 == l1_refframe % 2) ? ref_pic1->top_field : ref_pic1->bottom_field;
            }

            int td = clip3(-128, 127, ref_pic1->poc - ref_pic0->poc);
            if (td == 0 || ref_pic1->is_long_term || ref_pic0->is_long_term) {
                weight0 = 32;
                weight1 = 32;
            } else {
                int poc = !shr.MbaffFrameFlag || !mb.mb_field_decoding_flag ? shr.PicOrderCnt :
                          mb.mbAddrX % 2 == 0 ? shr.TopFieldOrderCnt : shr.BottomFieldOrderCnt;
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
                mb_pred[i] = (px_t) clip1(color_clip, result + offset);
            } else
                mb_pred[i] = (px_t) (((block_l0[j][i] + block_l1[j][i]) + 1) >> 1);
        }
        mb_pred += 16;
    }
}

void InterPrediction::get_block_luma(
    storable_picture* curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y,
    px_t block[16][16], int comp, mb_t& mb)
{
    shr_t& shr = this->sets.slice->header;
    int max_imgpel_value = (1 << (comp > 0 ? this->sets.sps->BitDepthC : this->sets.sps->BitDepthY)) - 1;

    if (curr_ref->no_ref) {
        memset(block[0], max_imgpel_value, block_size_y * block_size_x * sizeof(px_t));
        return;
    }

    int shift_x  = this->sets.pic->iLumaStride;
    int maxold_x = this->sets.pic->size_x - 1;
    int maxold_y = (mb.mb_field_decoding_flag) ? (this->sets.pic->size_y >> 1) - 1 : this->sets.pic->size_y - 1;

    px_t **cur_imgY = (this->sets.sps->separate_colour_plane_flag && shr.colour_plane_id > PLANE_Y) ?
                        curr_ref->imgUV[shr.colour_plane_id-1] : 
                        comp ? curr_ref->imgUV[comp - 1] : curr_ref->imgY;
    int xFracL = (x_pos & 3);
    int yFracL = (y_pos & 3);
    int xIntL  = clip3(-18, maxold_x + 2, x_pos >> 2);
    int yIntL  = clip3(-10, maxold_y + 2, y_pos >> 2);

    int tmp_res[16+5][16+5];

    if (xFracL == 0 && yFracL == 0) {
        for (int y = 0; y < block_size_y; y++)
            memcpy(&block[y][0], &cur_imgY[yIntL + y][xIntL], 16 * sizeof(px_t));
    } else if (xFracL == 0 || yFracL == 0) {
        if (yFracL == 0)
            shift_x = 1;

        int img_pos_x = xIntL - (yFracL == 0 ? 2 : 0);
        int img_pos_y = yIntL - (xFracL == 0 ? 2 : 0);
        int blk_pos_x = xIntL + (xFracL == 3 && yFracL == 0 ? 1 : 0);
        int blk_pos_y = yIntL + (xFracL == 0 && yFracL == 3 ? 1 : 0);

        for (int y = 0; y < block_size_y; y++) {
            for (int x = 0; x < block_size_x; x++) {
                px_t p0 = cur_imgY[img_pos_y + y][img_pos_x + x + shift_x * 0];
                px_t p1 = cur_imgY[img_pos_y + y][img_pos_x + x + shift_x * 1];
                px_t p2 = cur_imgY[img_pos_y + y][img_pos_x + x + shift_x * 2];
                px_t p3 = cur_imgY[img_pos_y + y][img_pos_x + x + shift_x * 3];
                px_t p4 = cur_imgY[img_pos_y + y][img_pos_x + x + shift_x * 4];
                px_t p5 = cur_imgY[img_pos_y + y][img_pos_x + x + shift_x * 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                block[y][x] = (px_t) clip1(max_imgpel_value, (result + 16) >> 5);

                if (xFracL == 1 || xFracL == 3 || yFracL == 1 || yFracL == 3) {
                    px_t q0 = cur_imgY[blk_pos_y + y][blk_pos_x + x];
                    block[y][x] = (px_t) ((block[y][x] + q0 + 1) >> 1);
                }
            }
        }
    } else if (xFracL == 2 || yFracL == 2) {
        if (xFracL == 2)
            shift_x = 1;
        int step_x = xFracL != 2 && yFracL == 2 ? 1 : 0;
        int step_y = xFracL == 2 || yFracL != 2 ? 1 : 0;

        int img_pos_x = xIntL - 2;
        int img_pos_y = yIntL - 2;
        int blk_pos_x = yFracL == 2 ? (xFracL == 3 ? 3 : 2) : 0;
        int blk_pos_y = xFracL == 2 ? (yFracL == 3 ? 3 : 2) : 0;

        for (int j = 0; j < block_size_y + (xFracL == 2 ? 5 : 0); ++j) {
            for (int i = 0; i < block_size_x + (yFracL == 2 ? 5 : 0); ++i) {
                px_t p0 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 0];
                px_t p1 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 1];
                px_t p2 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 2];
                px_t p3 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 3];
                px_t p4 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 4];
                px_t p5 = cur_imgY[img_pos_y + j][img_pos_x + i + shift_x * 5];
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
                block[j][i] = (px_t) clip1(max_imgpel_value, (result + 512) >> 10);

                if (xFracL == 1 || xFracL == 3 || yFracL == 1 || yFracL == 3) {
                    px_t q0 = clip1(max_imgpel_value, (tmp_res[blk_pos_y + j][blk_pos_x + i] + 16) >> 5);
                    block[j][i] = (px_t) ((block[j][i] + q0 + 1) >> 1);
                }
            }
        }
    } else {
        int img_pos_x = xIntL - 2;
        int img_pos_y = yIntL + (yFracL == 3 ? 1 : 0);
        int blk_pos_x = xIntL + (xFracL == 3 ? 1 : 0);
        int blk_pos_y = yIntL - 2;

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                px_t p0 = cur_imgY[img_pos_y + j][img_pos_x + i + 0];
                px_t p1 = cur_imgY[img_pos_y + j][img_pos_x + i + 1];
                px_t p2 = cur_imgY[img_pos_y + j][img_pos_x + i + 2];
                px_t p3 = cur_imgY[img_pos_y + j][img_pos_x + i + 3];
                px_t p4 = cur_imgY[img_pos_y + j][img_pos_x + i + 4];
                px_t p5 = cur_imgY[img_pos_y + j][img_pos_x + i + 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                block[j][i] = (px_t) clip1(max_imgpel_value, (result + 16) >> 5);
            }
        }

        for (int j = 0; j < block_size_y; ++j) {
            for (int i = 0; i < block_size_x; ++i) {
                px_t p0 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 0];
                px_t p1 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 1];
                px_t p2 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 2];
                px_t p3 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 3];
                px_t p4 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 4];
                px_t p5 = cur_imgY[blk_pos_y + j][blk_pos_x + i + shift_x * 5];
                int result = (p0 - 5 * p1 + 20 * p2 + 20 * p3 - 5 * p4 + p5);
                px_t q0 = (px_t) clip1(max_imgpel_value, (result + 16) >> 5);

                block[j][i] = (px_t) ((block[j][i] + q0 + 1) >> 1);
            }
        }
    }
}

void InterPrediction::get_block_chroma(
    storable_picture* curr_ref, int x_pos, int y_pos,
    int maxold_x, int maxold_y, int partWidthC, int partHeightC,
    px_t predPartLXC[16][16], int comp, mb_t& mb)
{
    px_t** imgUV = curr_ref->imgUV[comp - 1];

    px_t no_ref_value = (px_t)(1 << (this->sets.sps->BitDepthC - 1));

    int mvLX[2] = {x_pos, y_pos};

    if (this->sets.sps->ChromaArrayType == 1) {
        slice_t& slice = *this->sets.slice;
        shr_t& shr = slice.header;

        if (!shr.MbaffFrameFlag && shr.field_pic_flag) {
            if (shr.structure != curr_ref->slice.structure)
                mvLX[1] += (!shr.bottom_field_flag ? -2 : 2);
        }
        if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag) {
            if (mb.mbAddrX % 2 == 0 && curr_ref->slice.structure == BOTTOM_FIELD)
                mvLX[1] -= 2;
            if (mb.mbAddrX % 2 == 1 && curr_ref->slice.structure == TOP_FIELD)
                mvLX[1] += 2;
        }
    }

    int mvCX[2] = {mvLX[0], mvLX[1]};

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y * (this->sets.sps->ChromaArrayType == 2 ? 2 : 1);

    assert(partWidthC <= iChromaPadX && partHeightC <= iChromaPadY);

    if (curr_ref->no_ref) {
        memset(predPartLXC, no_ref_value, partHeightC * partWidthC * sizeof(px_t));
        return;
    }

    int shiftpel_y = this->sets.sps->ChromaArrayType == 1 ? 3 : 2;
    int xAL    = clip3(-iChromaPadX, maxold_x, mvCX[0] >> 3); //16
    int yAL    = clip3(-iChromaPadY, maxold_y, mvCX[1] >> shiftpel_y); //8
    //int xAL = mvCX[0] >> 3;
    //int yAL = mvCX[1] >> shiftpel_y;
    int xFracC = (mvCX[0] & 7);
    int yFracC = (this->sets.sps->ChromaArrayType == 1) ? (mvCX[1] & 7) : (mvCX[1] & 3) << 1;

    if (xFracC == 0 && yFracC == 0) {
        for (int yC = 0; yC < partHeightC; yC++)
            memcpy(&predPartLXC[yC][0], &imgUV[yAL + yC][xAL], 16 * sizeof(px_t));
        return;
    }

    uint32_t PicWidthInSamplesC  = this->sets.sps->PicWidthInSamplesC;
    uint32_t PicHeightInSamplesC = this->sets.slice->header.PicHeightInSamplesC;
    uint32_t refPicHeightEffectiveC =
        !this->sets.slice->header.MbaffFrameFlag || !mb.mb_field_decoding_flag ?
        PicHeightInSamplesC : PicHeightInSamplesC / 2;

    for (int yC = 0; yC < partHeightC; yC++) {
        for (int xC = 0; xC < partWidthC; xC++) {
            int xIntC = xAL + xC;
            int yIntC = yAL + yC;
            int xAC = clip3<int>(0, PicWidthInSamplesC - 1, xIntC);
            int xBC = clip3<int>(0, PicWidthInSamplesC - 1, xIntC + 1);
            int xCC = clip3<int>(0, PicWidthInSamplesC - 1, xIntC);
            int xDC = clip3<int>(0, PicWidthInSamplesC - 1, xIntC + 1);
            int yAC = clip3<int>(0, refPicHeightEffectiveC - 1, yIntC);
            int yBC = clip3<int>(0, refPicHeightEffectiveC - 1, yIntC);
            int yCC = clip3<int>(0, refPicHeightEffectiveC - 1, yIntC + 1);
            int yDC = clip3<int>(0, refPicHeightEffectiveC - 1, yIntC + 1);
            px_t A = imgUV[yAC][xAC];
            px_t B = imgUV[yBC][xBC];
            px_t C = imgUV[yCC][xCC];
            px_t D = imgUV[yDC][xDC];
            predPartLXC[yC][xC] = ((8 - xFracC) * (8 - yFracC) * A +
                                   xFracC * (8 - yFracC) * B +
                                   (8 - xFracC) * yFracC * C +
                                   xFracC * yFracC * D + 32) >> 6;
        }
    }
}


void InterPrediction::check_motion_vector_range(mb_t& mb, const mv_t *mv)
{
    int max_vmv_r;
    const shr_t& shr = this->sets.slice->header;
    if (this->sets.sps->level_idc <= 10)
        max_vmv_r = 64 * 4;
    else if (this->sets.sps->level_idc <= 20)
        max_vmv_r = 128 * 4;
    else if (this->sets.sps->level_idc <= 30)
        max_vmv_r = 256 * 4;
    else
        max_vmv_r = 512 * 4; // 512 pixels in quarter pixels
    max_vmv_r = shr.field_pic_flag || mb.mb_field_decoding_flag ? max_vmv_r >> 1 : max_vmv_r;

    if (mv->mv_x > 8191 || mv->mv_x < -8192)
        fprintf(stderr, "WARNING! Horizontal motion vector %d is out of allowed range {-8192, 8191} in picture %d, macroblock %d\n",
                mv->mv_x, this->sets.slice->p_Vid->number, mb.mbAddrX);
    if ((mv->mv_y > max_vmv_r - 1) || (mv->mv_y < -max_vmv_r))
        fprintf(stderr, "WARNING! Vertical motion vector %d is out of allowed range {%d, %d} in picture %d, macroblock %d\n",
                mv->mv_y, -max_vmv_r, max_vmv_r - 1,
                this->sets.slice->p_Vid->number, mb.mbAddrX);
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

void InterPrediction::inter_pred(mb_t& mb, int comp, int pred_dir, int i, int j, int block_size_x, int block_size_y)
{
    shr_t& shr = this->sets.slice->header;

    mv_t* mv_l0, *mv_l1;
    short ref_idx_l0, ref_idx_l1;
    storable_picture* RefPicList0, *RefPicList1;

    auto mv_info = &this->sets.pic->mv_info[mb.mb.y * 4 + j][mb.mb.x * 4 + i];
    if (pred_dir != 2) {
        mv_l0 = &mv_info->mv[pred_dir];
        ref_idx_l0 = mv_info->ref_idx[pred_dir];
        RefPicList0 = get_ref_pic(mb, this->sets.slice->RefPicList[pred_dir], ref_idx_l0);
    } else {
        mv_l0 = &mv_info->mv[LIST_0];
        mv_l1 = &mv_info->mv[LIST_1];
        ref_idx_l0 = mv_info->ref_idx[LIST_0];
        ref_idx_l1 = mv_info->ref_idx[LIST_1];
        RefPicList0 = get_ref_pic(mb, this->sets.slice->RefPicList[LIST_0], ref_idx_l0);
        RefPicList1 = get_ref_pic(mb, this->sets.slice->RefPicList[LIST_1], ref_idx_l1);
    }

    int block_y_aff;
    if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag)
        block_y_aff = (mb.mb.y * 4 - 4 * (mb.mbAddrX % 2)) / 2;
    else
        block_y_aff = mb.mb.y * 4;

    this->check_motion_vector_range(mb, mv_l0);
    int vec1_x = (mb.mb.x * 4 + i) * 16 + mv_l0->mv_x;
    int vec1_y = (block_y_aff + j) * 16 + mv_l0->mv_y;
    int vec2_x, vec2_y;
    if (pred_dir == 2) {
        this->check_motion_vector_range(mb, mv_l1);
        vec2_x = (mb.mb.x * 4 + i) * 16 + mv_l1->mv_x;
        vec2_y = (block_y_aff + j) * 16 + mv_l1->mv_y;
    }

    // vars for get_block_luma
    px_t tmp_block_l0[16][16];
    px_t tmp_block_l1[16][16];
    px_t tmp_block_l2[16][16];
    px_t tmp_block_l3[16][16];

    if (CheckVertMV(&mb, vec1_y, block_size_y)) {
        get_block_luma(RefPicList0, vec1_x, vec1_y   , block_size_x, 8             , tmp_block_l0    , comp, mb);
        get_block_luma(RefPicList0, vec1_x, vec1_y+32, block_size_x, block_size_y-8, tmp_block_l0 + 8, comp, mb);
    } else
        get_block_luma(RefPicList0, vec1_x, vec1_y   , block_size_x, block_size_y  , tmp_block_l0    , comp, mb);
    if (pred_dir == 2) {
        if (CheckVertMV(&mb, vec2_y, block_size_y)) {
            get_block_luma(RefPicList1, vec2_x, vec2_y   , block_size_x, 8             , tmp_block_l1    , comp, mb);
            get_block_luma(RefPicList1, vec2_x, vec2_y+32, block_size_x, block_size_y-8, tmp_block_l1 + 8, comp, mb);
        } else
            get_block_luma(RefPicList1, vec2_x, vec2_y   , block_size_x, block_size_y  , tmp_block_l1    , comp, mb);
    }

    if (pred_dir != 2)
        mc_prediction(&this->sets.slice->mb_pred[comp][j * 4][i * 4], tmp_block_l0, block_size_y, block_size_x, mb, comp, ref_idx_l0, pred_dir);
    else
        bi_prediction(&this->sets.slice->mb_pred[comp][j * 4][i * 4], tmp_block_l0, tmp_block_l1, block_size_y, block_size_x, mb, comp, ref_idx_l0, ref_idx_l1);

    if (this->sets.sps->chroma_format_idc != CHROMA_FORMAT_400 && this->sets.sps->chroma_format_idc != CHROMA_FORMAT_444) {
        int maxold_x = this->sets.pic->size_x_cr - 1;
        int maxold_y = mb.mb_field_decoding_flag ? (this->sets.pic->size_y_cr >> 1) - 1 : this->sets.pic->size_y_cr - 1;

        int ioff_cr = this->sets.sps->MbWidthC  == 16 ? i * 4 : i * 2;
        int joff_cr = this->sets.sps->MbHeightC == 16 ? j * 4 : j * 2;
        int partWidthC  = this->sets.sps->MbWidthC  == 16 ? block_size_x : block_size_x >> 1;
        int partHeightC = this->sets.sps->MbHeightC == 16 ? block_size_y : block_size_y >> 1;

        get_block_chroma(RefPicList0, vec1_x, vec1_y, maxold_x, maxold_y, partWidthC, partHeightC, tmp_block_l0, PLANE_U, mb);
        get_block_chroma(RefPicList0, vec1_x, vec1_y, maxold_x, maxold_y, partWidthC, partHeightC, tmp_block_l2, PLANE_V, mb);
        if (pred_dir == 2) {
            get_block_chroma(RefPicList1, vec2_x, vec2_y, maxold_x, maxold_y, partWidthC, partHeightC, tmp_block_l1, PLANE_U, mb);
            get_block_chroma(RefPicList1, vec2_x, vec2_y, maxold_x, maxold_y, partWidthC, partHeightC, tmp_block_l3, PLANE_V, mb);
        }

        if (pred_dir != 2) {
            mc_prediction(&this->sets.slice->mb_pred[1][joff_cr][ioff_cr], tmp_block_l0, partHeightC, partWidthC, mb, PLANE_U, ref_idx_l0, pred_dir);
            mc_prediction(&this->sets.slice->mb_pred[2][joff_cr][ioff_cr], tmp_block_l2, partHeightC, partWidthC, mb, PLANE_V, ref_idx_l0, pred_dir);
        } else {
            bi_prediction(&this->sets.slice->mb_pred[1][joff_cr][ioff_cr], tmp_block_l0, tmp_block_l1, partHeightC, partWidthC, mb, PLANE_U, ref_idx_l0, ref_idx_l1);
            bi_prediction(&this->sets.slice->mb_pred[2][joff_cr][ioff_cr], tmp_block_l2, tmp_block_l3, partHeightC, partWidthC, mb, PLANE_V, ref_idx_l0, ref_idx_l1);
        }
    }
}


}
}
