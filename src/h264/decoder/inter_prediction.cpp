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
    storable_picture* curr_ref, int x_pos, int y_pos, int partWidthL, int partHeightL,
    px_t partPredLXL[16][16], int comp, mb_t& mb)
{
    int max_imgpel_value = (1 << (comp > 0 ? this->sets.sps->BitDepthC : this->sets.sps->BitDepthY)) - 1;

    if (curr_ref->no_ref) {
        memset(partPredLXL[0], max_imgpel_value, partHeightL * partWidthL * sizeof(px_t));
        return;
    }

    int mvLX[2] = {x_pos, y_pos};

    int maxold_x = this->sets.pic->size_x - 1;
    int maxold_y = (mb.mb_field_decoding_flag) ? (this->sets.pic->size_y >> 1) - 1 : this->sets.pic->size_y - 1;

    shr_t& shr = this->sets.slice->header;
    px_t** cur_imgY = (this->sets.sps->separate_colour_plane_flag && shr.colour_plane_id > PLANE_Y) ?
                        curr_ref->imgUV[shr.colour_plane_id-1] : 
                        comp ? curr_ref->imgUV[comp - 1] : curr_ref->imgY;

    uint32_t PicWidthInSamplesL  = this->sets.sps->PicWidthInSamplesL;
    uint32_t PicHeightInSamplesL = this->sets.slice->header.PicHeightInSamplesL;
    uint32_t refPicHeightEffectiveL =
        !this->sets.slice->header.MbaffFrameFlag || !mb.mb_field_decoding_flag ?
        PicHeightInSamplesL : PicHeightInSamplesL / 2;

    int xAL    = clip3(-18, maxold_x + 2, mvLX[0] >> 2);
    int yAL    = clip3(-10, maxold_y + 2, mvLX[1] >> 2);
    int xFracL = (mvLX[0] & 3);
    int yFracL = (mvLX[1] & 3);

    int tmp_res[16+5][16+5];

    if (xFracL == 0 && yFracL == 0) {
        for (int yL = 0; yL < partHeightL; yL++)
            memcpy(&partPredLXL[yL][0], &cur_imgY[yAL + yL][xAL], 16 * sizeof(px_t));
    } else if (xFracL == 0 || yFracL == 0) {
        for (int yL = 0; yL < partHeightL; yL++) {
            for (int xL = 0; xL < partWidthL; xL++) {
                int xIntL = xAL + xL;
                int yIntL = yAL + yL;

                if (yFracL == 0) {
                    int xEL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL - 2);
                    int xFL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL - 1);
                    int xGL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 0);
                    int xHL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 1);
                    int xIL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 2);
                    int xJL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 3);
                    int yZL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 0);
                    px_t E = cur_imgY[yZL][xEL];
                    px_t F = cur_imgY[yZL][xFL];
                    px_t G = cur_imgY[yZL][xGL];
                    px_t H = cur_imgY[yZL][xHL];
                    px_t I = cur_imgY[yZL][xIL];
                    px_t J = cur_imgY[yZL][xJL];
                    int b1 = (E - 5 * F + 20 * G + 20 * H - 5 * I + J);
                    int b  = clip1(max_imgpel_value, (b1 + 16) >> 5);

                    if (xFracL % 2 == 0)
                        partPredLXL[yL][xL] = (px_t) b;
                    else {
                        int G = cur_imgY[yZL][xFracL == 1 ? xGL : xHL];
                        int a = (G + b + 1) >> 1;
                        partPredLXL[yL][xL] = (px_t) a;
                    }
                }

                if (xFracL == 0) {
                    int xZL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 0);
                    int yAL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL - 2);
                    int yCL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL - 1);
                    int yGL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 0);
                    int yML = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 1);
                    int yRL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 2);
                    int yTL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 3);
                    px_t A = cur_imgY[yAL][xZL];
                    px_t C = cur_imgY[yCL][xZL];
                    px_t G = cur_imgY[yGL][xZL];
                    px_t M = cur_imgY[yML][xZL];
                    px_t R = cur_imgY[yRL][xZL];
                    px_t T = cur_imgY[yTL][xZL];
                    int h1 = (A - 5 * C + 20 * G + 20 * M - 5 * R + T);
                    int h  = clip1(max_imgpel_value, (h1 + 16) >> 5);

                    if (yFracL % 2 == 0)
                        partPredLXL[yL][xL] = (px_t) h;
                    else {
                        int G = cur_imgY[yFracL == 1 ? yGL : yML][xZL];
                        int d = (G + h + 1) >> 1;
                        partPredLXL[yL][xL] = (px_t) d;
                    }
                }
            }
        }
    } else if (xFracL % 2 == 1 && yFracL % 2 == 1) {
        for (int yL = 0; yL < partHeightL; yL++) {
            for (int xL = 0; xL < partWidthL; xL++) {
                int xIntL = xAL + xL;
                int yIntL = yAL + yL;

                int xEL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL - 2);
                int xFL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL - 1);
                int xGL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 0);
                int xHL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 1);
                int xIL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 2);
                int xJL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + 3);
                int yZL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + (yFracL == 3 ? 1 : 0));
                px_t E = cur_imgY[yZL][xEL];
                px_t F = cur_imgY[yZL][xFL];
                px_t G = cur_imgY[yZL][xGL];
                px_t H = cur_imgY[yZL][xHL];
                px_t I = cur_imgY[yZL][xIL];
                px_t J = cur_imgY[yZL][xJL];
                int b1 = (E - 5 * F + 20 * G + 20 * H - 5 * I + J);
                int b  = clip1(max_imgpel_value, (b1 + 16) >> 5);

                int xZL = clip3<int>(0, PicWidthInSamplesL - 1, xIntL + (xFracL == 3 ? 1 : 0));
                int yAL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL - 2);
                int yCL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL - 1);
                int yGL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 0);
                int yML = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 1);
                int yRL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 2);
                int yTL = clip3<int>(0, refPicHeightEffectiveL - 1, yIntL + 3);
                px_t A = cur_imgY[yAL][xZL];
                px_t C = cur_imgY[yCL][xZL];
                     G = cur_imgY[yGL][xZL];
                px_t M = cur_imgY[yML][xZL];
                px_t R = cur_imgY[yRL][xZL];
                px_t T = cur_imgY[yTL][xZL];
                int h1 = (A - 5 * C + 20 * G + 20 * M - 5 * R + T);
                int h  = clip1(max_imgpel_value, (h1 + 16) >> 5);

                int e  = (b + h + 1) >> 1;
                partPredLXL[yL][xL] = (px_t) e;
            }
        }
    } else {
        int shift_x = xFracL == 2 ? 1 : this->sets.pic->iLumaStride;
        int step_x  = xFracL != 2 && yFracL == 2 ? 1 : 0;
        int step_y  = xFracL == 2 || yFracL != 2 ? 1 : 0;

        int img_pos_x = xAL - 2;
        int img_pos_y = yAL - 2;
        int blk_pos_x = yFracL == 2 ? (xFracL == 3 ? 3 : 2) : 0;
        int blk_pos_y = xFracL == 2 ? (yFracL == 3 ? 3 : 2) : 0;

        for (int yL = 0; yL < partHeightL + (xFracL == 2 ? 5 : 0); yL++) {
            for (int xL = 0; xL < partWidthL + (yFracL == 2 ? 5 : 0); xL++) {
                int xIntL = img_pos_x + xL;
                int yIntL = img_pos_y + yL;
                px_t A = cur_imgY[yIntL][xIntL + shift_x * 0];
                px_t C = cur_imgY[yIntL][xIntL + shift_x * 1];
                px_t G = cur_imgY[yIntL][xIntL + shift_x * 2];
                px_t M = cur_imgY[yIntL][xIntL + shift_x * 3];
                px_t R = cur_imgY[yIntL][xIntL + shift_x * 4];
                px_t T = cur_imgY[yIntL][xIntL + shift_x * 5];
                int h1 = (A - 5 * C + 20 * G + 20 * M - 5 * R + T);
                tmp_res[yL][xL] = h1;
            }
        }

        for (int yL = 0; yL < partHeightL; yL++) {
            for (int xL = 0; xL < partWidthL; xL++) {
                int aa = tmp_res[yL + step_y * 0][xL + step_x * 0];
                int bb = tmp_res[yL + step_y * 1][xL + step_x * 1];
                int b1 = tmp_res[yL + step_y * 2][xL + step_x * 2];
                int s1 = tmp_res[yL + step_y * 3][xL + step_x * 3];
                int gg = tmp_res[yL + step_y * 4][xL + step_x * 4];
                int hh = tmp_res[yL + step_y * 5][xL + step_x * 5];
                int j1 = (aa - 5 * bb + 20 * b1 + 20 * s1 - 5 * gg + hh);
                int j  = clip1(max_imgpel_value, (j1 + 512) >> 10);
                partPredLXL[yL][xL] = (px_t) j;

                if (xFracL == 1 || xFracL == 3 || yFracL == 1 || yFracL == 3) {
                    px_t q0 = clip1(max_imgpel_value, (tmp_res[blk_pos_y + yL][blk_pos_x + xL] + 16) >> 5);
                    partPredLXL[yL][xL] = (px_t) ((partPredLXL[yL][xL] + q0 + 1) >> 1);
                }
            }
        }
    }
}

void InterPrediction::get_block_chroma(storable_picture* pic,
    int mvLXX[2], int partWidthC, int partHeightC,
    px_t predPartLXC[16][16], int comp, mb_t& mb)
{
    int mvLX[2] = {mvLXX[0], mvLXX[1]};

    if (this->sets.sps->ChromaArrayType == 1) {
        slice_t& slice = *this->sets.slice;
        shr_t& shr = slice.header;

        if (!shr.MbaffFrameFlag && shr.field_pic_flag) {
            if (shr.structure != pic->slice.structure)
                mvLX[1] += (!shr.bottom_field_flag ? -2 : 2);
        }
        if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag) {
            if (mb.mbAddrX % 2 == 0 && pic->slice.structure == BOTTOM_FIELD)
                mvLX[1] -= 2;
            if (mb.mbAddrX % 2 == 1 && pic->slice.structure == TOP_FIELD)
                mvLX[1] += 2;
        }
    }

    int mvCX[2] = {mvLX[0], mvLX[1]};

    if (pic->no_ref) {
        px_t no_ref_value = (px_t)(1 << (this->sets.sps->BitDepthC - 1));
        memset(predPartLXC, no_ref_value, partHeightC * partWidthC * sizeof(px_t));
        return;
    }

    px_t** imgUV = pic->imgUV[comp - 1];

    uint32_t PicWidthInSamplesC  = this->sets.sps->PicWidthInSamplesC;
    uint32_t PicHeightInSamplesC = this->sets.slice->header.PicHeightInSamplesC;
    uint32_t refPicHeightEffectiveC =
        !this->sets.slice->header.MbaffFrameFlag || !mb.mb_field_decoding_flag ?
        PicHeightInSamplesC : PicHeightInSamplesC / 2;

    int xAL    = (mvCX[0] >> 3);
    int yAL    = (this->sets.sps->ChromaArrayType == 1) ? (mvCX[1] >> 3) : (mvCX[1] >> 2);
    int xFracC = (mvCX[0] & 7);
    int yFracC = (this->sets.sps->ChromaArrayType == 1) ? (mvCX[1] & 7) : (mvCX[1] & 3) << 1;

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
            predPartLXC[yC][xC] =
                ((8 - xFracC) * (8 - yFracC) * A + xFracC * (8 - yFracC) * B +
                 (8 - xFracC) * yFracC * C + xFracC * yFracC * D + 32) >> 6;
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

void InterPrediction::inter_pred(mb_t& mb, int comp, int pred_dir, int i, int j, int partWidthL, int partHeightL)
{
    shr_t& shr = this->sets.slice->header;

    mv_t* mv_l0, *mv_l1;
    short ref_idx_l0, ref_idx_l1;
    storable_picture* refPic0, *refPic1;

    auto mv_info = &this->sets.pic->mv_info[mb.mb.y * 4 + j][mb.mb.x * 4 + i];
    if (pred_dir != 2) {
        mv_l0 = &mv_info->mv[pred_dir];
        ref_idx_l0 = mv_info->ref_idx[pred_dir];
        refPic0 = get_ref_pic(mb, this->sets.slice->RefPicList[pred_dir], ref_idx_l0);
    } else {
        mv_l0 = &mv_info->mv[LIST_0];
        mv_l1 = &mv_info->mv[LIST_1];
        ref_idx_l0 = mv_info->ref_idx[LIST_0];
        ref_idx_l1 = mv_info->ref_idx[LIST_1];
        refPic0 = get_ref_pic(mb, this->sets.slice->RefPicList[LIST_0], ref_idx_l0);
        refPic1 = get_ref_pic(mb, this->sets.slice->RefPicList[LIST_1], ref_idx_l1);
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
    px_t partPredL0L[16][16];
    px_t partPredL1L[16][16];

    if (CheckVertMV(&mb, vec1_y, partHeightL)) {
        get_block_luma(refPic0, vec1_x, vec1_y   , partWidthL, 8            , partPredL0L  , comp, mb);
        get_block_luma(refPic0, vec1_x, vec1_y+32, partWidthL, partHeightL-8, partPredL0L+8, comp, mb);
    } else
        get_block_luma(refPic0, vec1_x, vec1_y   , partWidthL, partHeightL  , partPredL0L  , comp, mb);
    if (pred_dir == 2) {
        if (CheckVertMV(&mb, vec2_y, partHeightL)) {
            get_block_luma(refPic1, vec2_x, vec2_y   , partWidthL, 8            , partPredL1L  , comp, mb);
            get_block_luma(refPic1, vec2_x, vec2_y+32, partWidthL, partHeightL-8, partPredL1L+8, comp, mb);
        } else
            get_block_luma(refPic1, vec2_x, vec2_y   , partWidthL, partHeightL  , partPredL1L  , comp, mb);
    }

    if (pred_dir != 2)
        mc_prediction(&this->sets.slice->mb_pred[comp][j * 4][i * 4], partPredL0L, partHeightL, partWidthL, mb, comp, ref_idx_l0, pred_dir);
    else
        bi_prediction(&this->sets.slice->mb_pred[comp][j * 4][i * 4], partPredL0L, partPredL1L, partHeightL, partWidthL, mb, comp, ref_idx_l0, ref_idx_l1);

    if (this->sets.sps->ChromaArrayType == 1 || this->sets.sps->ChromaArrayType == 2) {
        px_t partPredL0Cb[16][16];
        px_t partPredL1Cb[16][16];
        px_t partPredL0Cr[16][16];
        px_t partPredL1Cr[16][16];

        int mvL0[2] = {vec1_x, vec1_y};
        int mvL1[2] = {vec2_x, vec2_y};

        int partWidthC  = partWidthL  / this->sets.sps->SubWidthC;
        int partHeightC = partHeightL / this->sets.sps->SubHeightC;
        get_block_chroma(refPic0, mvL0, partWidthC, partHeightC, partPredL0Cb, PLANE_U, mb);
        get_block_chroma(refPic0, mvL0, partWidthC, partHeightC, partPredL0Cr, PLANE_V, mb);
        if (pred_dir == 2) {
            get_block_chroma(refPic1, mvL1, partWidthC, partHeightC, partPredL1Cb, PLANE_U, mb);
            get_block_chroma(refPic1, mvL1, partWidthC, partHeightC, partPredL1Cr, PLANE_V, mb);
        }

        int xO = i * 4 / this->sets.sps->SubWidthC;
        int yO = j * 4 / this->sets.sps->SubHeightC;
        if (pred_dir != 2) {
            mc_prediction(&this->sets.slice->mb_pred[1][yO][xO], partPredL0Cb, partHeightC, partWidthC, mb, PLANE_U, ref_idx_l0, pred_dir);
            mc_prediction(&this->sets.slice->mb_pred[2][yO][xO], partPredL0Cr, partHeightC, partWidthC, mb, PLANE_V, ref_idx_l0, pred_dir);
        } else {
            bi_prediction(&this->sets.slice->mb_pred[1][yO][xO], partPredL0Cb, partPredL1Cb, partHeightC, partWidthC, mb, PLANE_U, ref_idx_l0, ref_idx_l1);
            bi_prediction(&this->sets.slice->mb_pred[2][yO][xO], partPredL0Cr, partPredL1Cr, partHeightC, partWidthC, mb, PLANE_V, ref_idx_l0, ref_idx_l1);
        }
    }
}


}
}
