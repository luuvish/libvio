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
 *  File      : transform.cpp
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
#include "transform.h"
#include "image.h"

#include "intra_prediction.h"


namespace vio  {
namespace h264 {


static int Flat_4x4_16[16] = {
    16, 16, 16, 16,
    16, 16, 16, 16,
    16, 16, 16, 16,
    16, 16, 16, 16
};

static int Flat_8x8_16[64] = {
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16
};

// Table 7-3 Specification of default scaling lists Default_4x4_Intra and Default_4x4_Inter

static int Default_4x4_Intra[16] = {
     6, 13, 20, 28,
    13, 20, 28, 32,
    20, 28, 32, 37,
    28, 32, 37, 42
};

static int Default_4x4_Inter[16] = {
    10, 14, 20, 24,
    14, 20, 24, 27,
    20, 24, 27, 30,
    24, 27, 30, 34
};

// Table 7-4 Specification of default scaling lists Default_8x8_Intra and Default_8x8_Inter

static int Default_8x8_Intra[64] = {
     6, 10, 13, 16, 18, 23, 25, 27,
    10, 11, 16, 18, 23, 25, 27, 29,
    13, 16, 18, 23, 25, 27, 29, 31,
    16, 18, 23, 25, 27, 29, 31, 33,
    18, 23, 25, 27, 29, 31, 33, 36,
    23, 25, 27, 29, 31, 33, 36, 38,
    25, 27, 29, 31, 33, 36, 38, 40,
    27, 29, 31, 33, 36, 38, 40, 42
};

static int Default_8x8_Inter[64] = {
     9, 13, 15, 17, 19, 21, 22, 24,
    13, 13, 17, 19, 21, 22, 24, 25,
    15, 17, 19, 21, 22, 24, 25, 27,
    17, 19, 21, 22, 24, 25, 27, 28,
    19, 21, 22, 24, 25, 27, 28, 30,
    21, 22, 24, 25, 27, 28, 30, 32,
    22, 24, 25, 27, 28, 30, 32, 33,
    24, 25, 27, 28, 30, 32, 33, 35
};


static const int dequant_coef[6][4][4] = {
    {{ 10, 13, 10, 13 },
     { 13, 16, 13, 16 },
     { 10, 13, 10, 13 },
     { 13, 16, 13, 16 }},
    {{ 11, 14, 11, 14 },
     { 14, 18, 14, 18 },
     { 11, 14, 11, 14 },
     { 14, 18, 14, 18 }},
    {{ 13, 16, 13, 16 },
     { 16, 20, 16, 20 },
     { 13, 16, 13, 16 },
     { 16, 20, 16, 20 }},
    {{ 14, 18, 14, 18 },
     { 18, 23, 18, 23 },
     { 14, 18, 14, 18 },
     { 18, 23, 18, 23 }},
    {{ 16, 20, 16, 20 },
     { 20, 25, 20, 25 },
     { 16, 20, 16, 20 },
     { 20, 25, 20, 25 }},
    {{ 18, 23, 18, 23 },
     { 23, 29, 23, 29 },
     { 18, 23, 18, 23 },
     { 23, 29, 23, 29 }}
};

// exported variables
static const int dequant_coef8[6][8][8] = {
    {{ 20, 19, 25, 19, 20, 19, 25, 19 },
     { 19, 18, 24, 18, 19, 18, 24, 18 },
     { 25, 24, 32, 24, 25, 24, 32, 24 },
     { 19, 18, 24, 18, 19, 18, 24, 18 },
     { 20, 19, 25, 19, 20, 19, 25, 19 },
     { 19, 18, 24, 18, 19, 18, 24, 18 },
     { 25, 24, 32, 24, 25, 24, 32, 24 },
     { 19, 18, 24, 18, 19, 18, 24, 18 }},
    {{ 22, 21, 28, 21, 22, 21, 28, 21 },
     { 21, 19, 26, 19, 21, 19, 26, 19 },
     { 28, 26, 35, 26, 28, 26, 35, 26 },
     { 21, 19, 26, 19, 21, 19, 26, 19 },
     { 22, 21, 28, 21, 22, 21, 28, 21 },
     { 21, 19, 26, 19, 21, 19, 26, 19 },
     { 28, 26, 35, 26, 28, 26, 35, 26 },
     { 21, 19, 26, 19, 21, 19, 26, 19 }},
    {{ 26, 24, 33, 24, 26, 24, 33, 24 },
     { 24, 23, 31, 23, 24, 23, 31, 23 },
     { 33, 31, 42, 31, 33, 31, 42, 31 },
     { 24, 23, 31, 23, 24, 23, 31, 23 },
     { 26, 24, 33, 24, 26, 24, 33, 24 },
     { 24, 23, 31, 23, 24, 23, 31, 23 },
     { 33, 31, 42, 31, 33, 31, 42, 31 },
     { 24, 23, 31, 23, 24, 23, 31, 23 }},
    {{ 28, 26, 35, 26, 28, 26, 35, 26 },
     { 26, 25, 33, 25, 26, 25, 33, 25 },
     { 35, 33, 45, 33, 35, 33, 45, 33 },
     { 26, 25, 33, 25, 26, 25, 33, 25 },
     { 28, 26, 35, 26, 28, 26, 35, 26 },
     { 26, 25, 33, 25, 26, 25, 33, 25 },
     { 35, 33, 45, 33, 35, 33, 45, 33 },
     { 26, 25, 33, 25, 26, 25, 33, 25 }},
    {{ 32, 30, 40, 30, 32, 30, 40, 30 },
     { 30, 28, 38, 28, 30, 28, 38, 28 },
     { 40, 38, 51, 38, 40, 38, 51, 38 },
     { 30, 28, 38, 28, 30, 28, 38, 28 },
     { 32, 30, 40, 30, 32, 30, 40, 30 },
     { 30, 28, 38, 28, 30, 28, 38, 28 },
     { 40, 38, 51, 38, 40, 38, 51, 38 },
     { 30, 28, 38, 28, 30, 28, 38, 28 }},
    {{ 36, 34, 46, 34, 36, 34, 46, 34 },
     { 34, 32, 43, 32, 34, 32, 43, 32 },
     { 46, 43, 58, 43, 46, 43, 58, 43 },
     { 34, 32, 43, 32, 34, 32, 43, 32 },
     { 36, 34, 46, 34, 36, 34, 46, 34 },
     { 34, 32, 43, 32, 34, 32, 43, 32 },
     { 46, 43, 58, 43, 46, 43, 58, 43 },
     { 34, 32, 43, 32, 34, 32, 43, 32 }}
};


void Transform::init(slice_t& slice)
{
    sps_t* sps = slice.active_sps;
    pps_t* pps = slice.active_pps;

    if (!pps->pic_scaling_matrix_present_flag &&
        !sps->seq_scaling_matrix_present_flag) {
        for (int i = 0; i < 12; i++)
            this->qmatrix[i] = (i < 6) ? Flat_4x4_16 : Flat_8x8_16;
    } else {
        int n_ScalingList = (sps->chroma_format_idc != YUV444) ? 8 : 12;
        if (sps->seq_scaling_matrix_present_flag) { // check sps first
            for (int i = 0; i < n_ScalingList; i++) {
                if (i < 6) {
                    if (!sps->seq_scaling_list_present_flag[i]) { // fall-back rule A
                        if (i == 0)
                            this->qmatrix[i] = Default_4x4_Intra;
                        else if (i == 3)
                            this->qmatrix[i] = Default_4x4_Inter;
                        else
                            this->qmatrix[i] = this->qmatrix[i-1];
                    } else {
                        if (sps->UseDefaultScalingMatrix4x4Flag[i])
                            this->qmatrix[i] = (i < 3) ? Default_4x4_Intra : Default_4x4_Inter;
                        else
                            this->qmatrix[i] = sps->ScalingList4x4[i];
                    }
                } else {
                    if (!sps->seq_scaling_list_present_flag[i]) { // fall-back rule A
                        if (i == 6)
                            this->qmatrix[i] = Default_8x8_Intra;
                        else if (i == 7)
                            this->qmatrix[i] = Default_8x8_Inter;
                        else
                            this->qmatrix[i] = this->qmatrix[i-2];
                    } else {
                        if (sps->UseDefaultScalingMatrix8x8Flag[i-6])
                            this->qmatrix[i] = (i == 6 || i == 8 || i == 10) ? Default_8x8_Intra : Default_8x8_Inter;
                        else
                            this->qmatrix[i] = sps->ScalingList8x8[i-6];
                    }
                }
            }
        }

        if (pps->pic_scaling_matrix_present_flag) { // then check pps
            for (int i = 0; i < n_ScalingList; i++) {
                if (i < 6) {
                    if (!pps->pic_scaling_list_present_flag[i]) { // fall-back rule B
                        if (i == 0) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                this->qmatrix[i] = Default_4x4_Intra;
                        } else if (i == 3) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                this->qmatrix[i] = Default_4x4_Inter;
                        } else
                            this->qmatrix[i] = this->qmatrix[i-1];
                    } else {
                        if (pps->UseDefaultScalingMatrix4x4Flag[i])
                            this->qmatrix[i] = (i < 3) ? Default_4x4_Intra : Default_4x4_Inter;
                        else
                            this->qmatrix[i] = pps->ScalingList4x4[i];
                    }
                } else {
                    if (!pps->pic_scaling_list_present_flag[i]) { // fall-back rule B
                        if (i == 6) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                this->qmatrix[i] = Default_8x8_Intra;
                        } else if (i == 7) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                this->qmatrix[i] = Default_8x8_Inter;
                        } else  
                            this->qmatrix[i] = this->qmatrix[i-2];
                    } else {
                        if (pps->UseDefaultScalingMatrix8x8Flag[i-6])
                            this->qmatrix[i] = (i == 6 || i == 8 || i == 10) ? Default_8x8_Intra : Default_8x8_Inter;
                        else
                            this->qmatrix[i] = pps->ScalingList8x8[i-6];
                    }
                }
            }
        }
    }

    this->set_quant(slice);
}

void Transform::set_quant(slice_t& slice)
{
    sps_t* sps = slice.active_sps;
    pps_t* pps = slice.active_pps;

    for (int k = 0; k < 6; k++) {
        for (int j = 0; j < 4; j++) {
            for (int i = 0; i < 4; i++) {
                this->InvLevelScale4x4_Intra[0][k][j][i] = dequant_coef[k][j][i] * this->qmatrix[0][4 * j + i];
                this->InvLevelScale4x4_Intra[1][k][j][i] = dequant_coef[k][j][i] * this->qmatrix[1][4 * j + i];
                this->InvLevelScale4x4_Intra[2][k][j][i] = dequant_coef[k][j][i] * this->qmatrix[2][4 * j + i];
                this->InvLevelScale4x4_Inter[0][k][j][i] = dequant_coef[k][j][i] * this->qmatrix[3][4 * j + i];
                this->InvLevelScale4x4_Inter[1][k][j][i] = dequant_coef[k][j][i] * this->qmatrix[4][4 * j + i];
                this->InvLevelScale4x4_Inter[2][k][j][i] = dequant_coef[k][j][i] * this->qmatrix[5][4 * j + i];
            }
        }
    }

    if (!pps->transform_8x8_mode_flag)
        return;

    for (int k = 0; k < 6; k++) {
        for (int j = 0; j < 8; j++) {
            for (int i = 0; i < 8; i++) {
                this->InvLevelScale8x8_Intra[0][k][j][i] = dequant_coef8[k][j][i] * this->qmatrix[6][8 * j + i];
                this->InvLevelScale8x8_Inter[0][k][j][i] = dequant_coef8[k][j][i] * this->qmatrix[7][8 * j + i];
            }
        }
    }

    if (sps->chroma_format_idc != 3)
        return;

    for (int k = 0; k < 6; k++) {
        for (int j = 0; j < 8; j++) {
            for (int i = 0; i < 8; i++) {
                this->InvLevelScale8x8_Intra[1][k][j][i] = dequant_coef8[k][j][i] * this->qmatrix[ 8][8 * j + i];
                this->InvLevelScale8x8_Inter[1][k][j][i] = dequant_coef8[k][j][i] * this->qmatrix[ 9][8 * j + i];
                this->InvLevelScale8x8_Intra[2][k][j][i] = dequant_coef8[k][j][i] * this->qmatrix[10][8 * j + i];
                this->InvLevelScale8x8_Inter[2][k][j][i] = dequant_coef8[k][j][i] * this->qmatrix[11][8 * j + i];
            }
        }
    }
}


// Table 8-13 Specification of mapping of idx to Cij for zig-zag and field scan

static const uint8_t ZIGZAG_SCAN_4x4[2][16][2] = {
    {{ 0, 0 }, { 1, 0 }, { 0, 1 }, { 0, 2 },
     { 1, 1 }, { 2, 0 }, { 3, 0 }, { 2, 1 },
     { 1, 2 }, { 0, 3 }, { 1, 3 }, { 2, 2 },
     { 3, 1 }, { 3, 2 }, { 2, 3 }, { 3, 3 }},
    {{ 0, 0 }, { 0, 1 }, { 1, 0 }, { 0, 2 },
     { 0, 3 }, { 1, 1 }, { 1, 2 }, { 1, 3 },
     { 2, 0 }, { 2, 1 }, { 2, 2 }, { 2, 3 },
     { 3, 0 }, { 3, 1 }, { 3, 2 }, { 3, 3 }}
};

// Table 8-14 Specification of mapping of idx to Cij for 8x8 zig-zag and 8x8 field scan

static const uint8_t ZIGZAG_SCAN_8x8[2][64][2] = {
    {{ 0, 0 }, { 1, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 2, 0 }, { 3, 0 }, { 2, 1 },
     { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 3 }, { 2, 2 }, { 3, 1 }, { 4, 0 }, { 5, 0 },
     { 4, 1 }, { 3, 2 }, { 2, 3 }, { 1, 4 }, { 0, 5 }, { 0, 6 }, { 1, 5 }, { 2, 4 },
     { 3, 3 }, { 4, 2 }, { 5, 1 }, { 6, 0 }, { 7, 0 }, { 6, 1 }, { 5, 2 }, { 4, 3 },
     { 3, 4 }, { 2, 5 }, { 1, 6 }, { 0, 7 }, { 1, 7 }, { 2, 6 }, { 3, 5 }, { 4, 4 },
     { 5, 3 }, { 6, 2 }, { 7, 1 }, { 7, 2 }, { 6, 3 }, { 5, 4 }, { 4, 5 }, { 3, 6 },
     { 2, 7 }, { 3, 7 }, { 4, 6 }, { 5, 5 }, { 6, 4 }, { 7, 3 }, { 7, 4 }, { 6, 5 },
     { 5, 6 }, { 4, 7 }, { 5, 7 }, { 6, 6 }, { 7, 5 }, { 7, 6 }, { 6, 7 }, { 7, 7 }},
    {{ 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 }, { 0, 3 }, { 0, 4 }, { 1, 2 },
     { 2, 0 }, { 1, 3 }, { 0, 5 }, { 0, 6 }, { 0, 7 }, { 1, 4 }, { 2, 1 }, { 3, 0 },
     { 2, 2 }, { 1, 5 }, { 1, 6 }, { 1, 7 }, { 2, 3 }, { 3, 1 }, { 4, 0 }, { 3, 2 },
     { 2, 4 }, { 2, 5 }, { 2, 6 }, { 2, 7 }, { 3, 3 }, { 4, 1 }, { 5, 0 }, { 4, 2 },
     { 3, 4 }, { 3, 5 }, { 3, 6 }, { 3, 7 }, { 4, 3 }, { 5, 1 }, { 6, 0 }, { 5, 2 },
     { 4, 4 }, { 4, 5 }, { 4, 6 }, { 4, 7 }, { 5, 3 }, { 6, 1 }, { 6, 2 }, { 5, 4 },
     { 5, 5 }, { 5, 6 }, { 5, 7 }, { 6, 3 }, { 7, 0 }, { 7, 1 }, { 6, 4 }, { 6, 5 },
     { 6, 6 }, { 6, 7 }, { 7, 2 }, { 7, 3 }, { 7, 4 }, { 7, 5 }, { 7, 6 }, { 7, 7 }}
};

pos_t Transform::inverse_scan_luma_dc(mb_t* mb, int run)
{
    slice_t* slice = mb->p_Slice;

    bool field = slice->field_pic_flag || mb->mb_field_decoding_flag;
    const uint8_t (*zigzag_scan_4x4)[2] = ZIGZAG_SCAN_4x4[field];

    return {zigzag_scan_4x4[run][0], zigzag_scan_4x4[run][1]};
}

pos_t Transform::inverse_scan_luma_ac(mb_t* mb, int run)
{
    slice_t* slice = mb->p_Slice;

    bool field = slice->field_pic_flag || mb->mb_field_decoding_flag;
    const uint8_t (*zigzag_scan_4x4)[2] = ZIGZAG_SCAN_4x4[field];
    const uint8_t (*zigzag_scan_8x8)[2] = ZIGZAG_SCAN_8x8[field];

    if (!mb->transform_size_8x8_flag)
        return {zigzag_scan_4x4[run][0], zigzag_scan_4x4[run][1]};
    else
        return {zigzag_scan_8x8[run][0], zigzag_scan_8x8[run][1]};
}

pos_t Transform::inverse_scan_chroma_dc(mb_t* mb, int run)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    if (sps->ChromaArrayType == 1)
        return {run % 2, run / 2};
    if (sps->ChromaArrayType == 2)
        return {ZIGZAG_SCAN_4x4[1][run][0], ZIGZAG_SCAN_4x4[1][run][1]};
    return {0, 0};
}

pos_t Transform::inverse_scan_chroma_ac(mb_t* mb, int run)
{
    slice_t* slice = mb->p_Slice;

    bool field = slice->field_pic_flag || mb->mb_field_decoding_flag;
    const uint8_t (*zigzag_scan_4x4)[2] = ZIGZAG_SCAN_4x4[field];

    return {zigzag_scan_4x4[run][0], zigzag_scan_4x4[run][1]};
}


static inline int rshift_rnd_sf(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}

int Transform::inverse_quantize(mb_t* mb, bool uv, ColorPlane pl, int i0, int j0, int levarr)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int qp_per = mb->qp_scaled[pl] / 6;
    int qp_rem = mb->qp_scaled[pl] % 6;
    int transform_pl = sps->separate_colour_plane_flag ? slice->colour_plane_id : pl;

    if (uv) {
        int (*InvLevelScale4x4)[4] = mb->is_intra_block ?
            this->InvLevelScale4x4_Intra[pl][qp_rem] :
            this->InvLevelScale4x4_Inter[pl][qp_rem];
        levarr = rshift_rnd_sf((levarr * InvLevelScale4x4[j0][i0]) << qp_per, 4);
    } else if (!mb->transform_size_8x8_flag) {
        int (*InvLevelScale4x4)[4] = mb->is_intra_block ?
            this->InvLevelScale4x4_Intra[transform_pl][qp_rem] :
            this->InvLevelScale4x4_Inter[transform_pl][qp_rem];
        levarr = rshift_rnd_sf((levarr * InvLevelScale4x4[j0][i0]) << qp_per, 4);
    } else {
        int (*InvLevelScale8x8)[8] = mb->is_intra_block ?
            this->InvLevelScale8x8_Intra[transform_pl][qp_rem] :
            this->InvLevelScale8x8_Inter[transform_pl][qp_rem];
        levarr = rshift_rnd_sf((levarr * InvLevelScale8x8[j0][i0]) << qp_per, 6);
    }

    return levarr;
}


void Transform::coeff_luma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    const pos_t& pos = inverse_scan_luma_dc(mb, runarr);
    this->cof[pl][pos.y * 4][pos.x * 4] = levarr;
}

void Transform::coeff_luma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    if (!mb->transform_size_8x8_flag)
        mb->cbp_blks[pl] |= ((uint64_t)0x01 << (y0 * 4 + x0));
    else
        mb->cbp_blks[pl] |= ((uint64_t)0x33 << (y0 * 4 + x0));

    const pos_t& pos = inverse_scan_luma_ac(mb, runarr);
    if (!mb->TransformBypassModeFlag)
        levarr = this->inverse_quantize(mb, false, pl, pos.x, pos.y, levarr);
    this->cof[pl][y0 * 4 + pos.y][x0 * 4 + pos.x] = levarr;
}

void Transform::coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    const pos_t& pos = inverse_scan_chroma_dc(mb, runarr);
    this->cof[pl][pos.y * 4][pos.x * 4] = levarr;
}

void Transform::coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    const pos_t& pos = inverse_scan_chroma_ac(mb, runarr);
    if (!mb->TransformBypassModeFlag)
        levarr = this->inverse_quantize(mb, true, pl, pos.x, pos.y, levarr);
    this->cof[pl][y0 * 4 + pos.y][x0 * 4 + pos.x] = levarr;
}



void Transform::ihadamard_2x2(int c[2][2], int f[2][2])
{
    int e[2][2];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 2; ++i) {
        int ci0 = c[i][0];
        int ci1 = c[i][1];

        e[i][0] = ci0 + ci1;
        e[i][1] = ci0 - ci1;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 2; ++j) {
        int g0j = e[0][j];
        int g1j = e[1][j];
        
        f[0][j] = g0j + g1j;
        f[1][j] = g0j - g1j;
    }
}

void Transform::ihadamard_2x4(int c[4][2], int f[4][2])
{
    int e[4][2];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 4; ++i) {
        int ci0 = c[i][0];
        int ci1 = c[i][1];

        e[i][0] = ci0 + ci1;
        e[i][1] = ci0 - ci1;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 2; ++j) {
        int g0j = e[0][j];
        int g1j = e[1][j];
        int g3j = e[2][j];
        int g4j = e[3][j];

        int h0j = g0j + g3j;
        int h1j = g0j - g3j;
        int h2j = g1j - g4j;
        int h3j = g1j + g4j;
        
        f[0][j] = h0j + h3j;
        f[1][j] = h1j + h2j;
        f[2][j] = h1j - h2j;
        f[3][j] = h0j - h3j;
    }
}

void Transform::ihadamard_4x4(int c[4][4], int f[4][4])
{
    int e[4][4];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 4; ++i) {
        int ci0 = c[i][0];
        int ci1 = c[i][1];
        int ci2 = c[i][2];
        int ci3 = c[i][3];

        int di0 = ci0 + ci2;
        int di1 = ci0 - ci2;
        int di2 = ci1 - ci3;
        int di3 = ci1 + ci3;

        e[i][0] = di0 + di3;
        e[i][1] = di1 + di2;
        e[i][2] = di1 - di2;
        e[i][3] = di0 - di3;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 4; ++j) {
        int g0j = e[0][j];
        int g1j = e[1][j];
        int g2j = e[2][j];
        int g3j = e[3][j];

        int h0j = g0j + g2j;
        int h1j = g0j - g2j;
        int h2j = g1j - g3j;
        int h3j = g1j + g3j;

        f[0][j] = h0j + h3j;
        f[1][j] = h1j + h2j;
        f[2][j] = h1j - h2j;
        f[3][j] = h0j - h3j;
    }
}

void Transform::forward_4x4(int p[16][16], int c[16][16], int pos_y, int pos_x)
{
    int f[4][4];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 4; ++i) {
        int pi0 = p[pos_y + i][pos_x + 0];
        int pi1 = p[pos_y + i][pos_x + 1];
        int pi2 = p[pos_y + i][pos_x + 2];
        int pi3 = p[pos_y + i][pos_x + 3];

        int ei0 = pi0 + pi3;
        int ei1 = pi1 + pi2;
        int ei2 = pi1 - pi2;
        int ei3 = pi0 - pi3;

        f[i][0] = ei0 + ei1;
        f[i][1] = ei2 + (ei3 << 1);
        f[i][2] = ei0 - ei1;
        f[i][3] = ei3 - (ei2 << 1);
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 4; ++j) {
        int f0j = f[0][j];
        int f1j = f[1][j];
        int f2j = f[2][j];
        int f3j = f[3][j];

        int g0j = f0j + f3j;
        int g1j = f1j + f2j;
        int g2j = f1j - f2j;
        int g3j = f0j - f3j;

        c[pos_y + 0][pos_x + j] = g0j + g1j;
        c[pos_y + 1][pos_x + j] = g2j + (g3j << 1);
        c[pos_y + 2][pos_x + j] = g0j - g1j;
        c[pos_y + 3][pos_x + j] = g3j - (g2j << 1);
    }
}

void Transform::inverse_4x4(int d[16][16], int r[16][16], int pos_y, int pos_x)
{
    int f[4][4];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 4; ++i) {
        int di0 = d[pos_y + i][pos_x + 0];
        int di1 = d[pos_y + i][pos_x + 1];
        int di2 = d[pos_y + i][pos_x + 2];
        int di3 = d[pos_y + i][pos_x + 3];

        int ei0 = di0 + di2;
        int ei1 = di0 - di2;
        int ei2 = (di1 >> 1) - di3;
        int ei3 = di1 + (di3 >> 1);

        f[i][0] = ei0 + ei3;
        f[i][1] = ei1 + ei2;
        f[i][2] = ei1 - ei2;
        f[i][3] = ei0 - ei3;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 4; ++j) {
        int f0j = f[0][j];
        int f1j = f[1][j];
        int f2j = f[2][j];
        int f3j = f[3][j];

        int g0j = f0j + f2j;
        int g1j = f0j - f2j;
        int g2j = (f1j >> 1) - f3j;
        int g3j = f1j + (f3j >> 1);

        int h0j = g0j + g3j;
        int h1j = g1j + g2j;
        int h2j = g1j - g2j;
        int h3j = g0j - g3j;

        r[pos_y + 0][pos_x + j] = (h0j + (1 << 5)) >> 6;
        r[pos_y + 1][pos_x + j] = (h1j + (1 << 5)) >> 6;
        r[pos_y + 2][pos_x + j] = (h2j + (1 << 5)) >> 6;
        r[pos_y + 3][pos_x + j] = (h3j + (1 << 5)) >> 6;
    }
}

void Transform::inverse_8x8(int d[16][16], int r[16][16], int pos_y, int pos_x)
{
    int g[8][8];

    // horizontal row of scaled transform coefficients
    for (int i = 0; i < 8; ++i) {
        int di0 = d[pos_y + i][pos_x + 0];
        int di1 = d[pos_y + i][pos_x + 1];
        int di2 = d[pos_y + i][pos_x + 2];
        int di3 = d[pos_y + i][pos_x + 3];
        int di4 = d[pos_y + i][pos_x + 4];
        int di5 = d[pos_y + i][pos_x + 5];
        int di6 = d[pos_y + i][pos_x + 6];
        int di7 = d[pos_y + i][pos_x + 7];

        int ei0 = di0 + di4;
        int ei1 = -di3 + di5 - di7 - (di7 >> 1);
        int ei2 = di0 - di4;
        int ei3 = di1 + di7 - di3 - (di3 >> 1);
        int ei4 = (di2 >> 1) - di6;
        int ei5 = -di1 + di7 + di5 + (di5 >> 1);
        int ei6 = di2 + (di6 >> 1);
        int ei7 = di3 + di5 + di1 + (di1 >> 1);

        int fi0 = ei0 + ei6;
        int fi1 = ei1 + (ei7 >> 2);
        int fi2 = ei2 + ei4;
        int fi3 = ei3 + (ei5 >> 2);
        int fi4 = ei2 - ei4;
        int fi5 = (ei3 >> 2) - ei5;
        int fi6 = ei0 - ei6;
        int fi7 = ei7 - (ei1 >> 2);

        g[i][0] = fi0 + fi7;
        g[i][1] = fi2 + fi5;
        g[i][2] = fi4 + fi3;
        g[i][3] = fi6 + fi1;
        g[i][4] = fi6 - fi1;
        g[i][5] = fi4 - fi3;
        g[i][6] = fi2 - fi5;
        g[i][7] = fi0 - fi7;
    }

    // vertical column of the resulting matrix
    for (int j = 0; j < 8; ++j) {
        int g0j = g[0][j];
        int g1j = g[1][j];
        int g2j = g[2][j];
        int g3j = g[3][j];
        int g4j = g[4][j];
        int g5j = g[5][j];
        int g6j = g[6][j];
        int g7j = g[7][j];

        int h0j = g0j + g4j;
        int h1j = -g3j + g5j - g7j - (g7j >> 1);
        int h2j = g0j - g4j;
        int h3j = g1j + g7j - g3j - (g3j >> 1);
        int h4j = (g2j >> 1) - g6j;
        int h5j = -g1j + g7j + g5j + (g5j >> 1);
        int h6j = g2j + (g6j >> 1);
        int h7j = g3j + g5j + g1j + (g1j >> 1);

        int k0j = h0j + h6j;
        int k1j = h1j + (h7j >> 2);
        int k2j = h2j + h4j;
        int k3j = h3j + (h5j >> 2);
        int k4j = h2j - h4j;
        int k5j = (h3j >> 2) - h5j;
        int k6j = h0j - h6j;
        int k7j = h7j - (h1j >> 2);

        int m0j = k0j + k7j;
        int m1j = k2j + k5j;
        int m2j = k4j + k3j;
        int m3j = k6j + k1j;
        int m4j = k6j - k1j;
        int m5j = k4j - k3j;
        int m6j = k2j - k5j;
        int m7j = k0j - k7j;

        r[pos_y + 0][pos_x + j] = (m0j + (1 << 5)) >> 6;
        r[pos_y + 1][pos_x + j] = (m1j + (1 << 5)) >> 6;
        r[pos_y + 2][pos_x + j] = (m2j + (1 << 5)) >> 6;
        r[pos_y + 3][pos_x + j] = (m3j + (1 << 5)) >> 6;
        r[pos_y + 4][pos_x + j] = (m4j + (1 << 5)) >> 6;
        r[pos_y + 5][pos_x + j] = (m5j + (1 << 5)) >> 6;
        r[pos_y + 6][pos_x + j] = (m6j + (1 << 5)) >> 6;
        r[pos_y + 7][pos_x + j] = (m7j + (1 << 5)) >> 6;
    }
}


void Transform::bypass_4x4(int r[16][16], int f[16][16], int ioff, int joff, uint8_t pred_mode)
{
    if (pred_mode == IntraPrediction::Intra_4x4_Vertical) {
        for (int i = 0; i < 4; ++i) {
            f[joff + 0][ioff + i] = r[joff + 0][ioff + i];
            for (int j = 1; j < 4; ++j)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j - 1][ioff + i];
        }
    } else if (pred_mode == IntraPrediction::Intra_4x4_Horizontal) {
        for (int j = 0; j < 4; ++j) {
            f[joff + j][ioff + 0] = r[joff + j][ioff + 0];
            for (int i = 1; i < 4; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j][ioff + i - 1];
        }
    } else {
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i];
        }
    }
}

void Transform::bypass_8x8(int r[16][16], int f[16][16], int ioff, int joff, uint8_t pred_mode)
{
    if (pred_mode == IntraPrediction::Intra_8x8_Vertical) {
        for (int i = 0; i < 8; ++i) {
            f[joff + 0][ioff + i] = r[joff + 0][ioff + i];
            for (int j = 1; j < 8; ++j)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j - 1][ioff + i];
        }
    } else if (pred_mode == IntraPrediction::Intra_8x8_Horizontal) {
        for (int j = 0; j < 8; ++j) {
            f[joff + j][ioff + 0] = r[joff + j][ioff + 0];
            for (int i = 1; i < 8; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i] + f[joff + j][ioff + i - 1];
        }
    } else {
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 8; ++i)
                f[joff + j][ioff + i] = r[joff + j][ioff + i];
        }
    }
}

void Transform::bypass_16x16(int r[16][16], int f[16][16], int ioff, int joff, uint8_t pred_mode)
{
    if (pred_mode == IntraPrediction::Intra_16x16_Vertical) {
        for (int i = 0; i < 16; ++i) {
            f[0][i] = r[0][i];
            for (int j = 1; j < 16; ++j)
                f[j][i] = r[j][i] + f[j - 1][i];
        }
    } else if (pred_mode == IntraPrediction::Intra_16x16_Horizontal) {
        for (int j = 0; j < 16; ++j) {
            f[j][0] = r[j][0];
            for (int i = 1; i < 16; ++i)
                f[j][i] = r[j][i] + f[j][i - 1];
        }
    } else {
        for (int j = 0; j < 16; ++j) {
            for (int i = 0; i < 16; ++i)
                f[j][i] = r[j][i];
        }
    }
}

void Transform::bypass_chroma(int r[16][16], int f[16][16], int nW, int nH, uint8_t pred_mode)
{
    if (pred_mode == IntraPrediction::Intra_Chroma_Vertical) {
        for (int i = 0; i < nW; ++i) {
            f[0][i] = r[0][i];
            for (int j = 1; j < nH; ++j)
                f[j][i] = r[j][i] + f[j - 1][i];
        }
    } else if (pred_mode == IntraPrediction::Intra_Chroma_Horizontal) {
        for (int j = 0; j < nH; ++j) {
            f[j][0] = r[j][0];
            for (int i = 1; i < nW; ++i)
                f[j][i] = r[j][i] + f[j][i - 1];
        }
    } else {
        for (int j = 0; j < nH; ++j) {
            for (int i = 0; i < nW; ++i)
                f[j][i] = r[j][i];
        }
    }
}


void Transform::transform_luma_dc(mb_t* mb, ColorPlane pl)
{
    if (!mb->TransformBypassModeFlag) {
        slice_t* slice = mb->p_Slice;
        sps_t* sps = slice->active_sps;

        int transform_pl = sps->separate_colour_plane_flag ? PLANE_Y : pl;
        int (*cof)[16] = this->cof[transform_pl];
        int (*dcY)[16] = this->cof[transform_pl];
        int c[4][4];
        int f[4][4];

        int qP = mb->qp_scaled[transform_pl];
        int scale = this->InvLevelScale4x4_Intra[pl][qP % 6][0][0];

        for (int i = 0; i < 16; i += 4) {
            for (int j = 0; j < 16; j += 4)
                c[i / 4][j / 4] = cof[i][j];
        }

        this->ihadamard_4x4(c, f);

        for (int i = 0; i < 16; i += 4) {
            for (int j = 0; j < 16; j += 4) {
                if (qP >= 36)
                    dcY[i][j] = (f[i / 4][j / 4] * scale) << (qP / 6 - 6);
                else
                    dcY[i][j] = (f[i / 4][j / 4] * scale + (1 << (5 - qP / 6))) >> (6 - qP / 6);
            }
        }
    }
}

void Transform::transform_chroma_dc(mb_t* mb, ColorPlane pl)
{
    if (!mb->TransformBypassModeFlag) {
        slice_t* slice = mb->p_Slice;
        sps_t* sps = slice->active_sps;

        bool smb = (slice->slice_type == SP_slice && !mb->is_intra_block) ||
                   (slice->slice_type == SI_slice && mb->mb_type == SI4MB);
        int (*cof)[16] = this->cof[pl];
        int (*dcC)[16] = this->cof[pl];

        int qP = mb->qp_scaled[pl];
        int scale = mb->is_intra_block ?
            this->InvLevelScale4x4_Intra[pl][qP % 6][0][0] :
            this->InvLevelScale4x4_Inter[pl][qP % 6][0][0];

        if (sps->ChromaArrayType == 1 && !smb) {
            int c[2][2];
            int f[2][2];
            for (int i = 0; i < sps->MbHeightC; i += 4) {
                for (int j = 0; j < sps->MbWidthC; j += 4)
                    c[i / 4][j / 4] = cof[i][j];
            }

            this->ihadamard_2x2(c, f);

            for (int i = 0; i < sps->MbHeightC; i += 4) {
                for (int j = 0; j < sps->MbWidthC; j += 4)
                    dcC[i][j] = ((f[i / 4][j / 4] * scale) << (qP / 6)) >> 5;
            }
        }
        if (sps->ChromaArrayType == 2) {
            int c[4][2];
            int f[4][2];
            for (int i = 0; i < sps->MbHeightC; i += 4) {
                for (int j = 0; j < sps->MbWidthC; j += 4)
                    c[i / 4][j / 4] = cof[i][j];
            }

            this->ihadamard_2x4(c, f);

            for (int i = 0; i < sps->MbHeightC; i += 4) {
                for (int j = 0; j < sps->MbWidthC; j += 4) {
                    if (qP >= 36)
                        dcC[i][j] = (f[i / 4][j / 4] * scale) << (qP / 6 - 6);
                    else
                        dcC[i][j] = (f[i / 4][j / 4] * scale + (1 << (5 - qP / 6))) >> (6 - qP / 6);
                }
            }
        }
    }
}


void Transform::construction(mb_t* mb, ColorPlane pl, int ioff, int joff, int nW, int nH)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;
    int block8x8 = (joff / 8) * 2 + (ioff / 8);

    int (*mb_rres)[16] = this->mb_rres[pl];
    imgpel** mb_pred = slice->mb_pred[pl];
    imgpel (*mb_rec)[16] = this->mb_rec[pl];

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        for (int j = 0; j < nH; ++j) {
            for (int i = 0; i < nW; ++i)
                mb_rec[joff + j][ioff + i] = (imgpel) clip1(max_pel_value_comp, mb_rres[joff + j][ioff + i] + mb_pred[joff + j][ioff + i]);
        }
    } else {
        for (int j = 0; j < nH; ++j)
            memcpy(&mb_rec[joff + j][ioff], &mb_pred[joff + j][ioff], nW * sizeof(imgpel));
    }

    for (int j = 0; j < nH; ++j)
        memcpy(&curr_img[mb->mb.y * 16 + joff + j][mb->mb.x * 16 + ioff], &mb_rec[joff + j][ioff], nW * sizeof (imgpel));
}

void Transform::construction_16x16(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int (*mb_rres)[16] = this->mb_rres[pl];
    imgpel **mb_pred = slice->mb_pred[pl];
    imgpel (*mb_rec)[16] = this->mb_rec [pl];

    for (int j = 0; j < 16; j++) {
        for (int i = 0; i < 16; i++)
            mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_rres[j][i] + mb_pred[j][i]);
    }

    for (int j = 0; j < 16; ++j)
        memcpy(&curr_img[mb->mb.y * 16 + joff + j][mb->mb.x * 16 + ioff], &mb_rec[joff + j][ioff], 16 * sizeof (imgpel));
}

void Transform::construction_chroma(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    int (*mb_rres)[16] = this->mb_rres[pl];
    imgpel **mb_pred = slice->mb_pred[pl];
    imgpel (*mb_rec)[16] = this->mb_rec [pl];

    for (int j = 0; j < sps->MbHeightC; j++) {
        for (int i = 0; i < sps->MbWidthC; i++)
            mb_rec[j][i] = (imgpel) clip1(max_pel_value_comp, mb_rres[j][i] + mb_pred[j][i]);
    }

    for (int joff = 0; joff < sps->MbHeightC; joff += 4) {
        for (int ioff = 0; ioff < sps->MbWidthC; ioff += 4)
            for (int j = 0; j < 4; ++j)
                memcpy(&curr_img[mb->mb.y * sps->MbHeightC + joff + j][mb->mb.x * sps->MbWidthC + ioff],
                       &mb_rec[joff + j][ioff], 4 * sizeof (imgpel));
    }
}

void Transform::inverse_transform_4x4(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    int block8x8 = (joff / 8) * 2 + (ioff / 8);

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        int i4x4 = ((joff / 4) / 2) * 8 + ((joff / 4) % 2) * 2 +
                   ((ioff / 4) / 2) * 4 + ((ioff / 4) % 2);
        uint8_t pred_mode = mb->Intra4x4PredMode[i4x4];
        if (mb->TransformBypassModeFlag)
            this->bypass_4x4(this->cof[pl], this->mb_rres[pl], ioff, joff, pred_mode);
        else
            this->inverse_4x4(this->cof[pl], this->mb_rres[pl], joff, ioff);
    }

    this->construction(mb, pl, ioff, joff, 4, 4);
}

void Transform::inverse_transform_8x8(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    int block8x8 = (joff / 8) * 2 + (ioff / 8);

    if (mb->CodedBlockPatternLuma & (1 << block8x8)) {
        uint8_t pred_mode = mb->Intra8x8PredMode[block8x8];
        if (mb->TransformBypassModeFlag)
            this->bypass_8x8(this->cof[pl], this->mb_rres[pl], ioff, joff, pred_mode);
        else
            this->inverse_8x8(this->cof[pl], this->mb_rres[pl], joff, ioff);
    }

    this->construction(mb, pl, ioff, joff, 8, 8);
}

void Transform::inverse_transform_16x16(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    uint8_t pred_mode = mb->Intra16x16PredMode;
    if (mb->TransformBypassModeFlag)
        this->bypass_16x16(this->cof[pl], this->mb_rres[pl], ioff, joff, pred_mode);
    else {
        for (int j = 0; j < 16; j += 4) {
            for (int i = 0; i < 16; i += 4)
                this->inverse_4x4(this->cof[pl], this->mb_rres[pl], j, i);
        }
    }

    this->construction_16x16(mb, pl, ioff, joff);
}

void Transform::inverse_transform_chroma(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    uint8_t pred_mode = mb->intra_chroma_pred_mode;
    if (mb->TransformBypassModeFlag)
        this->bypass_chroma(this->cof[pl], this->mb_rres[pl], sps->MbWidthC, sps->MbHeightC, pred_mode);
    else {
        for (int joff = 0; joff < sps->MbHeightC; joff += 4) {
            for (int ioff = 0; ioff < sps->MbWidthC; ioff += 4)
                this->inverse_4x4(this->cof[pl], this->mb_rres[pl], joff, ioff);
        }
    }

    this->construction_chroma(mb, pl, 0, 0);
}

void Transform::inverse_transform_inter(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;

    if (mb->CodedBlockPatternLuma) {
        if (!mb->transform_size_8x8_flag) {
            for (int y = 0; y < 16; y += 4) {
                for (int x = 0; x < 16; x += 4)
                    this->inverse_transform_4x4(mb, pl, x, y);
            }
        } else {
            for (int y = 0; y < 16; y += 8) {
                for (int x = 0; x < 16; x += 8)
                    this->inverse_transform_8x8(mb, pl, x, y);
            }
        }
    } else {
        for (int j = 0; j < 16; ++j)
            memcpy(&curr_img[mb->mb.y * 16 + j][mb->mb.x * 16], &slice->mb_pred[pl][j][0], 16 * sizeof (imgpel));
    }

    if (mb->CodedBlockPatternLuma)
        slice->parser.is_reset_coeff = false;

    if (sps->chroma_format_idc == YUV400 || sps->chroma_format_idc == YUV444)
        return;

    for (int uv = 0; uv < 2; ++uv) {
        imgpel **curUV = &dec_picture->imgUV[uv][mb->mb.y * sps->MbHeightC]; 
        imgpel **mb_pred = slice->mb_pred[uv + 1];

        if (mb->CodedBlockPatternChroma)
            this->inverse_transform_chroma(mb, (ColorPlane)(uv + 1));
        else {
            for (int j = 0; j < sps->MbHeightC; ++j)
                memcpy(&curUV[j][mb->mb.x * sps->MbWidthC], &mb_pred[j][0], sps->MbWidthC * sizeof (imgpel));
        }
    }

    if (mb->CodedBlockPatternChroma)
        slice->parser.is_reset_coeff_cr = false;
}


static const uint8_t A[4][4] = {
    { 16, 20, 16, 20 },
    { 20, 25, 20, 25 },
    { 16, 20, 16, 20 },
    { 20, 25, 20, 25 }
};

static const uint16_t LevelScale2[6][4][4] = {
    {{ 13107,  8066, 13107,  8066 },
     {  8066,  5243,  8066,  5243 },
     { 13107,  8066, 13107,  8066 },
     {  8066,  5243,  8066,  5243 }},
    {{ 11916,  7490, 11916,  7490 },
     {  7490,  4660,  7490,  4660 },
     { 11916,  7490, 11916,  7490 },
     {  7490,  4660,  7490,  4660 }},
    {{ 10082,  6554, 10082,  6554 },
     {  6554,  4194,  6554,  4194 },
     { 10082,  6554, 10082,  6554 },
     {  6554,  4194,  6554,  4194 }},
    {{  9362,  5825,  9362,  5825 },
     {  5825,  3647,  5825,  3647 },
     {  9362,  5825,  9362,  5825 },
     {  5825,  3647,  5825,  3647 }},
    {{  8192,  5243,  8192,  5243 },
     {  5243,  3355,  5243,  3355 },
     {  8192,  5243,  8192,  5243 },
     {  5243,  3355,  5243,  3355 }},
    {{  7282,  4559,  7282,  4559 },
     {  4559,  2893,  4559,  2893 },
     {  7282,  4559,  7282,  4559 },
     {  4559,  2893,  4559,  2893 }}
};

void Transform::itrans_sp(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int QpY = (slice->slice_type == SI_slice) ? slice->QsY : slice->SliceQpY;
    int QsY = slice->QsY;

    int    (*cof    )[16] = this->cof    [pl];
    int    (*mb_rres)[16] = this->mb_rres[pl];
    imgpel (*mb_rec )[16] = this->mb_rec [pl];
    int max_pel_value_comp = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;

    const int (*InvLevelScale4x4)  [4] = dequant_coef[QpY % 6];
    const int (*InvLevelScale4x4SP)[4] = dequant_coef[QsY % 6];  

    int p[16][16];
    int c[16][16];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j)
            p[joff + i][ioff + j] = slice->mb_pred[pl][joff + i][ioff + j];
    }

    this->forward_4x4(p, c, joff, ioff);

    int crij, cpij, csij, cij;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            //crij = (cof[joff + j][ioff + i] >> (QpY / 6)) / InvLevelScale4x4[j][i];
            crij = cof[joff + j][ioff + i];
            cpij = c  [joff + j][ioff + i];
            if (slice->sp_for_switch_flag || slice->slice_type == SI_slice) {
                csij = sign(cpij) * ((abs(cpij) * LevelScale2[QsY % 6][j][i] + (1 << (14 + QsY / 6))) >> (15 + QsY / 6));
                cij  = crij + csij;
            } else {
                //c[j][i] += ((icof * InvLevelScale4x4[j][i] * A[j][i] << (QpY / 6)) >> 6);
                csij = cpij + (((crij * InvLevelScale4x4[j][i] * A[j][i]) << (QpY / 6)) >> 10);
                cij  = sign(csij) * ((abs(csij) * LevelScale2[QsY % 6][j][i] + (1 << (14 + QsY / 6))) >> (15 + QsY / 6));
            }

            //cof[joff + j][ioff + i] = cij * InvLevelScale4x4SP[j][i] << (QsY / 6);
            if (QsY >= 24)
                cof[joff + j][ioff + i] = (cij * InvLevelScale4x4SP[j][i]) << (QsY / 6 - 4);
            else
                cof[joff + j][ioff + i] = (cij * InvLevelScale4x4SP[j][i] + (1 << (3 - QsY / 6))) >> (4 - QsY / 6);
        }
    }

    this->inverse_4x4(cof, mb_rres, joff, ioff);

    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i)
            mb_rec[joff + j][ioff + i] = (imgpel)clip1(max_pel_value_comp, mb_rres[joff + j][ioff + i]);
    }
}


void Transform::itrans_sp_cr(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    imgpel** mb_pred = slice->mb_pred[pl];
    int (*cof)[16] = this->cof[pl];

    int QpC = mb->QpC[pl - 1];
    int QsC = mb->QsC[pl - 1];

    const int (*InvLevelScale4x4)  [4] = dequant_coef[QpC % 6];
    const int (*InvLevelScale4x4SP)[4] = dequant_coef[QsC % 6];

    int p[16][16];
    int c[16][16];
    for (int j = 0; j < sps->MbHeightC; ++j) {
        for (int i = 0; i < sps->MbWidthC; ++i)
            p[j][i] = mb_pred[j][i];
            //mb_pred[j][i] = 0;
    }

    for (int j = 0; j < sps->MbHeightC; j += 4) {
        for (int i = 0; i < sps->MbWidthC; i += 4)
            this->forward_4x4(p, c, j, i);
    }

    int mp1[2][2];
    //     2X2 transform of DC coeffs.
    mp1[0][0] = (c[0][0] + c[4][0] + c[0][4] + c[4][4]);
    mp1[0][1] = (c[0][0] - c[4][0] + c[0][4] - c[4][4]);
    mp1[1][0] = (c[0][0] + c[4][0] - c[0][4] - c[4][4]);
    mp1[1][1] = (c[0][0] - c[4][0] - c[0][4] + c[4][4]);

    int crij, cpij, csij, cij;
    for (int n2 = 0; n2 < 2; ++n2) {
        for (int n1 = 0; n1 < 2; ++n1) {
            crij = cof[n2 * 4][n1 * 4];
            cpij = mp1[n2][n1];
            if (slice->sp_for_switch_flag || slice->slice_type == SI_slice) {
                csij = (sign(cpij) * (abs(cpij) * LevelScale2[QsC][0][0] + (1 << (15 + QsC / 6)))) >> (16 + QsC / 6);
                cij  = csij + cpij;
            } else {
                //csij = cpij + (((crij * InvLevelScale4x4[0][0] * A[0][0]) << qp_per) >> 5);
                csij = cpij + (((crij * InvLevelScale4x4[0][0] * A[0][0]) << (QpC / 6)) >> 9);
                cij  = (sign(csij) * (abs(csij) * LevelScale2[QsC][0][0] + (1 << (15 + QsC / 6)))) >> (16 + QsC / 6);
            }
            mp1[n2][n1] = cij * InvLevelScale4x4SP[0][0] << (QpC / 6);
        }
    }

    for (int n2 = 0; n2 < sps->MbHeightC; n2 += 4) {
        for (int n1 = 0; n1 < sps->MbWidthC; n1 += 4) {
            for (int j = 0; j < 4; ++j) {
                for (int i = 0; i < 4; ++i) {
                    crij = cof[n2 + j][n1 + i];
                    cpij = p[n2 + j][n1 + i];
                    //icof = (cof[n2 + j][n1 + i] >> qp_per) / dequant_coef[qp_rem][j][i];
                    if (slice->sp_for_switch_flag || slice->slice_type == SI_slice) {
                        csij = (sign(cpij) * (abs(cpij) * LevelScale2[QsC % 6][j][i] + (1 << (14 + QsC / 6)))) >> (15 + QsC / 6);
                        cij  = crij + csij;
                    } else {
                        csij = cpij + ((crij * InvLevelScale4x4[j][i] * A[j][i] << (QpC / 6)) >> 9);
                        cij  = (sign(csij) * (abs(csij) * LevelScale2[QsC % 6][j][i] + (1 << (14 + QsC / 6)))) >> (15 + QsC / 6);
                    }
                    cof[n2 + j][n1 + i] = cij * InvLevelScale4x4SP[j][i] << (QpC / 6);
                }
            }
        }
    }

    cof[0][0] = (mp1[0][0] + mp1[0][1] + mp1[1][0] + mp1[1][1]) >> 1;
    cof[0][4] = (mp1[0][0] + mp1[0][1] - mp1[1][0] - mp1[1][1]) >> 1;
    cof[4][0] = (mp1[0][0] - mp1[0][1] + mp1[1][0] - mp1[1][1]) >> 1;
    cof[4][4] = (mp1[0][0] - mp1[0][1] - mp1[1][0] + mp1[1][1]) >> 1;
}

void Transform::inverse_transform_sp(mb_t* mb, ColorPlane pl)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;
    imgpel** curr_img = pl ? dec_picture->imgUV[pl - 1] : dec_picture->imgY;

    if (!mb->transform_size_8x8_flag) {
        for (int y = 0; y < 16; y += 4) {
            for (int x = 0; x < 16; x += 4)
                this->itrans_sp(mb, pl, x, y);
        }
    } else {
        for (int y = 0; y < 16; y += 8) {
            for (int x = 0; x < 16; x += 8)
                this->inverse_transform_8x8(mb, pl, x, y);
        }
    }

    for (int j = 0; j < 16; ++j)
        memcpy(&curr_img[mb->mb.y * 16 + j][mb->mb.x * 16], &this->mb_rec[pl][j][0], 16 * sizeof (imgpel));

    slice->parser.is_reset_coeff = false;

    if (sps->chroma_format_idc == YUV400 || sps->chroma_format_idc == YUV444)
        return;

    for (int uv = 0; uv < 2; ++uv) {
        this->itrans_sp_cr(mb, (ColorPlane)(uv + 1));
        this->inverse_transform_chroma(mb, (ColorPlane)(uv + 1));
    }

    slice->parser.is_reset_coeff_cr = false;
}


}
}
