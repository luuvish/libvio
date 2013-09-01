/*!
 ***********************************************************************
 * \file read_comp_cavlc.c
 *
 * \brief
 *     Read Coefficient Components (CAVLC version)
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

#include "global.h"
#include "slice.h"
#include "data_partition.h"
#include "macroblock.h"
#include "mb_read.h"
#include "transform.h"
#include "neighbour.h"

#include "mb_read_syntax.h"


#define IS_I16MB(MB) ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)

static inline int rshift_rnd_sf(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}


// Table 9-5 coeff_token mapping to TotalCoeff(coeff_token) and TrailingOnes(coeff_token)
static const uint8_t coeff_token_length[5][4][17] = {
    //  0 <= nC < 2
    { {  1,  6,  8,  9, 10, 11, 13, 13, 13, 14, 14, 15, 15, 16, 16, 16, 16 },
      {  0,  2,  6,  8,  9, 10, 11, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16 },
      {  0,  0,  3,  7,  8,  9, 10, 11, 13, 13, 14, 14, 15, 15, 16, 16, 16 },
      {  0,  0,  0,  5,  6,  7,  8,  9, 10, 11, 13, 14, 14, 15, 15, 16, 16 } },
    // 2 <= nC < 4
    { {  2,  6,  6,  7,  8,  8,  9, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14 },
      {  0,  2,  5,  6,  6,  7,  8,  9, 11, 11, 12, 12, 13, 13, 14, 14, 14 },
      {  0,  0,  3,  6,  6,  7,  8,  9, 11, 11, 12, 12, 13, 13, 13, 14, 14 },
      {  0,  0,  0,  4,  4,  5,  6,  6,  7,  9, 11, 11, 12, 13, 13, 13, 14 } },
    // 4 <= nC < 8
    { {  4,  6,  6,  6,  7,  7,  7,  7,  8,  8,  9,  9,  9, 10, 10, 10, 10 },
      {  0,  4,  5,  5,  5,  5,  6,  6,  7,  8,  8,  9,  9,  9, 10, 10, 10 },
      {  0,  0,  4,  5,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9, 10, 10, 10 },
      {  0,  0,  0,  4,  4,  4,  4,  4,  5,  6,  7,  8,  8,  9, 10, 10, 10 } },
    // nC == -1 
    { {  2,  6,  6,  6,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  1,  6,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  3,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 } },
    // nC == -2 
    { {  1,  7,  7,  9,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  2,  7,  7,  9, 10, 11, 12, 12,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  3,  7,  7,  9, 10, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  5,  6,  7,  7, 10, 11,  0,  0,  0,  0,  0,  0,  0,  0 } }
};

static const uint8_t coeff_token_code[5][4][17] = {
    // 0 <= nC < 2
    { {  1,  5,  7,  7,  7,  7, 15, 11,  8, 15, 11, 15, 11, 15, 11,  7,  4 },
      {  0,  1,  4,  6,  6,  6,  6, 14, 10, 14, 10, 14, 10,  1, 14, 10,  6 },
      {  0,  0,  1,  5,  5,  5,  5,  5, 13,  9, 13,  9, 13,  9, 13,  9,  5 },
      {  0,  0,  0,  3,  3,  4,  4,  4,  4,  4, 12, 12,  8, 12,  8, 12,  8 } },
    // 2 <= nC < 4
    { {  3, 11,  7,  7,  7,  4,  7, 15, 11, 15, 11,  8, 15, 11,  7,  9,  7 },
      {  0,  2,  7, 10,  6,  6,  6,  6, 14, 10, 14, 10, 14, 10, 11,  8,  6 },
      {  0,  0,  3,  9,  5,  5,  5,  5, 13,  9, 13,  9, 13,  9,  6, 10,  5 },
      {  0,  0,  0,  5,  4,  6,  8,  4,  4,  4, 12,  8, 12, 12,  8,  1,  4 } },
    // 4 <= nC < 8
    { { 15, 15, 11,  8, 15, 11,  9,  8, 15, 11, 15, 11,  8, 13,  9,  5,  1 },
      {  0, 14, 15, 12, 10,  8, 14, 10, 14, 14, 10, 14, 10,  7, 12,  8,  4 },
      {  0,  0, 13, 14, 11,  9, 13,  9, 13, 10, 13,  9, 13,  9, 11,  7,  3 },
      {  0,  0,  0, 12, 11, 10,  9,  8, 13, 12, 12, 12,  8, 12, 10,  6,  2 } },
    // nC == -1
    { {  1,  7,  4,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  1,  6,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 } },
    // nC == -2
    { {  1, 15, 14,  7,  6,  7,  7,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  1, 13, 12,  5,  6,  6,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  1, 11, 10,  4,  5,  5,  4,  0,  0,  0,  0,  0,  0,  0,  0 },
      {  0,  0,  0,  1,  1,  9,  8,  4,  4,  0,  0,  0,  0,  0,  0,  0,  0 } }
};

// Table 9-7 total_zeros tables for 4x4 blocks with tzVlcIndex 1 to 7
// Table 9-8 total_zeros tables for 4x4 blocks with tzVlcIndex 8 to 15
// Table 9-9 total_zeros tables for chroma DC 2x2 and 2x4 blocks
static const uint8_t total_zeros_length[3][15][16] = {
    // YUV420
    { { 1, 2, 3, 3 },
      { 1, 2, 2 },
      { 1, 1 } },
    // YUV422
    { { 1, 3, 3, 4, 4, 4, 5, 5 },
      { 3, 2, 3, 3, 3, 3, 3 },
      { 3, 3, 2, 2, 3, 3 },
      { 3, 2, 2, 2, 3 },
      { 2, 2, 2, 2 },
      { 2, 2, 1 },
      { 1, 1 } },
    // YUV444
    { { 1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9 },
      { 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6 },
      { 4, 3, 3, 3, 4, 4, 3, 3, 4, 5, 5, 6, 5, 6 },
      { 5, 3, 4, 4, 3, 3, 3, 4, 3, 4, 5, 5, 5 },
      { 4, 4, 4, 3, 3, 3, 3, 3, 4, 5, 4, 5 },
      { 6, 5, 3, 3, 3, 3, 3, 3, 4, 3, 6 },
      { 6, 5, 3, 3, 3, 2, 3, 4, 3, 6 },
      { 6, 4, 5, 3, 2, 2, 3, 3, 6 },
      { 6, 6, 4, 2, 2, 3, 2, 5 },
      { 5, 5, 3, 2, 2, 2, 4 },
      { 4, 4, 3, 3, 1, 3 },
      { 4, 4, 2, 1, 3 },
      { 3, 3, 1, 2 },
      { 2, 2, 1 },
      { 1, 1 } }
};

static const uint8_t total_zeros_code[3][15][16] = {
    // YUV420
    { { 1, 1 , 1 , 0 },
      { 1, 1 , 0 },
      { 1, 0 } },
    // YUV422
    { { 1, 2, 3, 2, 3, 1, 1, 0 },
      { 0, 1, 1, 4, 5, 6, 7 },
      { 0, 1, 1, 2, 6, 7 },
      { 6, 0, 1, 2, 7 },
      { 0, 1, 2, 3 },
      { 0, 1, 1 },
      { 0, 1 } },
    // YUV444
    { { 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1 },
      { 7, 6, 5, 4, 3, 5, 4, 3, 2, 3, 2, 3, 2, 1, 0 },
      { 5, 7, 6, 5, 4, 3, 4, 3, 2, 3, 2, 1, 1, 0 },
      { 3, 7, 5, 4, 6, 5, 4, 3, 3, 2, 2, 1, 0 },
      { 5, 4, 3, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
      { 1, 1, 7, 6, 5, 4, 3, 2, 1, 1, 0 },
      { 1, 1, 5, 4, 3, 3, 2, 1, 1, 0 },
      { 1, 1, 1, 3, 3, 2, 2, 1, 0 },
      { 1, 0, 1, 3, 2, 1, 1, 1 },
      { 1, 0, 1, 3, 2, 1, 1 },
      { 0, 1, 1, 2, 1, 3 },
      { 0, 1, 1, 1, 1 },
      { 0, 1, 1, 1 },
      { 0, 1, 1 },
      { 0, 1 } }
};

// Table 9-10 Tables for run_before
static const uint8_t run_before_length[15][16] = {
    { 1, 1 },
    { 1, 2, 2 },
    { 2, 2, 2, 2 },
    { 2, 2, 2, 3, 3 },
    { 2, 2, 3, 3, 3, 3 },
    { 2, 3, 3, 3, 3, 3, 3 },
    { 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
};

static const uint8_t run_before_code[15][16] = {
    { 1, 0 },
    { 1, 1, 0 },
    { 3, 2, 1, 0 },
    { 3, 2, 1, 1, 0 },
    { 3, 2, 3, 2, 1, 0 },
    { 3, 0, 1, 3, 2, 5, 4 },
    { 7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
};


uint8_t macroblock_t::parse_coeff_token(int nC)
{
    slice_t* slice = this->p_Slice;
    data_partition_t* dp = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];

    if (nC >= 8) {
        int code = dp->read_bits(6);
        int TotalCoeff   = (code >> 2);
        int TrailingOnes = (code & 3);
        if (TotalCoeff == 0 && TrailingOnes == 3)
            TrailingOnes = 0;
        else
            TotalCoeff++;
        return (TotalCoeff << 2) | (TrailingOnes);
    }

    int tab = (nC == -2) ? 4 : (nC == -1) ? 3 : (nC < 2) ? 0 : (nC < 4) ? 1 : (nC < 8) ? 2 : 5;

    for (int TrailingOnes = 0; TrailingOnes < 4; TrailingOnes++) {
        for (int TotalCoeff = 0; TotalCoeff < 17; TotalCoeff++) {
            int length = coeff_token_length[tab][TrailingOnes][TotalCoeff];
            int code   = coeff_token_code  [tab][TrailingOnes][TotalCoeff];
            if (length > 0 && dp->next_bits(length) == code) {
                dp->read_bits(length);
                return (TotalCoeff << 2) | (TrailingOnes);
            }
        }
    }

    assert(false);
    return -1;
}

/*
static int16_t parse_level(data_partition_t *currStream, uint8_t level_prefix, uint8_t suffixLength)
{
    int level, sign;

    if (suffixLength == 0) {
        if (level_prefix < 14) {
            sign  = level_prefix & 1;
            level = (level_prefix >> 1) + 1;
        } else if (level_prefix == 14) {
            // escape
            int level_suffix = currStream->u(4);
            sign  = (level_suffix & 1);
            level = (level_suffix >> 1) + 8;
        } else {
            // escape
            int level_suffix = currStream->u(level_prefix - 3);
            sign  = (level_suffix & 1);
            level = (level_suffix >> 1) + (1 << (level_prefix - 4)) - 2047 + 15;
        }
    } else {
        if (level_prefix < 15) {
            int level_suffix = currStream->u(suffixLength);
            sign  = (level_suffix & 1);
            level = (level_suffix >> 1) + (level_prefix << (suffixLength - 1)) + 1;
        } else { // escape
            int level_suffix = currStream->u(level_prefix - 3);
            sign  = (level_suffix & 1);
            level = (level_suffix >> 1) + (1 << (level_prefix - 4)) - 2047 + (15 << (suffixLength - 1));
        }
    }

    return sign ? -level : level;
}
*/

uint8_t macroblock_t::parse_total_zeros(int yuv, int tzVlcIndex)
{
    slice_t* slice = this->p_Slice;
    data_partition_t* dp = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];

    int tab = tzVlcIndex - 1;

    for (int total_zeros = 0; total_zeros < 16; total_zeros++) {
        int length = total_zeros_length[yuv][tab][total_zeros];
        int code   = total_zeros_code  [yuv][tab][total_zeros];
        if (length > 0 && dp->next_bits(length) == code) {
            dp->read_bits(length);
            return total_zeros;
        }
    }

    assert(false);
    return -1;
}

uint8_t macroblock_t::parse_run_before(uint8_t zerosLeft)
{
    slice_t* slice = this->p_Slice;
    data_partition_t* dp = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];

    int tab = imin(zerosLeft, 7) - 1;

    for (int run_before = 0; run_before < 16; run_before++) {
        int length = run_before_length[tab][run_before];
        int code   = run_before_code  [tab][run_before];
        if (length > 0 && dp->next_bits(length) == code) {
            dp->read_bits(length);
            return run_before;
        }
    }

    assert(false);
    return -1;
}


void macroblock_t::read_coeff_4x4_CAVLC(int maxNumCoeff, int nC,
                                        int levelVal[16], int runVal[16], int *number_coefficients)
{
    slice_t* slice = this->p_Slice;
    data_partition_t* dp = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];

    memset(levelVal, 0, maxNumCoeff * sizeof(int));
    memset(runVal,   0, maxNumCoeff * sizeof(int));

    uint8_t coeff_token  = this->parse_coeff_token(nC);
    uint8_t TotalCoeff   = coeff_token >> 2;
    uint8_t TrailingOnes = coeff_token % 4;

    if (TotalCoeff > 0) {
        int suffixLength = (TotalCoeff > 10 && TrailingOnes < 3 ? 1 : 0);

        if (TrailingOnes) {
            int code = dp->f(TrailingOnes);
            int ntr = TrailingOnes;
            for (int i = TotalCoeff - 1; i > TotalCoeff - 1 - TrailingOnes; i--) {
                int trailing_ones_sign_flag = (code >> (--ntr)) & 1;
                levelVal[i] = 1 - 2 * trailing_ones_sign_flag;
            }
        }

        for (int i = TotalCoeff - 1 - TrailingOnes; i >= 0; i--) {
            int level_prefix, level_suffix;
            int levelSuffixSize, levelCode;

            int leadingZeroBits = -1;
            for (int b = 0; !b; leadingZeroBits++)
                b = dp->read_bits(1);
            level_prefix = leadingZeroBits;

            levelSuffixSize = (level_prefix == 14 && suffixLength == 0) ? 4 :
                              (level_prefix >= 15) ? level_prefix - 3 : suffixLength;
            if (levelSuffixSize > 0)
                level_suffix = dp->u(levelSuffixSize);
            else
                level_suffix = 0;

            levelCode = (imin(15, level_prefix) << suffixLength) + level_suffix;
            if (level_prefix >= 15 && suffixLength == 0)
                levelCode += 15;
            if (level_prefix >= 16)
                levelCode += (1 << (level_prefix - 3)) - 4096;
            if (i == TotalCoeff - 1 - TrailingOnes && TrailingOnes < 3)
                levelCode += 2;

            if ((levelCode % 2) == 0)
                levelVal[i] = (levelCode + 2) >> 1;
            else
                levelVal[i] = (-levelCode - 1) >> 1;

            if (suffixLength == 0)
                suffixLength = 1;
            if (iabs(levelVal[i]) > (3 << (suffixLength - 1)) && suffixLength < 6)
                suffixLength++;
        }

        int zerosLeft = 0;
        if (TotalCoeff < maxNumCoeff) {
            int yuv = maxNumCoeff == 4 ? 0 : maxNumCoeff == 8 ? 1 : 2;
            zerosLeft = this->parse_total_zeros(yuv, TotalCoeff);
        }

        for (int i = TotalCoeff - 1; i > 0; i--) {
//        for (i = 0; i < TotalCoeff - 1; i++) {
            if (zerosLeft > 0)
                runVal[i] = this->parse_run_before(zerosLeft);
            else
                runVal[i] = 0;
            zerosLeft -= runVal[i];
        }
        runVal[0] = zerosLeft;
//        runVal[TotalCoeff - 1] = zerosLeft;
    }

    *number_coefficients = TotalCoeff;
}


static void read_tc_luma(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    const byte (*pos_scan8x8)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN8x8 : FIELD_SCAN8x8;

    if (IS_I16MB(currMB) && !currMB->dpl_flag) {
        //int16_t i16x16DClevel[16];
        //currMB->residual_block_cavlc(i16x16DClevel, 0, 15, 16, pl, 0, 0);

        int levelVal[16], runVal[16], numcoeff;
        int nC = predict_nnz(currMB, pl, 0 * 4, 0 * 4);
        currMB->read_coeff_4x4_CAVLC(16, nC, levelVal, runVal, &numcoeff);

        int coeffNum = -1;
    //    for (int k = numcoeff - 1; k >= 0; k--) {
        for (int k = 0; k < numcoeff; ++k) {
            if (levelVal[k] != 0) {
                coeffNum += runVal[k] + 1;
                //coeffLevel[startIdx + coeffNum] = levelVal[k];
                int i0 = pos_scan4x4[coeffNum][0];
                int j0 = pos_scan4x4[coeffNum][1];
                currSlice->cof[pl][j0 << 2][i0 << 2] = levelVal[k];
            }
        }

        if (!currMB->TransformBypassModeFlag)
            itrans_2(currMB, pl);
    }

    int qp_per = currMB->qp_scaled[pl] / 6;
    int qp_rem = currMB->qp_scaled[pl] % 6;
    int transform_pl = sps->separate_colour_plane_flag ? currSlice->colour_plane_id : pl;
    int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
        currSlice->InvLevelScale4x4_Intra[transform_pl][qp_rem] :
        currSlice->InvLevelScale4x4_Inter[transform_pl][qp_rem];
    int (*InvLevelScale8x8)[8] = currMB->is_intra_block ?
        currSlice->InvLevelScale8x8_Intra[transform_pl][qp_rem] :
        currSlice->InvLevelScale8x8_Inter[transform_pl][qp_rem];

    int start_scan = IS_I16MB(currMB) ? 1 : 0;

    for (int i8x8 = 0; i8x8 < 4; i8x8++) {
        int block_x = (i8x8 % 2) * 2;
        int block_y = (i8x8 / 2) * 2;

        for (int i4x4 = 0; i4x4 < 4; i4x4++) {
            if (currMB->cbp & (1 << i8x8)) {
                int i = block_x + (i4x4 % 2);
                int j = block_y + (i4x4 / 2);

                int levarr[16] = {0}, runarr[16] = {0}, numcoeff;
                int nC = predict_nnz(currMB, pl, i * 4, j * 4);
                currMB->read_coeff_4x4_CAVLC(16 - start_scan, nC, levarr, runarr, &numcoeff);
                currMB->nz_coeff[pl][j][i] = numcoeff;

                int coef_ctr = start_scan - 1;

                for (int k = 0; k < numcoeff; ++k) {
                //for (int k = numcoeff - 1; k >= 0; k--) {
                    if (levarr[k] != 0) {
                        coef_ctr += runarr[k] + 1;

                        if (!currMB->transform_size_8x8_flag) {
                            currMB->s_cbp[pl].blk |= (int64)(0x01 << (j * 4 + i));
                            int i0 = pos_scan4x4[coef_ctr][0];
                            int j0 = pos_scan4x4[coef_ctr][1];

                            if (!currMB->TransformBypassModeFlag)
                                currSlice->cof[pl][j * 4 + j0][i * 4 + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0]) << qp_per, 4);
                            else
                                currSlice->cof[pl][j * 4 + j0][i * 4 + i0] = levarr[k];
                        } else {
                            currMB->s_cbp[pl].blk |= (int64)(0x33 << (block_y * 4 + block_x));
                            int i0 = pos_scan8x8[coef_ctr * 4 + i4x4][0];
                            int j0 = pos_scan8x8[coef_ctr * 4 + i4x4][1];

                            if (!currMB->TransformBypassModeFlag)
                                currSlice->cof[pl][block_y * 4 + j0][block_x * 4 + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale8x8[j0][i0]) << qp_per, 6);
                            else
                                currSlice->cof[pl][block_y * 4 + j0][block_x * 4 + i0] = levarr[k];
                        }
                    }
                }
            }
        }
        if (!(currMB->cbp & (1 << i8x8))) {
            currMB->nz_coeff[pl][block_y + 0][block_x + 0] = 0;
            currMB->nz_coeff[pl][block_y + 0][block_x + 1] = 0;
            currMB->nz_coeff[pl][block_y + 1][block_x + 0] = 0;
            currMB->nz_coeff[pl][block_y + 1][block_x + 1] = 0;
        }
    }
}

static void read_tc_chroma(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int NumC8x8 = 4 / (sps->SubWidthC * sps->SubHeightC);

    if (currMB->cbp > 15) {
        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            int levelVal[16], runVal[16], numcoeff;
            int nC = sps->ChromaArrayType == 1 ? -1 : sps->ChromaArrayType == 2 ? -2 : 0;
            currMB->read_coeff_4x4_CAVLC(NumC8x8 * 4, nC, levelVal, runVal, &numcoeff);

            int coeffNum = -1;
            for (int k = 0; k < numcoeff; ++k) {
            //for (int k = numcoeff - 1; k >= 0; k--) {
                int i0, j0;
                if (levelVal[k] != 0) {
                    coeffNum += runVal[k] + 1;
                    if (sps->ChromaArrayType == 1) {
                        i0 = coeffNum % 2;
                        j0 = coeffNum / 2;
                        currMB->s_cbp[0].blk |= (int64)(0xf << (iCbCr * 4 + 16));
                    }
                    if (sps->ChromaArrayType == 2) {
                        i0 = FIELD_SCAN[coeffNum][0];
                        j0 = FIELD_SCAN[coeffNum][1];
                        currMB->s_cbp[0].blk |= (int64)(0xff << (iCbCr * 8 + 16));
                    }
                    currSlice->cof[iCbCr + 1][j0 * 4][i0 * 4] = levelVal[k];
                }
            }

            if (sps->ChromaArrayType == 1) {
                int smb = (currSlice->slice_type == SP_SLICE && !currMB->is_intra_block) ||
                          (currSlice->slice_type == SI_SLICE && currMB->mb_type == SI4MB);
                if (!smb && !currMB->TransformBypassModeFlag)
                    itrans_420(currMB, (ColorPlane)(iCbCr + 1));
            }
            if (sps->ChromaArrayType == 2) {
                if (!currMB->TransformBypassModeFlag)
                    itrans_422(currMB, (ColorPlane)(iCbCr + 1));
            }
        }
    }

    if (currMB->cbp > 31) {
        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
            int qp_per_uv = currMB->qp_scaled[iCbCr + 1] / 6;
            int qp_rem_uv = currMB->qp_scaled[iCbCr + 1] % 6;

            for (int i8x8 = 0; i8x8 < NumC8x8; i8x8++) {
                int (*InvLevelScale4x4)[4] = NULL;
                if (!currMB->TransformBypassModeFlag)
                    InvLevelScale4x4 = currMB->is_intra_block ?
                        currSlice->InvLevelScale4x4_Intra[iCbCr + 1][qp_rem_uv] :
                        currSlice->InvLevelScale4x4_Inter[iCbCr + 1][qp_rem_uv];

                for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                    int i = (i4x4 % 2);
                    int j = (i4x4 / 2) + (i8x8 * 2);

                    int levarr[16], runarr[16], numcoeff;
                    int nC = predict_nnz(currMB, iCbCr + 1, i * 4, j * 4);
                    currMB->read_coeff_4x4_CAVLC(15, nC, levarr, runarr, &numcoeff);
                    currMB->nz_coeff[iCbCr + 1][j][i] = numcoeff;

                    int coef_ctr = 0;
                    //for (int k = numcoeff - 1; k >= 0; k--) {
                    for (int k = 0; k < numcoeff; ++k) {
                        if (levarr[k] != 0) {
                            currMB->s_cbp[0].blk |= (int64)(0x1 << (i8x8 * 4 + i4x4 + 16));
                            coef_ctr += runarr[k] + 1;
                            int i0 = pos_scan4x4[coef_ctr][0];
                            int j0 = pos_scan4x4[coef_ctr][1];

                            if (!currMB->TransformBypassModeFlag)
                                currSlice->cof[iCbCr + 1][j * 4 + j0][i * 4 + i0] = rshift_rnd_sf((levarr[k] * InvLevelScale4x4[j0][i0]) << qp_per_uv, 4);
                            else
                                currSlice->cof[iCbCr + 1][j * 4 + j0][i * 4 + i0] = levarr[k];
                        }
                    }
                }
            }
        }
    } else {
        memset(currMB->nz_coeff[1][0], 0, 2 * 16 * sizeof(uint8_t));
    }
}


void macroblock_t::read_CBP_and_coeffs_from_NAL_CAVLC()
{
    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;

    read_tc_luma(this, PLANE_Y);
    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422)
        read_tc_chroma(this);
    if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
        read_tc_luma(this, PLANE_U);
        read_tc_luma(this, PLANE_V);
    }
}

void macroblock_t::residual_block_cavlc(int16_t coeffLevel[16], uint8_t startIdx, uint8_t endIdx,
                                        uint8_t maxNumCoeff, ColorPlane pl, int bx, int by)
{
    slice_t *slice = this->p_Slice;

    const byte (*pos_scan4x4)[2] = !slice->field_pic_flag && !this->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;

    int levelVal[16], runVal[16], numcoeff;
    int nC = predict_nnz(this, pl, bx * 4, by * 4);
    //this->nz_coeff[pl][by][bx] = 0;
    this->read_coeff_4x4_CAVLC(16, nC, levelVal, runVal, &numcoeff);
    this->nz_coeff[pl][by][bx] = numcoeff;

    int coeffNum = -1;
//    for (int k = numcoeff - 1; k >= 0; k--) {
    for (int k = 0; k < numcoeff; ++k) {
        if (levelVal[k] != 0) {
            coeffNum += runVal[k] + 1;
            coeffLevel[startIdx + coeffNum] = levelVal[k];
            int i0 = pos_scan4x4[coeffNum][0];
            int j0 = pos_scan4x4[coeffNum][1];
            slice->cof[pl][j0 << 2][i0 << 2] = levelVal[k];
        }
    }
}
