#include <functional>

#include "global.h"
#include "slice.h"
#include "bitstream_cabac.h"
#include "data_partition.h"
#include "macroblock.h"
#include "quantization.h"
#include "transform.h"
#include "neighbour.h"


using vio::h264::cabac_context_t;
using vio::h264::cabac_engine_t;


#define IS_I16MB(MB) ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)


// Table 9-5 coeff_token mapping to TotalCoeff(coeff_token) and TrailingOnes(coeff_token)
static const uint8_t coeff_token_length[5][4][17] = {
    //  0 <= nC < 2
    {{  1,  6,  8,  9, 10, 11, 13, 13, 13, 14, 14, 15, 15, 16, 16, 16, 16 },
     {  0,  2,  6,  8,  9, 10, 11, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16 },
     {  0,  0,  3,  7,  8,  9, 10, 11, 13, 13, 14, 14, 15, 15, 16, 16, 16 },
     {  0,  0,  0,  5,  6,  7,  8,  9, 10, 11, 13, 14, 14, 15, 15, 16, 16 }},
    // 2 <= nC < 4
    {{  2,  6,  6,  7,  8,  8,  9, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14 },
     {  0,  2,  5,  6,  6,  7,  8,  9, 11, 11, 12, 12, 13, 13, 14, 14, 14 },
     {  0,  0,  3,  6,  6,  7,  8,  9, 11, 11, 12, 12, 13, 13, 13, 14, 14 },
     {  0,  0,  0,  4,  4,  5,  6,  6,  7,  9, 11, 11, 12, 13, 13, 13, 14 }},
    // 4 <= nC < 8
    {{  4,  6,  6,  6,  7,  7,  7,  7,  8,  8,  9,  9,  9, 10, 10, 10, 10 },
     {  0,  4,  5,  5,  5,  5,  6,  6,  7,  8,  8,  9,  9,  9, 10, 10, 10 },
     {  0,  0,  4,  5,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9, 10, 10, 10 },
     {  0,  0,  0,  4,  4,  4,  4,  4,  5,  6,  7,  8,  8,  9, 10, 10, 10 }},
    // nC == -1 
    {{  2,  6,  6,  6,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  1,  6,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  3,  7,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  6,  7,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }},
    // nC == -2 
    {{  1,  7,  7,  9,  9, 10, 11, 12, 13,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  2,  7,  7,  9, 10, 11, 12, 12,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  3,  7,  7,  9, 10, 11, 12,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  5,  6,  7,  7, 10, 11,  0,  0,  0,  0,  0,  0,  0,  0 }}
};

static const uint8_t coeff_token_code[5][4][17] = {
    // 0 <= nC < 2
    {{  1,  5,  7,  7,  7,  7, 15, 11,  8, 15, 11, 15, 11, 15, 11,  7,  4 },
     {  0,  1,  4,  6,  6,  6,  6, 14, 10, 14, 10, 14, 10,  1, 14, 10,  6 },
     {  0,  0,  1,  5,  5,  5,  5,  5, 13,  9, 13,  9, 13,  9, 13,  9,  5 },
     {  0,  0,  0,  3,  3,  4,  4,  4,  4,  4, 12, 12,  8, 12,  8, 12,  8 }},
    // 2 <= nC < 4
    {{  3, 11,  7,  7,  7,  4,  7, 15, 11, 15, 11,  8, 15, 11,  7,  9,  7 },
     {  0,  2,  7, 10,  6,  6,  6,  6, 14, 10, 14, 10, 14, 10, 11,  8,  6 },
     {  0,  0,  3,  9,  5,  5,  5,  5, 13,  9, 13,  9, 13,  9,  6, 10,  5 },
     {  0,  0,  0,  5,  4,  6,  8,  4,  4,  4, 12,  8, 12, 12,  8,  1,  4 }},
    // 4 <= nC < 8
    {{ 15, 15, 11,  8, 15, 11,  9,  8, 15, 11, 15, 11,  8, 13,  9,  5,  1 },
     {  0, 14, 15, 12, 10,  8, 14, 10, 14, 14, 10, 14, 10,  7, 12,  8,  4 },
     {  0,  0, 13, 14, 11,  9, 13,  9, 13, 10, 13,  9, 13,  9, 11,  7,  3 },
     {  0,  0,  0, 12, 11, 10,  9,  8, 13, 12, 12, 12,  8, 12, 10,  6,  2 }},
    // nC == -1
    {{  1,  7,  4,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  1,  6,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  1,  2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }},
    // nC == -2
    {{  1, 15, 14,  7,  6,  7,  7,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  1, 13, 12,  5,  6,  6,  6,  5,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  1, 11, 10,  4,  5,  5,  4,  0,  0,  0,  0,  0,  0,  0,  0 },
     {  0,  0,  0,  1,  1,  9,  8,  4,  4,  0,  0,  0,  0,  0,  0,  0,  0 }}
};

// Table 9-7 total_zeros tables for 4x4 blocks with tzVlcIndex 1 to 7
// Table 9-8 total_zeros tables for 4x4 blocks with tzVlcIndex 8 to 15
// Table 9-9 total_zeros tables for chroma DC 2x2 and 2x4 blocks
static const uint8_t total_zeros_length[3][15][16] = {
    // YUV420
    {{ 1, 2, 3, 3 },
     { 1, 2, 2 },
     { 1, 1 }},
    // YUV422
    {{ 1, 3, 3, 4, 4, 4, 5, 5 },
     { 3, 2, 3, 3, 3, 3, 3 },
     { 3, 3, 2, 2, 3, 3 },
     { 3, 2, 2, 2, 3 },
     { 2, 2, 2, 2 },
     { 2, 2, 1 },
     { 1, 1 }},
    // YUV444
    {{ 1, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 9 },
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
     { 1, 1 }}
};

static const uint8_t total_zeros_code[3][15][16] = {
    // YUV420
    {{ 1, 1 , 1 , 0 },
     { 1, 1 , 0 },
     { 1, 0 }},
    // YUV422
    {{ 1, 2, 3, 2, 3, 1, 1, 0 },
     { 0, 1, 1, 4, 5, 6, 7 },
     { 0, 1, 1, 2, 6, 7 },
     { 6, 0, 1, 2, 7 },
     { 0, 1, 2, 3 },
     { 0, 1, 1 },
     { 0, 1 }},
    // YUV444
    {{ 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 1 },
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
     { 0, 1 }}
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

    int tab = min<int>(zerosLeft, 7) - 1;

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

            levelCode = (min(15, level_prefix) << suffixLength) + level_suffix;
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
            if (abs(levelVal[i]) > (3 << (suffixLength - 1)) && suffixLength < 6)
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


typedef enum {
    LUMA_16DC     =  0, // ctxBlockCat =  0
    LUMA_16AC     =  1, // ctxBlockCat =  1
    LUMA_4x4      =  2, // ctxBlockCat =  2
    CHROMA_DC     =  3, // ctxBlockCat =  3
    CHROMA_AC     =  4, // ctxBlockCat =  4
    LUMA_8x8      =  5, // ctxBlockCat =  5 =
    CHROMA_DC_2x4 =  6, // ctxBlockCat =
    CB_16DC       =  7, // ctxBlockCat =  6
    CB_16AC       =  8, // ctxBlockCat =  7
    CB_4x4        =  9, // ctxBlockCat =  8
    CB_8x8        = 10, // ctxBlockCat =  9 =
    CR_16DC       = 11, // ctxBlockCat = 10
    CR_16AC       = 12, // ctxBlockCat = 11
    CR_4x4        = 13, // ctxBlockCat = 12
    CR_8x8        = 14  // ctxBlockCat = 13 =
} CABACBlockTypes;

static int read_and_store_CBP_block_bit(mb_t* mb, int type)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int y_dc = (type==LUMA_16DC || type==CB_16DC || type==CR_16DC); 
    int u_dc = ((type==CHROMA_DC || type==CHROMA_DC_2x4) && !mb->is_v_block);
    int v_dc = ((type==CHROMA_DC || type==CHROMA_DC_2x4) &&  mb->is_v_block);
    int y_ac = (type==LUMA_16AC || type==LUMA_8x8 || type==LUMA_4x4) ||
               (type==CB_16AC || type==CB_8x8 || type==CB_4x4) ||
               (type==CR_16AC || type==CR_8x8 || type==CR_4x4);
    int u_ac = (type==CHROMA_AC && !mb->is_v_block);
    int v_ac = (type==CHROMA_AC &&  mb->is_v_block);

    int size_8x8_flag = (type == LUMA_8x8 || type == CB_8x8 || type == CR_8x8);
    int pl = (type == CB_8x8 || type == CB_4x4 || type == CB_16AC || type == CB_16DC) ? 1 :
             (type == CR_8x8 || type == CR_4x4 || type == CR_16AC || type == CR_16DC) ? 2 : 0;

    int i = (y_ac || u_ac || v_ac ? mb->subblock_x : 0);
    int j = (y_ac || u_ac || v_ac ? mb->subblock_y : 0);
    int bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35);

    int bit_pos_a = 0;
    int bit_pos_b = 0;

    int mb_size[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };
    PixelPos block_a, block_b;
    get4x4Neighbour(mb, i - 1, j, mb_size[y_dc || y_ac ? IS_LUMA : IS_CHROMA], &block_a);
    get4x4Neighbour(mb, i, j - 1, mb_size[y_dc || y_ac ? IS_LUMA : IS_CHROMA], &block_b);
    if (y_ac || u_ac || v_ac) {
        if (block_a.available)
            bit_pos_a = 4 * block_a.y + block_a.x;
        if (block_b.available)
            bit_pos_b = 4 * block_b.y + block_b.x;
    }

    int condTermFlagA = (mb->is_intra_block ? 1 : 0);
    int condTermFlagB = (mb->is_intra_block ? 1 : 0);
    if ((sps->separate_colour_plane_flag || sps->chroma_format_idc != YUV444) && type != LUMA_8x8) {
        if (block_a.available) {
            mb_t* mb_a = &slice->mb_data[block_a.mb_addr];
            if (mb_a->mb_type == IPCM)
                condTermFlagA = 1;
            else
                condTermFlagA = (mb_a->s_cbp[0].bits >> (bit + bit_pos_a)) & 1;
        }
        if (block_b.available) {
            mb_t* mb_b = &slice->mb_data[block_b.mb_addr];
            if (mb_b->mb_type == IPCM)
                condTermFlagB = 1;
            else
                condTermFlagB = (mb_b->s_cbp[0].bits >> (bit + bit_pos_b)) & 1;
        }
    } else if (sps->chroma_format_idc == YUV444) {
        if (block_a.available) {
            mb_t* mb_a = &slice->mb_data[block_a.mb_addr];
            if (!size_8x8_flag || mb_a->transform_size_8x8_flag) {
                if (mb_a->mb_type == IPCM)
                    condTermFlagA = 1;
                else
                    condTermFlagA = (mb_a->s_cbp[pl].bits >> (bit + bit_pos_a)) & 1;
            }
        }
        if (block_b.available) {
            mb_t* mb_b = &slice->mb_data[block_b.mb_addr];
            if (!size_8x8_flag || mb_b->transform_size_8x8_flag) {
                if (mb_b->mb_type == IPCM)
                    condTermFlagB = 1;
                else
                    condTermFlagB = (mb_b->s_cbp[pl].bits >> (bit + bit_pos_b)) & 1;
            }
        }
    }
    int ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    return ctxIdxInc;
}

static void update_cbp(mb_t* mb, int type)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int y_dc = (type==LUMA_16DC || type==CB_16DC || type==CR_16DC); 
    int u_dc = ((type==CHROMA_DC || type==CHROMA_DC_2x4) && !mb->is_v_block);
    int v_dc = ((type==CHROMA_DC || type==CHROMA_DC_2x4) &&  mb->is_v_block);
    int y_ac = (type==LUMA_16AC || type==LUMA_8x8 || type==LUMA_4x4) ||
               (type==CB_16AC || type==CB_8x8 || type==CB_4x4) ||
               (type==CR_16AC || type==CR_8x8 || type==CR_4x4);
    int u_ac = (type==CHROMA_AC && !mb->is_v_block);
    int v_ac = (type==CHROMA_AC &&  mb->is_v_block);

    int pl = (type == CB_8x8 || type == CB_4x4 || type == CB_16AC || type == CB_16DC) ? 1 :
             (type == CR_8x8 || type == CR_4x4 || type == CR_16AC || type == CR_16DC) ? 2 : 0;
    int temp_pl = sps->chroma_format_idc == YUV444 ? pl : 0;

    int i = (y_ac || u_ac || v_ac ? mb->subblock_x : 0);
    int j = (y_ac || u_ac || v_ac ? mb->subblock_y : 0);
    int bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35)
            + (y_ac || u_ac || v_ac ? j + i / 4 : 0);

    int cbp = (type == LUMA_8x8 || type == CB_8x8 || type == CR_8x8) ? 0x33 : 0x01;

    mb->s_cbp[temp_pl].bits |= ((int64_t)cbp << bit);
}

static uint32_t unary_exp_golomb_level_decode(cabac_engine_t* dep_dp, cabac_context_t* ctx)
{
    const uint32_t cMax = 13;

    uint32_t bins = -1;
    bool b;
    for (b = 1; b && (bins + 1 < cMax); ++bins)
        b = dep_dp->decode_decision(ctx);
    if (!b)
        return bins;

    uint32_t k = 0;
    while (dep_dp->decode_bypass())
        bins += (1 << k++);
    while (k--)
        bins += (dep_dp->decode_bypass() << k);

    return bins + 1;
}

struct syntax_element_t {
    int value1;  //!< numerical value of syntax element
    int value2;  //!< for blocked symbols, e.g. run/level
    int context; //!< CABAC context
};

static void readRunLevel_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    static const short maxpos[22] = {
        15, 14, 15,  3, 14, 63,  7,
        15, 14, 15, 63, 15, 14, 15, 63
    };
    static const short c1isdc[22] = {
         1,  0,  1,  1,  0,  1,  1,
         1,  0,  1,  1,  1,  0,  1,  1
    };
    static const short max_c2[22] = {
         4,  4,  4,  3,  4,  4,  3,
         4,  4,  4,  4,  4,  4,  4,  4
    }; // 9

    static const short type2ctx_bcbp[22] = {
         0,  1,  4,  5,  6,  2,  5,
        10, 11, 14, 12, 16, 17, 20, 18
    };

    static const short type2ctx_map[22] = {
         0,  1,  5,  6,  7,  2,  6,
        10, 11, 15, 12, 16, 17, 21, 18
    }; // 8
    static const short type2ctx_one[22] = {
         0,  1,  4,  5,  6,  2,  5,
        10, 11, 14, 12, 16, 17, 20, 18
    }; // 7

    //===== position -> ctx for MAP =====
    //--- zig-zag scan ----
    static const uint8_t pos2ctx_map8x8[] = {
         0,  1,  2,  3,  4,  5,  5,  4,  4,  3,  3,  4,  4,  4,  5,  5,
         4,  4,  4,  4,  3,  3,  6,  7,  7,  7,  8,  9, 10,  9,  8,  7,
         7,  6, 11, 12, 13, 11,  6,  7,  8,  9, 14, 10,  9,  8,  6, 11,
        12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map8x8i[] = {
        0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  7,  7,  7,  8,  4,  5,
        6,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 11, 12, 11,
        9,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 13, 13,  9,
        9, 10, 10,  8, 13, 13,  9,  9, 10, 10, 14, 14, 14, 14, 14, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map4x4[] = {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map2x4c[] = {
        0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2
    }; // 15 CTX

    static const uint8_t *pos2ctx_map[2][22] = {
        { pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map4x4,  
          pos2ctx_map4x4,  pos2ctx_map8x8,
          pos2ctx_map2x4c,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map8x8,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map8x8 },
        { pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map8x8i,
          pos2ctx_map2x4c,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map8x8i,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map8x8i }
    };

    //===== position -> ctx for LAST =====
    static const uint8_t pos2ctx_last8x8[] = {
        0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
        5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8
    }; //  9 CTX
    static const uint8_t pos2ctx_last4x4[] = {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
    }; // 15 CTX
    static const uint8_t pos2ctx_last2x4c[] = {
        0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2
    }; // 15 CTX

    static const uint8_t *pos2ctx_last[22] = {
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last4x4,  pos2ctx_last8x8,  
        pos2ctx_last2x4c,
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last4x4,  pos2ctx_last8x8,
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last4x4,  pos2ctx_last8x8
    };

    slice_t *currSlice = currMB->p_Slice;
    sps_t* sps = currSlice->active_sps;

    bool field = currSlice->field_pic_flag || currMB->mb_field_decoding_flag;
    const byte *pos2ctx_Map  = pos2ctx_map [field][se->context];
    const byte *pos2ctx_Last = pos2ctx_last[se->context];
    cabac_context_t *map_ctx  = currSlice->mot_ctx->map_contexts [field][type2ctx_map[se->context]];
    cabac_context_t *last_ctx = currSlice->mot_ctx->last_contexts[field][type2ctx_map[se->context]];
    cabac_context_t *one_ctx = currSlice->mot_ctx->one_contexts[type2ctx_one[se->context]];
    cabac_context_t *abs_ctx = currSlice->mot_ctx->abs_contexts[type2ctx_one[se->context]];
    const short max_type = max_c2[se->context];

    static int coeff_pos = 0;
    static int coded_block_flag = -1;
    static int coeff_val[64];

    if (coded_block_flag < 0) {
        coded_block_flag = 1; // always one for 8x8 mode
        if (sps->chroma_format_idc == YUV444 || se->context != LUMA_8x8) {
            cabac_context_t* ctx = currSlice->mot_ctx->bcbp_contexts[type2ctx_bcbp[se->context]];
            int ctxIdxInc = read_and_store_CBP_block_bit(currMB, se->context);

            coded_block_flag = dep_dp->decode_decision(ctx + ctxIdxInc);
        }
        if (coded_block_flag)
            update_cbp(currMB, se->context);

        if (coded_block_flag) {
            int i;
            int *coeff = coeff_val;
            int i0     = 0;
            int i1     = maxpos[se->context];
            coded_block_flag = 0;

            if (!c1isdc[se->context]) {
                ++i0;
                ++i1;
            }

            for (i = i0; i < i1; ++i) {
                bool significant_coeff_flag = dep_dp->decode_decision(map_ctx + pos2ctx_Map[i]);
                *(coeff++) = significant_coeff_flag;
                coded_block_flag += significant_coeff_flag;
                if (significant_coeff_flag) {
                    bool last_significant_coeff_flag = dep_dp->decode_decision(last_ctx + pos2ctx_Last[i]);
                    if (last_significant_coeff_flag) {
                        memset(coeff, 0, (i1 - i) * sizeof(int));
                        i = i1 + 1;
                        break;
                    }
                }
            }
            if (i <= i1) {
                bool significant_coeff_flag = 1;
                *(coeff++) = significant_coeff_flag;
                coded_block_flag += significant_coeff_flag;
            }

            i = maxpos[se->context];
            int *cof = coeff_val + i;
            int c1 = 1;
            int c2 = 0;

            for (; i >= 0; i--) {
                if (*cof) {
                    if (dep_dp->decode_decision(one_ctx + c1)) {
                        *cof += unary_exp_golomb_level_decode(dep_dp, abs_ctx + c2) + 1;
                        c2 = min<int>(++c2, max_type);
                        c1 = 0;
                    } else if (c1)
                        c1 = min<int>(++c1, 4);

                    if (dep_dp->decode_bypass())
                        *cof = - *cof;
                }
                cof--;
            }
        }
    }

    se->value1 = 0;
    se->value2 = 0;
    if (coded_block_flag) {
        while (coeff_val[coeff_pos] == 0) {
            ++coeff_pos;
            ++se->value2;
        }
        se->value1 = coeff_val[coeff_pos++];
    }
    if ((coded_block_flag)-- == 0) 
        coeff_pos = 0;
}


// residual_block_cavlc(coeffLevel, startIdx, endIdx, maxNumCoeff)
void macroblock_t::residual_block_cavlc(int16_t coeffLevel[16], uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                        ColorPlane pl, bool chroma, bool ac, int blkIdx)
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    int i = chroma ? blkIdx % 2 : ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int j = chroma ? blkIdx / 2 : ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    int start_scan = 0;
    if (!chroma && ac)
        start_scan = IS_I16MB(this) ? 1 : 0;

    int levelVal[16], runVal[16], numcoeff;
    int nC;
    if (chroma && !ac)
        nC = sps->ChromaArrayType == 1 ? -1 : sps->ChromaArrayType == 2 ? -2 : 0;
    else
        nC = predict_nnz(this, pl, i * 4, j * 4);

    this->read_coeff_4x4_CAVLC(maxNumCoeff - start_scan, nC, levelVal, runVal, &numcoeff);

    if (ac)
        this->nz_coeff[pl][j][i] = numcoeff;

    int coeffNum = -1;
    if (!chroma && ac)
        coeffNum = start_scan - 1;
    if (chroma && ac)
        coeffNum = 0;
    //for (int k = numcoeff - 1; k >= 0; k--) {
    for (int k = 0; k < numcoeff; ++k) {
        if (levelVal[k] != 0) {
            coeffNum += runVal[k] + 1;
            //coeffLevel[start_scan + coeffNum] = levelVal[k];
            if (!chroma) {
                if (!ac)
                    quantization.coeff_luma_dc(this, pl, i, j, coeffNum, levelVal[k]);
                else {
                    int x0 = !this->transform_size_8x8_flag ? i : (i & ~1);
                    int y0 = !this->transform_size_8x8_flag ? j : (j & ~1);
                    int c0 = !this->transform_size_8x8_flag ? coeffNum : coeffNum * 4 + (blkIdx % 4);
                    quantization.coeff_luma_ac(this, pl, x0, y0, c0, levelVal[k]);
                }
            } else {
                if (!ac)
                    quantization.coeff_chroma_dc(this, pl, i, j, coeffNum, levelVal[k]);
                else
                    quantization.coeff_chroma_ac(this, pl, i, j, coeffNum, levelVal[k]);
            }
        }
    }
}

// residual_block_cabac(coeffLevel, startIdx, endIdx, maxNumCoeff)
void macroblock_t::residual_block_cabac(int16_t coeffLevel[16], uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                        ColorPlane pl, bool chroma, bool ac, int blkIdx)
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    int i = chroma ? blkIdx % 2 : ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int j = chroma ? blkIdx / 2 : ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    data_partition_t* dp = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];
    syntax_element_t currSE;
    if (!chroma) {
        if (!ac)
            currSE.context = LUMA_16DC;
        else {
            if (pl == PLANE_Y || sps->separate_colour_plane_flag)
                currSE.context = this->transform_size_8x8_flag ? LUMA_8x8 : IS_I16MB(this) ? LUMA_16AC : LUMA_4x4;
            else if (pl == PLANE_U)
                currSE.context = this->transform_size_8x8_flag ? CB_8x8 : IS_I16MB(this) ? CB_16AC : CB_4x4;
            else
                currSE.context = this->transform_size_8x8_flag ? CR_8x8 : IS_I16MB(this) ? CR_16AC : CR_4x4;
        }
    } else {
        if (!ac)
            currSE.context = sps->ChromaArrayType == 1 ? CHROMA_DC : CHROMA_DC_2x4;
        else
            currSE.context = CHROMA_AC;
    }
    if (chroma)
        this->is_v_block = pl - 1;
    if (ac) {
        this->subblock_x = i * 4;
        this->subblock_y = j * 4;
    }

    int start_scan = 0;
    if (!chroma && ac)
        start_scan = IS_I16MB(this) ? 1 : 0;

    int coef_ctr = -1;
    if (!chroma && ac)
        coef_ctr = start_scan - 1;
    if (chroma && ac)
        coef_ctr = 0;
    int level = 1;
    for (int k = start_scan; k < maxNumCoeff + 1 && level != 0; ++k) {
        readRunLevel_CABAC(this, &currSE, &dp->de_cabac);
        level = currSE.value1;
        if (level != 0) {
            coef_ctr += currSE.value2 + 1;
            //coeffLevel[startIdx + coef_ctr] = level;
            if (!ac)
                assert(coef_ctr < maxNumCoeff);
            if (!chroma) {
                if (!ac)
                    quantization.coeff_luma_dc(this, pl, i, j, coef_ctr, level);
                else
                    quantization.coeff_luma_ac(this, pl, i, j, coef_ctr, level);
            } else {
                if (!ac)
                    quantization.coeff_chroma_dc(this, pl, i, j, coef_ctr, level);
                else
                    quantization.coeff_chroma_ac(this, pl, i, j, coef_ctr, level);
            }
        }
    }
}

// residual_luma(i16x16DClevel, i16x16AClevel, level4x4, level8x8, 0, 15)
void macroblock_t::residual_luma(ColorPlane pl)
{
    slice_t* slice = this->p_Slice;
    pps_t* pps = slice->active_pps;

    int CodedBlockPatternLuma = this->cbp & 15;

    auto residual_block = !pps->entropy_coding_mode_flag ?
        std::mem_fn(&macroblock_t::residual_block_cavlc) :
        std::mem_fn(&macroblock_t::residual_block_cabac);

    if (IS_I16MB(this) && !this->dpl_flag) {
        // residual_block(i16x16DClevel, 0, 15, 16);
        residual_block(this, nullptr, 0, 15, 16, pl, false, false, 0);

        transform.inverse_luma_dc(this, pl);
    }

    for (int i8x8 = 0; i8x8 < 4; i8x8++) {
        if (!this->transform_size_8x8_flag || !pps->entropy_coding_mode_flag) {
            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (CodedBlockPatternLuma & (1 << i8x8)) {
                    //if (IS_I16MB(this))
                        // residual_block(i16x16AClevel[i8x8 * 4 + i4x4], max(0, startIdx - 1), endIdx - 1, 15);
                    //else
                        // residual_block(level4x4[i8x8 * 4 + i4x4], startIdx, endIdx, 16);
                    residual_block(this, nullptr, 0, 0, 16, pl, false, true, i8x8 * 4 + i4x4);
                } else {
                    if (!pps->entropy_coding_mode_flag) {
                        int i = (i8x8 % 2) * 2 + (i4x4 % 2);
                        int j = (i8x8 / 2) * 2 + (i4x4 / 2);
                        this->nz_coeff[pl][j][i] = 0;
                    }
                }
            }
        } else if (CodedBlockPatternLuma & (1 << i8x8))
            // residual_block(level8x8[i8x8], 4 * startIdx, 4 * endIdx + 3, 64);
            residual_block(this, nullptr, 0, 0, 64, pl, false, true, i8x8 * 4);
        else {
            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (!pps->entropy_coding_mode_flag) {
                    int i = (i8x8 % 2) * 2 + (i4x4 % 2);
                    int j = (i8x8 / 2) * 2 + (i4x4 / 2);
                    this->nz_coeff[pl][j][i] = 0;
                }
            }
        }
    }
}

void macroblock_t::residual_chroma()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    int NumC8x8 = 4 / (sps->SubWidthC * sps->SubHeightC);

    int CodedBlockPatternChroma = this->cbp / 16;

    auto residual_block = !pps->entropy_coding_mode_flag ?
        std::mem_fn(&macroblock_t::residual_block_cavlc) :
        std::mem_fn(&macroblock_t::residual_block_cabac);

    if (CodedBlockPatternChroma & 3) {
        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            // residual_block(ChromaDCLevel[iCbCr], 0, 4 * NumC8x8 - 1, 4 * NumC8x8);
            residual_block(this, nullptr, 0, 4 * NumC8x8 - 1, 4 * NumC8x8, (ColorPlane)(iCbCr+1), true, false, 0);

            transform.inverse_chroma_dc(this, (ColorPlane)(iCbCr + 1));
        }
    }

    for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
        for (int i8x8 = 0; i8x8 < NumC8x8; i8x8++) {
            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (CodedBlockPatternChroma & 2) {
                    // residual_block(ChromaACLevel[iCbCr][i8x8 * 4 + i4x4], max(0, startIdx - 1), endIdx - 1, 15);
                    residual_block(this, nullptr, 0, 0, 15, (ColorPlane)(iCbCr+1), true, true, i8x8 * 4 + i4x4);
                } else {
                    if (!pps->entropy_coding_mode_flag) {
                        int i = (i4x4 % 2);
                        int j = (i4x4 / 2) + (i8x8 * 2);
                        this->nz_coeff[iCbCr + 1][j][i] = 0;
                    }
                }
            }
        }
    }
}

void macroblock_t::residual()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    this->residual_luma(PLANE_Y);
    if (sps->ChromaArrayType == 1 || sps->ChromaArrayType == 2)
        this->residual_chroma();
    else if (sps->ChromaArrayType == 3) {
        this->residual_luma(PLANE_U);
        this->residual_luma(PLANE_V);
    }
}
