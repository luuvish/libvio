#include <functional>

#include "global.h"
#include "slice.h"
#include "macroblock.h"

#include "neighbour.h"


namespace vio  {
namespace h264 {


Parser::Residual::Residual(mb_t& _mb) :
    sps { *_mb.p_Slice->active_sps },
    pps { *_mb.p_Slice->active_pps },
    slice { *_mb.p_Slice },
    mb { _mb },
    se { _mb }
{
}

Parser::Residual::~Residual()
{
}


/*
static int16_t parse_level(InterpreterRbsp *currStream, uint8_t level_prefix, uint8_t suffixLength)
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

void Parser::Residual::residual_block_cavlc(uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                            ColorPlane pl, bool chroma, bool ac, int blkIdx)
{
    InterpreterRbsp* dp = &slice.parser.partArr[slice.parser.dp_mode ? (mb.is_intra_block ? 1 : 2) : 0];

    int i = chroma ? blkIdx % 2 : ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int j = chroma ? blkIdx / 2 : ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    int levelVal[16], runVal[16];
    int nC;
    if (chroma && !ac)
        nC = sps.ChromaArrayType == 1 ? -1 : sps.ChromaArrayType == 2 ? -2 : 0;
    else
        nC = slice.neighbour.predict_nnz(&mb, pl, i * 4, j * 4);

    uint8_t coeff_token  = se.coeff_token(nC);
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
            zerosLeft = se.total_zeros(yuv, TotalCoeff);
        }

        for (int i = TotalCoeff - 1; i > 0; i--) {
//        for (i = 0; i < TotalCoeff - 1; i++) {
            if (zerosLeft > 0)
                runVal[i] = se.run_before(zerosLeft);
            else
                runVal[i] = 0;
            zerosLeft -= runVal[i];
        }
        runVal[0] = zerosLeft;
//        runVal[TotalCoeff - 1] = zerosLeft;
    }

    if (ac)
        mb.nz_coeff[pl][j][i] = TotalCoeff;

    int coeffNum = startIdx - 1;
    //for (int k = TotalCoeff - 1; k >= 0; k--) {
    for (int k = 0; k < TotalCoeff; ++k) {
        if (levelVal[k] != 0) {
            coeffNum += runVal[k] + 1;
            //coeffLevel[start_scan + coeffNum] = levelVal[k];
            if (!chroma) {
                if (!ac)
                    slice.decoder.coeff_luma_dc(&mb, pl, i, j, coeffNum, levelVal[k]);
                else {
                    int x0 = !mb.transform_size_8x8_flag ? i : (i & ~1);
                    int y0 = !mb.transform_size_8x8_flag ? j : (j & ~1);
                    int c0 = !mb.transform_size_8x8_flag ? coeffNum : coeffNum * 4 + (blkIdx % 4);
                    slice.decoder.coeff_luma_ac(&mb, pl, x0, y0, c0, levelVal[k]);
                }
            } else {
                if (!ac)
                    slice.decoder.coeff_chroma_dc(&mb, pl, i, j, coeffNum, levelVal[k]);
                else
                    slice.decoder.coeff_chroma_ac(&mb, pl, i, j, coeffNum, levelVal[k]);
            }
        }
    }
}




// Table 9-43 Mapping of scanning position to ctxIdxInc for ctxBlockCat == 5, 9, or 13

static const uint8_t pos2ctx_map8x8[] = {
     0,  1,  2,  3,  4,  5,  5,  4,  4,  3,  3,  4,  4,  4,  5,  5,
     4,  4,  4,  4,  3,  3,  6,  7,  7,  7,  8,  9, 10,  9,  8,  7,
     7,  6, 11, 12, 13, 11,  6,  7,  8,  9, 14, 10,  9,  8,  6, 11,
    12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 14
};
static const uint8_t pos2ctx_map8x8i[] = {
    0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  7,  7,  7,  8,  4,  5,
    6,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 11, 12, 11,
    9,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 13, 13,  9,
    9, 10, 10,  8, 13, 13,  9,  9, 10, 10, 14, 14, 14, 14, 14, 14
};
static const uint8_t pos2ctx_last8x8[] = {
    0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
    3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
    5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8
};

static const uint8_t pos2ctx_map4x4[] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14
}; // 15 CTX
static const uint8_t pos2ctx_map2x4c[] = {
    0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2
}; // 15 CTX
//===== position -> ctx for LAST =====
static const uint8_t pos2ctx_last4x4[] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
}; // 15 CTX
static const uint8_t pos2ctx_last2x4c[] = {
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

// Table 9-40 Assignment of ctxIdxBlockCatOffset to ctxBlockCat for syntax elements coded_block_flag,
//            significant_coeff_flag, last_significant_coeff_flag, and coeff_abs_level_minus1

static const uint8_t ctxBlockCat[2][4][14] = {
    {{  0,  0,  0,  0,  0, 44, 20, 20, 20, 44, 32, 32, 32, 44 },  // coded_block_flag
     {  0,  0,  0,  0,  0, 61, 76, 76, 76,164,120,120,120,179 },  // significant_coeff_flag
     {  0,  0,  0,  0,  0, 61, 76, 76, 76,164,120,120,120,179 },  // last_significant_coeff_flag
     {  0,  0,  0,  0,  0, 44, 20, 20, 20, 44, 32, 32, 32, 44 }}, // coeff_abs_level_minus1
    {{  0,  4,  8, 12, 16,  0,  0,  4,  8,  4,  0,  4,  8,  8 },  // coded_block_flag
     {  0, 15, 29, 44, 47,  0,  0, 15, 29,  0,  0, 15, 29,  0 },  // significant_coeff_flag
     {  0, 15, 29, 44, 47,  0,  0, 15, 29,  0,  0, 15, 29,  0 },  // last_significant_coeff_flag
     {  0, 10, 20, 30, 39,  0,  0, 10, 20,  0,  0, 10, 20,  0 }}  // coeff_abs_level_minus1
};

static const short type2ctx_bcbp[22] = {
     0,  4,  8, 12, 16, 44, 12,
    20, 24, 28, 48, 32, 36, 40, 52
};
static const short type2ctx_map[22] = {
      0,  15,  29,  44,  47,  61,  44,
     76,  91, 105, 164, 120, 135, 149, 179
};
static const short type2ctx_one[22] = {
      0,  10,  20,  30,  40,  50,  30,
     80,  90, 100,  60, 110, 120, 130,  70
};

// Table 9-42 Specification of ctxBlockCat for the different blocks

typedef enum {
    LUMA_16DC     =  0, // ctxBlockCat =  0
    LUMA_16AC     =  1, // ctxBlockCat =  1
    LUMA_4x4      =  2, // ctxBlockCat =  2
    CHROMA_DC     =  3, // ctxBlockCat =  3
    CHROMA_AC     =  4, // ctxBlockCat =  4
    LUMA_8x8      =  5, // ctxBlockCat =  5 =
    CHROMA_DC_2x4 =  6, // ctxBlockCat =
    //CB_16DC       =  7, // ctxBlockCat =  6
    CB_16AC       =  8, // ctxBlockCat =  7
    CB_4x4        =  9, // ctxBlockCat =  8
    CB_8x8        = 10, // ctxBlockCat =  9 =
    //CR_16DC       = 11, // ctxBlockCat = 10
    CR_16AC       = 12, // ctxBlockCat = 11
    CR_4x4        = 13, // ctxBlockCat = 12
    CR_8x8        = 14  // ctxBlockCat = 13 =
} CABACBlockTypes;

static uint32_t unary_exp_golomb_level_decode(cabac_engine_t* dep_dp, cabac_context_t* ctx, int* ctxIdxIncs)
{
    const uint32_t cMax = 13;

    if (!dep_dp->decode_decision(ctx + ctxIdxIncs[0]))
        return 0;

    uint32_t bins = -1;
    bool b;
    for (b = 1; b && (bins + 1 < cMax); ++bins)
        b = dep_dp->decode_decision(ctx + ctxIdxIncs[1]);
    if (!b)
        return bins + 1;

    uint32_t k = 0;
    while (dep_dp->decode_bypass())
        bins += (1 << k++);
    while (k--)
        bins += (dep_dp->decode_bypass() << k);

    return bins + 1 + 1;
}


void Parser::Residual::residual_block_cabac(uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                            ColorPlane pl, bool chroma, bool ac, int blkIdx)
{
    shr_t& shr = slice.header;

    cabac_engine_t& cabac = slice.parser.cabac[slice.parser.dp_mode ? (mb.is_intra_block ? 1 : 2) : 0];

    int context;
    if (!chroma) {
        if (!ac)
            context = LUMA_16DC;
        else {
            if (pl == PLANE_Y || sps.separate_colour_plane_flag)
                context = mb.transform_size_8x8_flag ? LUMA_8x8 : mb.mb_type == I_16x16 ? LUMA_16AC : LUMA_4x4;
            else if (pl == PLANE_U)
                context = mb.transform_size_8x8_flag ? CB_8x8 : mb.mb_type == I_16x16 ? CB_16AC : CB_4x4;
            else
                context = mb.transform_size_8x8_flag ? CR_8x8 : mb.mb_type == I_16x16 ? CR_16AC : CR_4x4;
        }
    } else {
        if (!ac)
            context = sps.ChromaArrayType == 1 ? CHROMA_DC : CHROMA_DC_2x4;
        else
            context = CHROMA_AC;
    }

    int coded_block_flag = 1; // always one for 8x8 mode
    if (sps.chroma_format_idc == CHROMA_FORMAT_444 || context != LUMA_8x8) {
        cabac_context_t* ctx = slice.parser.mot_ctx.bcbp_contexts + type2ctx_bcbp[context];
        int ctxIdxInc = coded_block_flag_ctxIdxInc(mb, pl, chroma, ac, blkIdx);

        coded_block_flag = cabac.decode_decision(ctx + ctxIdxInc);
    }
    if (coded_block_flag)
        update_coded_block_flag(&mb, pl, chroma, ac, blkIdx);

    if (!coded_block_flag)
        return;

    bool field = shr.field_pic_flag || mb.mb_field_decoding_flag;
    const uint8_t* pos2ctx_Map  = pos2ctx_map [field][context];
    const uint8_t* pos2ctx_Last = pos2ctx_last[context];
    cabac_context_t* map_ctx  = slice.parser.mot_ctx.map_contexts [field] + type2ctx_map[context];
    cabac_context_t* last_ctx = slice.parser.mot_ctx.last_contexts[field] + type2ctx_map[context];

    int coeff_val[64], *coeff;
    int numCoeff = endIdx + 1;
    int ii = startIdx;
    coeff = coeff_val;// + startIdx;

    ii = 0;
    while (ii < numCoeff - 1) {
        bool significant_coeff_flag = cabac.decode_decision(map_ctx + pos2ctx_Map[ii]);
        *(coeff++) = significant_coeff_flag;
        if (significant_coeff_flag) {
            bool last_significant_coeff_flag = cabac.decode_decision(last_ctx + pos2ctx_Last[ii]);
            if (last_significant_coeff_flag)
                numCoeff = ii + 1;
        }
        ii++;
    }

    cabac_context_t* one_ctx = slice.parser.mot_ctx.one_contexts + type2ctx_one[context];

    coeff = coeff_val + numCoeff - 1;
    *coeff = 1;

    int numDecodAbsLevelEq1 = 0;
    int numDecodAbsLevelGt1 = 0;

    int i = chroma ? blkIdx % 2 : ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int j = chroma ? blkIdx / 2 : ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    for (ii = numCoeff - 1 + startIdx; ii >= startIdx; ii--) {
        if (*coeff) {
            int ctxIdxIncs[] = {
                (numDecodAbsLevelGt1 != 0 ? 0 : min(4, 1 + numDecodAbsLevelEq1)),
                5 + min(4 - (ctxBlockCat == CHROMA_DC ? 1 : 0), numDecodAbsLevelGt1)
            };
            int32_t coeff_abs_level_minus1 = unary_exp_golomb_level_decode(&cabac, one_ctx, ctxIdxIncs);
            bool    coeff_sign_flag = cabac.decode_bypass();
            int32_t coeffLevel = (coeff_abs_level_minus1 + 1) * (1 - 2 * coeff_sign_flag);

            *coeff = coeffLevel;

            numDecodAbsLevelEq1 += (coeff_abs_level_minus1 == 0);
            numDecodAbsLevelGt1 += (coeff_abs_level_minus1 != 0);

            //if (!ac)
            //    assert(startIdx + ii < numCoeff);
            if (!chroma) {
                if (!ac)
                    slice.decoder.coeff_luma_dc(&mb, pl, i, j, ii, *coeff);
                else
                    slice.decoder.coeff_luma_ac(&mb, pl, i, j, ii, *coeff);
            } else {
                if (!ac)
                    slice.decoder.coeff_chroma_dc(&mb, pl, i, j, ii, *coeff);
                else
                    slice.decoder.coeff_chroma_ac(&mb, pl, i, j, ii, *coeff);
            }
        }
        coeff--;
    }
}

void Parser::Residual::residual_luma(ColorPlane pl)
{
    auto residual_block = !pps.entropy_coding_mode_flag ?
        std::mem_fn(&Parser::Residual::residual_block_cavlc) :
        std::mem_fn(&Parser::Residual::residual_block_cabac);

    if (mb.mb_type == I_16x16 && !mb.dpl_flag) {
        residual_block(this, LUMA_16DC, 0, 15, 16, pl, false, false, 0);

        slice.decoder.transform_luma_dc(&mb, pl);
    }

    for (int i8x8 = 0; i8x8 < 4; i8x8++) {
        if (!mb.transform_size_8x8_flag || !pps.entropy_coding_mode_flag) {
            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (mb.CodedBlockPatternLuma & (1 << i8x8)) {
                    if (mb.mb_type == I_16x16)
                        residual_block(this, LUMA_16AC, 1, 14, 15, pl, false, true, i8x8 * 4 + i4x4);
                    else
                        residual_block(this, LUMA_4x4, 0, 15, 16, pl, false, true, i8x8 * 4 + i4x4);
                } else {
                    if (!pps.entropy_coding_mode_flag) {
                        int i = (i8x8 % 2) * 2 + (i4x4 % 2);
                        int j = (i8x8 / 2) * 2 + (i4x4 / 2);
                        mb.nz_coeff[pl][j][i] = 0;
                    }
                }
            }
        } else if (mb.CodedBlockPatternLuma & (1 << i8x8))
            residual_block(this, LUMA_8x8, 0, 63, 64, pl, false, true, i8x8 * 4);
        else {
            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (!pps.entropy_coding_mode_flag) {
                    int i = (i8x8 % 2) * 2 + (i4x4 % 2);
                    int j = (i8x8 / 2) * 2 + (i4x4 / 2);
                    mb.nz_coeff[pl][j][i] = 0;
                }
            }
        }
    }
}

void Parser::Residual::residual_chroma()
{
    int NumC8x8 = 4 / (sps.SubWidthC * sps.SubHeightC);

    auto residual_block = !pps.entropy_coding_mode_flag ?
        std::mem_fn(&Parser::Residual::residual_block_cavlc) :
        std::mem_fn(&Parser::Residual::residual_block_cabac);

    if (mb.CodedBlockPatternChroma & 3) {
        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            residual_block(this, CHROMA_DC, 0, 4 * NumC8x8 - 1, 4 * NumC8x8, (ColorPlane)(iCbCr+1), true, false, 0);

            slice.decoder.transform_chroma_dc(&mb, (ColorPlane)(iCbCr + 1));
        }
    }

    for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
        for (int i8x8 = 0; i8x8 < NumC8x8; i8x8++) {
            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (mb.CodedBlockPatternChroma & 2)
                    residual_block(this, CHROMA_AC, 1, 14, 15, (ColorPlane)(iCbCr+1), true, true, i8x8 * 4 + i4x4);
                else {
                    if (!pps.entropy_coding_mode_flag) {
                        int i = (i4x4 % 2);
                        int j = (i4x4 / 2) + (i8x8 * 2);
                        mb.nz_coeff[iCbCr + 1][j][i] = 0;
                    }
                }
            }
        }
    }
}

void Parser::Residual::residual()
{
    this->residual_luma(PLANE_Y);
    if (sps.ChromaArrayType == 1 || sps.ChromaArrayType == 2)
        this->residual_chroma();
    else if (sps.ChromaArrayType == 3) {
        this->residual_luma(PLANE_U);
        this->residual_luma(PLANE_V);
    }
}


}
}
