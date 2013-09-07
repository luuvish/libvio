#include "global.h"
#include "slice.h"
#include "bitstream_cabac.h"
#include "data_partition.h"
#include "macroblock.h"
#include "neighbour.h"
#include "quantization.h"
#include "transform.h"

#include "mb_read_syntax.h"


using vio::h264::cabac_context_t;
using vio::h264::cabac_engine_t;


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


#define IS_I16MB(MB)    ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)

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

void macroblock_t::residual_block_cabac(int16_t coeffLevel[16], uint8_t startIdx, uint8_t endIdx,
                                        uint8_t maxNumCoeff, ColorPlane pl, int bx, int by)
{
}


void macroblock_t::residual_luma_cabac(ColorPlane pl)
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    int CodedBlockPatternLuma = this->cbp & 15;

    if (IS_I16MB(this) && !this->dpl_flag) {
        // residual_block(i16x16DClevel, 0, 15, 16);

        data_partition_t *dP = &slice->partArr[slice->dp_mode ? 1 : 0];
        syntax_element_t currSE;
        currSE.context = LUMA_16DC;

        int coef_ctr = -1;
        int level = 1;
        for (int k = 0; k < 17 && level != 0; ++k) {
            readRunLevel_CABAC(this, &currSE, &dP->de_cabac);
            level = currSE.value1;
            if (level != 0) {
                coef_ctr += currSE.value2 + 1;
                //coeffLevel[startIdx + coef_ctr] = level;
                quantization.coeff_luma_dc(this, (ColorPlane)0, 0, 0, coef_ctr, level);
            }
        }

        transform.inverse_luma_dc(this, pl);
    }

    data_partition_t *dP = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];
    syntax_element_t currSE;
    if (pl == PLANE_Y || sps->separate_colour_plane_flag)
        currSE.context = this->transform_size_8x8_flag ? LUMA_8x8 : IS_I16MB(this) ? LUMA_16AC : LUMA_4x4;
    else if (pl == PLANE_U)
        currSE.context = this->transform_size_8x8_flag ? CB_8x8 : IS_I16MB(this) ? CB_16AC : CB_4x4;
    else
        currSE.context = this->transform_size_8x8_flag ? CR_8x8 : IS_I16MB(this) ? CR_16AC : CR_4x4;

    int start_scan = IS_I16MB(this)? 1 : 0;
    int max4x4 = this->transform_size_8x8_flag ? 1 : 4;
    int maxK   = 16 * (5 - max4x4) + 1;

    for (int i8x8 = 0; i8x8 < 4; i8x8++) {
        if (CodedBlockPatternLuma & (1 << i8x8)) {
            for (int i4x4 = 0; i4x4 < max4x4; i4x4++) {
                //if (!this->transform_size_8x8_flag || !pps->entropy_coding_mode_flag)
                    //if (IS_I16MB(this))
                        // residual_block(i16x16AClevel[i8x8 * 4 + i4x4], max(0, startIdx - 1), endIdx - 1, 15);
                    //else
                        // residual_block(level4x4[i8x8 * 4 + i4x4], startIdx, endIdx, 16);
                //else
                    // residual_block(level8x8[i8x8], 4 * startIdx, 4 * endIdx + 3, 64);

                int i = (i8x8 % 2) * 2 + (i4x4 % 2);
                int j = (i8x8 / 2) * 2 + (i4x4 / 2);

                this->subblock_x = i * 4;
                this->subblock_y = j * 4;

                int coef_ctr = start_scan - 1;
                int level = 1;
                for (int k = start_scan; k < maxK && level != 0; ++k) {
                    readRunLevel_CABAC(this, &currSE, &dP->de_cabac);
                    level = currSE.value1;
                    if (level != 0) {
                        coef_ctr += currSE.value2 + 1;
                        quantization.coeff_luma_ac(this, pl, i, j, coef_ctr, level);
                    }
                }
            }
        }
    }
}

void macroblock_t::residual_chroma_cabac()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    int NumC8x8 = 4 / (sps->SubWidthC * sps->SubHeightC);

    int CodedBlockPatternChroma = this->cbp / 16;

    if (CodedBlockPatternChroma & 3) {      
        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            // residual_block(ChromaDCLevel[iCbCr], 0, 4 * NumC8x8 - 1, 4 * NumC8x8);

            data_partition_t *dP = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];
            syntax_element_t currSE;
            currSE.context = sps->ChromaArrayType == 1 ? CHROMA_DC : CHROMA_DC_2x4;
            this->is_v_block = iCbCr * 2;

            int coef_ctr = -1;
            int level = 1;
            for (int k = 0; k < NumC8x8 * 4 + 1 && level != 0; ++k) {
                readRunLevel_CABAC(this, &currSE, &dP->de_cabac);
                level = currSE.value1;
                if (level != 0) {
                    coef_ctr += currSE.value2 + 1;
                    assert(coef_ctr < NumC8x8 * 4);
                    quantization.coeff_chroma_dc(this, (ColorPlane)(iCbCr + 1), 0, 0, coef_ctr, level);
                }
            }

            transform.inverse_chroma_dc(this, (ColorPlane)(iCbCr + 1));
        }
    }

    if (CodedBlockPatternChroma & 2) {
        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            for (int i8x8 = 0; i8x8 < NumC8x8; i8x8++) {
                for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                    // residual_block(ChromaACLevel[iCbCr][i8x8 * 4 + i4x4], max(0, startIdx - 1), endIdx - 1, 15);

                    int i = (i4x4 % 2);
                    int j = (i4x4 / 2) + (i8x8 * 2);

                    data_partition_t *dP = &slice->partArr[slice->dp_mode ? (this->is_intra_block ? 1 : 2) : 0];
                    syntax_element_t currSE;
                    currSE.context = CHROMA_AC;
                    this->is_v_block = iCbCr;
                    this->subblock_x = i * 4;
                    this->subblock_y = j * 4;

                    int coef_ctr = 0;
                    int level = 1;
                    for (int k = 0; k < 16 && level != 0; ++k) {
                        readRunLevel_CABAC(this, &currSE, &dP->de_cabac);
                        level = currSE.value1;
                        if (level != 0) {
                            coef_ctr += currSE.value2 + 1;
                            quantization.coeff_chroma_ac(this, (ColorPlane)(iCbCr + 1), i, j, coef_ctr, level);
                        }
                    }
                }
            }
        }
    }
}
