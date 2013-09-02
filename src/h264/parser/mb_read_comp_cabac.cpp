#include "global.h"
#include "slice.h"
#include "bitstream_cabac.h"
#include "data_partition.h"
#include "macroblock.h"
#include "mb_read.h"
#include "neighbour.h"
#include "transform.h"

#include "mb_read_syntax.h"


// CABAC block types
typedef enum {
    LUMA_16DC     =   0,
    LUMA_16AC     =   1,
    LUMA_8x8      =   2,
    LUMA_8x4      =   3,
    LUMA_4x8      =   4,
    LUMA_4x4      =   5,
    CHROMA_DC     =   6,
    CHROMA_AC     =   7,
    CHROMA_DC_2x4 =   8,
    CHROMA_DC_4x4 =   9,
    CB_16DC       =  10,
    CB_16AC       =  11,
    CB_8x8        =  12,
    CB_8x4        =  13,
    CB_4x8        =  14,
    CB_4x4        =  15,
    CR_16DC       =  16,
    CR_16AC       =  17,
    CR_8x8        =  18,
    CR_8x4        =  19,
    CR_4x8        =  20,
    CR_4x4        =  21
} CABACBlockTypes;

struct syntax_element_t {
    int value1;  //!< numerical value of syntax element
    int value2;  //!< for blocked symbols, e.g. run/level
    int context; //!< CABAC context
};


#define IS_I16MB(MB)    ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)

static inline int rshift_rnd_sf(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}

static inline int get_bit(int64_t x, int n)
{
    return (int)(((x >> n) & 1));
}

static uint32_t unary_exp_golomb_level_decode(cabac_engine_t* dep_dp, cabac_context_t* ctx)
{
    const uint32_t exp_start = 13;
    uint32_t symbol = 0;
    uint32_t binary_symbol = 0;
    uint32_t l, k = 1;

    if (dep_dp->decode_decision(ctx) == 0)
        return 0;

    do {
        l = dep_dp->decode_decision(ctx);
        ++symbol;
        ++k;
    } while (l != 0 && k != exp_start);

    if (l == 0)
        return symbol;

    k = 0;
    do {
        l = dep_dp->decode_bypass();
        if (l == 1) {
            symbol += (1 << k);
            ++k;
        }
    } while (l != 0);

    while (k--) {
        if (dep_dp->decode_bypass() == 1)
            binary_symbol |= (1 << k);
    }

    return symbol + binary_symbol + 1;
}

static int read_and_store_CBP_block_bit(mb_t *currMB, cabac_engine_t *dep_dp, int type)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    cabac_contexts_t *mot_ctx = currSlice->mot_ctx;
    mb_t *mb_data = currSlice->mb_data;
    int y_ac        = (type==LUMA_16AC || type==LUMA_8x8 || type==LUMA_8x4 || type==LUMA_4x8 || type==LUMA_4x4
                      || type==CB_16AC || type==CB_8x8 || type==CB_8x4 || type==CB_4x8 || type==CB_4x4
                      || type==CR_16AC || type==CR_8x8 || type==CR_8x4 || type==CR_4x8 || type==CR_4x4);
    int y_dc        = (type==LUMA_16DC || type==CB_16DC || type==CR_16DC); 
    int u_ac        = (type==CHROMA_AC && !currMB->is_v_block);
    int v_ac        = (type==CHROMA_AC &&  currMB->is_v_block);
    int chroma_dc   = (type==CHROMA_DC || type==CHROMA_DC_2x4 || type==CHROMA_DC_4x4);
    int u_dc        = (chroma_dc && !currMB->is_v_block);
    int v_dc        = (chroma_dc &&  currMB->is_v_block);
    int i           = (y_ac || u_ac || v_ac ? currMB->subblock_x : 0);
    int j           = (y_ac || u_ac || v_ac ? currMB->subblock_y : 0);
    int bit         = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35);
    int default_bit = (currMB->is_intra_block ? 1 : 0);
    int upper_bit   = default_bit;
    int left_bit    = default_bit;
    int cbp_bit     = 1;  // always one for 8x8 mode
    int ctx;
    int bit_pos_a   = 0;
    int bit_pos_b   = 0;

    int size_8x8_flag = (type == LUMA_8x8 || type == CB_8x8 || type == CR_8x8);
    int pl = (type == CB_8x8 || type == CB_4x4 || type == CB_4x8 || type == CB_8x4 || type == CB_16AC || type == CB_16DC) ? 1 :
             (type == CR_8x8 || type == CR_4x4 || type == CR_4x8 || type == CR_8x4 || type == CR_16AC || type == CR_16DC) ? 2 : 0;

    int mb_size[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };

    PixelPos block_a, block_b;
    get4x4Neighbour(currMB, i - 1, j, mb_size[y_dc || y_ac ? IS_LUMA : IS_CHROMA], &block_a);
    get4x4Neighbour(currMB, i, j - 1, mb_size[y_dc || y_ac ? IS_LUMA : IS_CHROMA], &block_b);
    if (y_ac || u_ac || v_ac) {
        if (block_a.available)
            bit_pos_a = 4 * block_a.y + block_a.x;
        if (block_b.available)
            bit_pos_b = 4 * block_b.y + block_b.x;
    }

    if ((sps->separate_colour_plane_flag || sps->chroma_format_idc != YUV444) && type != LUMA_8x8) {
        if (block_b.available) {
            if (mb_data[block_b.mb_addr].mb_type == IPCM)
                upper_bit = 1;
            else
                upper_bit = get_bit(mb_data[block_b.mb_addr].s_cbp[0].bits, bit + bit_pos_b);
        }
        if (block_a.available) {
            if (mb_data[block_a.mb_addr].mb_type == IPCM)
                left_bit = 1;
            else
                left_bit = get_bit(mb_data[block_a.mb_addr].s_cbp[0].bits, bit + bit_pos_a);
        }
    } else if (sps->chroma_format_idc == YUV444) {
        if (block_b.available) {
            if (!size_8x8_flag || mb_data[block_b.mb_addr].transform_size_8x8_flag) {
                if (mb_data[block_b.mb_addr].mb_type == IPCM)
                    upper_bit = 1;
                else
                    upper_bit = get_bit(mb_data[block_b.mb_addr].s_cbp[pl].bits, bit + bit_pos_b);
            }
        }
        if (block_a.available) {
            if (!size_8x8_flag || mb_data[block_a.mb_addr].transform_size_8x8_flag) {
                if (mb_data[block_a.mb_addr].mb_type == IPCM)
                    left_bit = 1;
                else
                    left_bit = get_bit(mb_data[block_a.mb_addr].s_cbp[pl].bits, bit + bit_pos_a);
            }
        }
    }

    static const short type2ctx_bcbp[] = {
         0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 10,
        11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20
    };

    if (sps->chroma_format_idc == YUV444 || type != LUMA_8x8) {
        ctx = 2 * upper_bit + left_bit;
        cbp_bit = dep_dp->decode_decision(mot_ctx->bcbp_contexts[type2ctx_bcbp[type]] + ctx);
    }

    bit = y_dc ?  0 :
          y_ac ?  1 + j + (i >> 2) :
          u_dc ? 17 :
          v_dc ? 18 :
          u_ac ? 19 + j + (i >> 2) :
                 35 + j + (i >> 2);

    if (cbp_bit) {
        if (sps->chroma_format_idc == YUV444) {
            if (size_8x8_flag)
                currMB->s_cbp[pl].bits |= ((int64_t)0x33 << bit);
            else if (type == LUMA_8x4 || type == CB_8x4 || type == CR_8x4)
                currMB->s_cbp[pl].bits |= ((int64_t)0x03 << bit);
            else if (type == LUMA_4x8 || type == CB_4x8 || type == CR_4x8)
                currMB->s_cbp[pl].bits |= ((int64_t)0x11 << bit);
            else
                currMB->s_cbp[pl].bits |= ((int64_t)0x01 << bit);
        } else {
            if (type == LUMA_8x8)
                currMB->s_cbp[0].bits |= ((int64_t)0x33 << bit);
            else if (type == LUMA_8x4)
                currMB->s_cbp[0].bits |= ((int64_t)0x03 << bit);
            else if (type == LUMA_4x8)
                currMB->s_cbp[0].bits |= ((int64_t)0x11 << bit);
            else
                currMB->s_cbp[0].bits |= ((int64_t)0x01 << bit);
        }
    }
    return cbp_bit;
}

static void readRunLevel_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    static const short maxpos[] = {
        15, 14, 63, 31, 31, 15,  3, 14,  7, 15, 15,
        14, 63, 31, 31, 15, 15, 14, 63, 31, 31, 15
    };
    static const short c1isdc[] = {
         1,  0,  1,  1,  1,  1,  1,  0,  1,  1,  1,
         0,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1
    };
    static const short max_c2[] = {
         4,  4,  4,  4,  4,  4,  3,  4,  3,  3,  4,
         4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
    }; // 9

    static const short type2ctx_map[] = {
         0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 10,
        11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
    }; // 8
    static const short type2ctx_one[] = {
         0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 10,
        11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20
    }; // 7

    //===== position -> ctx for MAP =====
    //--- zig-zag scan ----
    static const uint8_t pos2ctx_map8x8[] = {
         0,  1,  2,  3,  4,  5,  5,  4,  4,  3,  3,  4,  4,  4,  5,  5,
         4,  4,  4,  4,  3,  3,  6,  7,  7,  7,  8,  9, 10,  9,  8,  7,
         7,  6, 11, 12, 13, 11,  6,  7,  8,  9, 14, 10,  9,  8,  6, 11,
        12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map8x4[] = {
        0,  1,  2,  3,  4,  5,  7,  8,  9, 10, 11,  9,  8,  6,  7,  8,
        9, 10, 11,  9,  8,  6, 12,  8,  9, 10, 11,  9, 13, 13, 14, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map4x4[] = {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map2x4c[] = {
        0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2
    }; // 15 CTX
    static const uint8_t pos2ctx_map4x4c[] = {
        0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2
    }; // 15 CTX

    //--- interlace scan ----
    //taken from ABT
    static const uint8_t pos2ctx_map8x8i[] = {
        0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  7,  7,  7,  8,  4,  5,
        6,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 11, 12, 11,
        9,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 13, 13,  9,
        9, 10, 10,  8, 13, 13,  9,  9, 10, 10, 14, 14, 14, 14, 14, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map8x4i[] = {
        0,  1,  2,  3,  4,  5,  6,  3,  4,  5,  6,  3,  4,  7,  6,  8,
        9,  7,  6,  8,  9, 10, 11, 12, 12, 10, 11, 13, 13, 14, 14, 14
    }; // 15 CTX
    static const uint8_t pos2ctx_map4x8i[] = {
        0,  1,  1,  1,  2,  3,  3,  4,  4,  4,  5,  6,  2,  7,  7,  8,
        8,  8,  5,  6,  9, 10, 10, 11, 11, 11, 12, 13, 13, 14, 14, 14
    }; // 15 CTX

    static const uint8_t *pos2ctx_map[2][22] = {
        { pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map8x8,  pos2ctx_map8x4,
          pos2ctx_map8x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map2x4c, pos2ctx_map4x4c, 
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map8x8,  pos2ctx_map8x4,
          pos2ctx_map8x4,  pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map8x8,  pos2ctx_map8x4,
          pos2ctx_map8x4,  pos2ctx_map4x4 },
        { pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map8x8i, pos2ctx_map8x4i,
          pos2ctx_map4x8i, pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map2x4c, pos2ctx_map4x4c,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map8x8i, pos2ctx_map8x4i,
          pos2ctx_map8x4i, pos2ctx_map4x4,
          pos2ctx_map4x4,  pos2ctx_map4x4,
          pos2ctx_map8x8i, pos2ctx_map8x4i,
          pos2ctx_map8x4i, pos2ctx_map4x4 }
    };

    //===== position -> ctx for LAST =====
    static const uint8_t pos2ctx_last8x8[] = {
        0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
        2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
        3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
        5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8
    }; //  9 CTX
    static const uint8_t pos2ctx_last8x4[] = {
        0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,
        3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8
    }; //  9 CTX
    static const uint8_t pos2ctx_last4x4[] = {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15
    }; // 15 CTX
    static const uint8_t pos2ctx_last2x4c[] = {
        0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2
    }; // 15 CTX
    static const uint8_t pos2ctx_last4x4c[] = {
        0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2
    }; // 15 CTX

    static const uint8_t *pos2ctx_last[22] = {
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last8x8,  pos2ctx_last8x4,
        pos2ctx_last8x4,  pos2ctx_last4x4,
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last2x4c, pos2ctx_last4x4c,
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last8x8,  pos2ctx_last8x4,
        pos2ctx_last8x4,  pos2ctx_last4x4,
        pos2ctx_last4x4,  pos2ctx_last4x4,
        pos2ctx_last8x8,  pos2ctx_last8x4,
        pos2ctx_last8x4,  pos2ctx_last4x4
    };

    slice_t *currSlice = currMB->p_Slice;

    bool field = (currSlice->field_pic_flag || currMB->mb_field_decoding_flag);
    const byte *pos2ctx_Map  = pos2ctx_map [field][se->context];
    const byte *pos2ctx_Last = pos2ctx_last[se->context];
    cabac_context_t *map_ctx  = currSlice->mot_ctx->map_contexts [field][type2ctx_map[se->context]];
    cabac_context_t *last_ctx = currSlice->mot_ctx->last_contexts[field][type2ctx_map[se->context]];
    cabac_context_t *one_ctx = currSlice->mot_ctx->one_contexts[type2ctx_one[se->context]];
    cabac_context_t *abs_ctx = currSlice->mot_ctx->abs_contexts[type2ctx_one[se->context]];
    const short max_type = max_c2[se->context];

    if (currSlice->coeff_ctr < 0) {
        currSlice->coeff_ctr = read_and_store_CBP_block_bit(currMB, dep_dp, se->context);

        if (currSlice->coeff_ctr) {
            int i;
            int *coeff = currSlice->coeff;
            int i0     = 0;
            int i1     = maxpos[se->context];
            currSlice->coeff_ctr = 0;

            if (!c1isdc[se->context]) {
                ++i0;
                ++i1;
            }

            for (i = i0; i < i1; ++i) {
                int c = dep_dp->decode_decision(map_ctx + pos2ctx_Map[i]);
                *(coeff++) = c;
                currSlice->coeff_ctr += c;
                if (c && dep_dp->decode_decision(last_ctx + pos2ctx_Last[i])) {
                    memset(coeff, 0, (i1 - i) * sizeof(int));
                    i = i1 + 1;
                    break;
                }
            }
            if (i <= i1) {
                int c = 1;
                *(coeff++) = c;
                currSlice->coeff_ctr += c;
            }

            i = maxpos[se->context];
            int *cof = currSlice->coeff + i;
            int c1 = 1;
            int c2 = 0;

            for (; i >= 0; i--) {
                if (*cof != 0) {
                    *cof += dep_dp->decode_decision(one_ctx + c1);

                    if (*cof == 2) {
                        *cof += unary_exp_golomb_level_decode(dep_dp, abs_ctx + c2);
                        c2 = imin(++c2, max_type);
                        c1 = 0;
                    } else if (c1)
                        c1 = imin(++c1, 4);

                    if (dep_dp->decode_bypass())
                        *cof = - *cof;
                }
                cof--;
            }
        }
    }

    se->value1 = 0;
    se->value2 = 0;
    if (currSlice->coeff_ctr) {
        while (currSlice->coeff[currSlice->pos] == 0) {
            ++currSlice->pos;
            ++se->value2;
        }
        se->value1 = currSlice->coeff[currSlice->pos++];
    }
    if ((currSlice->coeff_ctr)-- == 0) 
        currSlice->pos = 0;
}


static void read_tc_luma(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    const byte (*pos_scan8x8)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN8x8 : FIELD_SCAN8x8;

    if (IS_I16MB(currMB) && !currMB->dpl_flag) {
        data_partition_t *dP = &currSlice->partArr[currSlice->dp_mode ? 1 : 0];
        syntax_element_t currSE;
        currSE.context = LUMA_16DC;

        int coef_ctr = -1;
        int level = 1;
        for (int k = 0; k < 17 && level != 0; ++k) {
            readRunLevel_CABAC(currMB, &currSE, &dP->de_cabac);
            level = currSE.value1;
            if (level != 0) {
                coef_ctr += currSE.value2 + 1;
                //coeffLevel[startIdx + coef_ctr] = level;
                int i0 = pos_scan4x4[coef_ctr][0];
                int j0 = pos_scan4x4[coef_ctr][1];
                currSlice->cof[0][j0 * 4][i0 * 4] = level;
            }
        }

        if (!currMB->TransformBypassModeFlag)
            itrans_2(currMB, pl);
    }

    if (!currMB->transform_size_8x8_flag) {
        int qp_per = currMB->qp_scaled[pl] / 6;
        int qp_rem = currMB->qp_scaled[pl] % 6;
        int transform_pl = sps->separate_colour_plane_flag ? currSlice->colour_plane_id : pl;
        int (*InvLevelScale4x4)[4] = currMB->is_intra_block ?
            currSlice->InvLevelScale4x4_Intra[transform_pl][qp_rem] :
            currSlice->InvLevelScale4x4_Inter[transform_pl][qp_rem];

        syntax_element_t currSE;
        if (pl == PLANE_Y || sps->separate_colour_plane_flag)
            currSE.context = IS_I16MB(currMB) ? LUMA_16AC : LUMA_4x4;
        else if (pl == PLANE_U)
            currSE.context = IS_I16MB(currMB) ? CB_16AC : CB_4x4;
        else
            currSE.context = IS_I16MB(currMB) ? CR_16AC : CR_4x4;

        int start_scan = IS_I16MB (currMB)? 1 : 0;

        for (int i8x8 = 0; i8x8 < 4; i8x8++) {
            int block_x = (i8x8 % 2) * 2;
            int block_y = (i8x8 / 2) * 2;

            for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                if (currMB->cbp & (1 << i8x8)) {
                    int i = block_x + (i4x4 % 2);
                    int j = block_y + (i4x4 / 2);

                    currMB->subblock_x = i * 4;
                    currMB->subblock_y = j * 4;

                    int coef_ctr = start_scan - 1;
                    int level = 1;
                    for (int k = start_scan; k < 17 && level != 0; ++k) {
                        data_partition_t *dP = &currSlice->partArr[currSlice->dp_mode ? (currMB->is_intra_block ? 1 : 2) : 0];

                        readRunLevel_CABAC(currMB, &currSE, &dP->de_cabac);
                        level = currSE.value1;
                        if (level != 0) {
                            coef_ctr += currSE.value2 + 1;
                            currMB->s_cbp[pl].blk |= (int64_t)(0x01 << (j * 4 + i));
                            int i0 = pos_scan4x4[coef_ctr][0];
                            int j0 = pos_scan4x4[coef_ctr][1];

                            if (!currMB->TransformBypassModeFlag)
                                currSlice->cof[pl][j * 4 + j0][i * 4 + i0] = rshift_rnd_sf((level * InvLevelScale4x4[j0][i0]) << qp_per, 4);
                            else
                                currSlice->cof[pl][j * 4 + j0][i * 4 + i0] = level;
                        }
                    }
                }
            }
        }
    } else {
        int qp_per = currMB->qp_scaled[pl] / 6;
        int qp_rem = currMB->qp_scaled[pl] % 6;
        int transform_pl = sps->separate_colour_plane_flag ? currSlice->colour_plane_id : pl;
        int (*InvLevelScale8x8)[8] = currMB->is_intra_block ?
            currSlice->InvLevelScale8x8_Intra[transform_pl][qp_rem] :
            currSlice->InvLevelScale8x8_Inter[transform_pl][qp_rem];

        syntax_element_t currSE;
        if (pl == PLANE_Y || sps->separate_colour_plane_flag)
            currSE.context = LUMA_8x8;
        else if (pl == PLANE_U)
            currSE.context = CB_8x8;
        else
            currSE.context = CR_8x8;  

        for (int i8x8 = 0; i8x8 < 4; i8x8++) {
            int block_x = (i8x8 % 2) * 2;
            int block_y = (i8x8 / 2) * 2;

            if (currMB->cbp & (1 << i8x8)) {
                currMB->subblock_x = block_x * 4;
                currMB->subblock_y = block_y * 4;

                int coef_ctr = -1;
                int level = 1;
                for (int k = 0; k < 65 && level != 0; ++k) {
                    data_partition_t *dP = &currSlice->partArr[currSlice->dp_mode ? (currMB->is_intra_block ? 1 : 2) : 0];

                    readRunLevel_CABAC(currMB, &currSE, &dP->de_cabac);
                    level = currSE.value1;
                    if (level != 0) {
                        coef_ctr += currSE.value2 + 1;
                        currMB->s_cbp[pl].blk |= (int64_t)(0x33 << (block_y * 4 + block_x));
                        int i0 = pos_scan8x8[coef_ctr][0];
                        int j0 = pos_scan8x8[coef_ctr][1];

                        if (!currMB->TransformBypassModeFlag)
                            currSlice->cof[pl][block_y * 4 + j0][block_x * 4 + i0] = rshift_rnd_sf((level * InvLevelScale8x8[j0][i0]) << qp_per, 6);
                        else
                            currSlice->cof[pl][block_y * 4 + j0][block_x * 4 + i0] = level;
                    }
                }
            }
        }
    }
}

static void read_tc_chroma(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int NumC8x8 = 4 / (sps->SubWidthC * sps->SubHeightC);

    if (currMB->cbp > 15) {      
        data_partition_t *dP = &currSlice->partArr[currSlice->dp_mode ? (currMB->is_intra_block ? 1 : 2) : 0];
        syntax_element_t currSE;
        currSE.context = sps->ChromaArrayType == 1 ? CHROMA_DC : CHROMA_DC_2x4;

        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            currMB->is_v_block = iCbCr * 2;

            int coef_ctr = -1;
            int level = 1;
            for (int k = 0; k < NumC8x8 * 4 + 1 && level != 0; ++k) {
                readRunLevel_CABAC(currMB, &currSE, &dP->de_cabac);
                level = currSE.value1;

                int i0, j0;
                if (level != 0) {
                    coef_ctr += currSE.value2 + 1;
                    assert(coef_ctr < NumC8x8 * 4);
                    if (sps->ChromaArrayType == 1) {
                        i0 = coef_ctr % 2;
                        j0 = coef_ctr / 2;
                        currMB->s_cbp[0].blk |= (int64_t)(0xf << (iCbCr * 4 + 16));
                    }
                    if (sps->ChromaArrayType == 2) {
                        i0 = FIELD_SCAN[coef_ctr][0];
                        j0 = FIELD_SCAN[coef_ctr][1];
                        currMB->s_cbp[0].blk |= (int64_t)(0xff << (iCbCr * 8 + 16));
                    }
                    if (coef_ctr < NumC8x8 * 4)
                        currSlice->cof[iCbCr + 1][j0 * 4][i0 * 4] = level;
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
        data_partition_t *dP = &currSlice->partArr[currSlice->dp_mode ? (currMB->is_intra_block ? 1 : 2) : 0];
        syntax_element_t currSE;
        currSE.context = CHROMA_AC;

        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
            int qp_per_uv = currMB->qp_scaled[iCbCr + 1] / 6;
            int qp_rem_uv = currMB->qp_scaled[iCbCr + 1] % 6;

            for (int i8x8 = 0; i8x8 < NumC8x8; i8x8++) {
                currMB->is_v_block = iCbCr;

                int (*InvLevelScale4x4)[4] = NULL;
                if (!currMB->TransformBypassModeFlag)
                    InvLevelScale4x4 = currMB->is_intra_block ?
                        currSlice->InvLevelScale4x4_Intra[iCbCr + 1][qp_rem_uv] :
                        currSlice->InvLevelScale4x4_Inter[iCbCr + 1][qp_rem_uv];

                for (int i4x4 = 0; i4x4 < 4; i4x4++) {
                    int i = (i4x4 % 2);
                    int j = (i4x4 / 2) + (i8x8 * 2);

                    currMB->subblock_x = i * 4;
                    currMB->subblock_y = j * 4;

                    int coef_ctr = 0;
                    int level = 1;
                    for (int k = 0; k < 16 && level != 0; ++k) {
                        readRunLevel_CABAC(currMB, &currSE, &dP->de_cabac);
                        level = currSE.value1;

                        if (level != 0) {
                            currMB->s_cbp[0].blk |= (int64_t)(0x1 << (i8x8 * 4 + i4x4 + 16));
                            coef_ctr += currSE.value2 + 1;
                            int i0 = pos_scan4x4[coef_ctr][0];
                            int j0 = pos_scan4x4[coef_ctr][1];

                            if (!currMB->TransformBypassModeFlag)
                                currSlice->cof[iCbCr + 1][j * 4 + j0][i * 4 + i0] = rshift_rnd_sf((level * InvLevelScale4x4[j0][i0]) << qp_per_uv, 4);
                            else
                                currSlice->cof[iCbCr + 1][j * 4 + j0][i * 4 + i0] = level;
                        }
                    }
                }
            }
        }
    }
}


void macroblock_t::read_CBP_and_coeffs_from_NAL_CABAC()
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

void macroblock_t::residual_block_cabac(int16_t coeffLevel[16], uint8_t startIdx, uint8_t endIdx,
                                        uint8_t maxNumCoeff, ColorPlane pl, int bx, int by)
{
    slice_t *slice = this->p_Slice;

    data_partition_t *dP = &slice->partArr[slice->dp_mode ? 1 : 0];
    syntax_element_t currSE;
    currSE.context = LUMA_16DC;

    const byte (*pos_scan4x4)[2] = !slice->field_pic_flag && !this->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    int coef_ctr = -1;
    int level = 1;
    for (int k = 0; k < 17 && level != 0; ++k) {
        readRunLevel_CABAC(this, &currSE, &dP->de_cabac);
        level = currSE.value1;
        if (level != 0) {
            coef_ctr += currSE.value2 + 1;
            coeffLevel[startIdx + coef_ctr] = level;
            int i0 = pos_scan4x4[coef_ctr][0];
            int j0 = pos_scan4x4[coef_ctr][1];
            slice->cof[0][j0 << 2][i0 << 2] = level;
        }
    }

    if (!this->TransformBypassModeFlag)
        itrans_2(this, pl);
}
