
#include <math.h>

#include "global.h"
#include "slice.h"
#include "dpb.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "macroblock.h"
#include "mb_read.h"
#include "fmo.h"
#include "image.h"
#include "neighbour.h"
#include "transform.h"
#include "mv_prediction.h"
#include "intra_prediction.h"
#include "inter_prediction.h"

#include "mb_read_syntax.h"


static void readMB_typeInfo_CABAC_i_slice(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    uint8_t mb_type = 0;

    if (currSlice->slice_type == SI_slice) {
        int condTermFlagA = currMB->mb_left && currMB->mb_left->mb_type != SI4MB ? 1 : 0;
        int condTermFlagB = currMB->mb_up   && currMB->mb_up  ->mb_type != SI4MB ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;
        se->context = ctxIdxInc; // store context

        mb_type = dep_dp->decode_decision(ctx->mb_type_contexts[1] + ctxIdxInc);
    }

    if (currSlice->slice_type == I_slice || mb_type) {
        int condTermFlagA = 0;
        int condTermFlagB = 0;
        if (currSlice->slice_type == I_slice) {
            condTermFlagA = currMB->mb_left && currMB->mb_left->mb_type != I4MB && currMB->mb_left->mb_type != I8MB ? 1 : 0;
            condTermFlagB = currMB->mb_up   && currMB->mb_up  ->mb_type != I4MB && currMB->mb_up->mb_type   != I8MB ? 1 : 0;
        } else if (currSlice->slice_type == SI_slice) {
            condTermFlagA = currMB->mb_left && currMB->mb_left->mb_type != I4MB ? 1 : 0;
            condTermFlagB = currMB->mb_up   && currMB->mb_up  ->mb_type != I4MB ? 1 : 0;
        }
        int ctxIdxInc = condTermFlagA + condTermFlagB;
        se->context = ctxIdxInc; // store context

        if (!dep_dp->decode_decision(ctx->mb_type_contexts[0] + ctxIdxInc))
            mb_type += 0;
        else if (!dep_dp->decode_terminate()) {
            mb_type += 1;
            mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[0] + 4) * 12;
            if (dep_dp->decode_decision(ctx->mb_type_contexts[0] + 5))
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[0] + 6) * 4 + 4;
            mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[0] + 7) * 2;
            mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[0] + 8);
        } else
            mb_type += 25;
    }

    se->value1 = mb_type;
}

static void readMB_typeInfo_CABAC_p_slice(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    uint8_t mb_type = 0;

    if (dep_dp->decode_decision(ctx->mb_type_contexts[1] + 4)) {
        mb_type = 6;
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 7);
    } else if (dep_dp->decode_decision(ctx->mb_type_contexts[1] + 5)) {
        mb_type = 3;
        mb_type -= dep_dp->decode_decision(ctx->mb_type_contexts[1] + 7);
    } else {
        mb_type = 1;
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 6) * 3;
    }

    if (mb_type <= 6)
        mb_type += 0;
    else if (!dep_dp->decode_terminate()) {
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 8) * 12;
        if (dep_dp->decode_decision(ctx->mb_type_contexts[1] + 9))
            mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 9) * 4 + 4;
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 10) * 2;
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 10);
    } else
        mb_type = 31;

    se->value1 = mb_type;
}

static void readMB_typeInfo_CABAC_b_slice(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    uint8_t mb_type = 0;

    int condTermFlagA = currMB->mb_left && currMB->mb_left->mb_type != 0 ? 1 : 0;
    int condTermFlagB = currMB->mb_up   && currMB->mb_up  ->mb_type != 0 ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    if (dep_dp->decode_decision(ctx->mb_type_contexts[2] + ctxIdxInc)) {
        if (dep_dp->decode_decision(ctx->mb_type_contexts[2] + 4)) {
            if (dep_dp->decode_decision(ctx->mb_type_contexts[2] + 5)) {
                mb_type = 12;
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6) * 8; 
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6) * 4;
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6) * 2;

                if (mb_type == 24)  
                    mb_type = 11;
                else if (mb_type == 26)  
                    mb_type = 22;
                else {
                    if (mb_type == 22)
                        mb_type = 23;
                    mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6);
                }
            } else {
                mb_type = 3;
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6) * 4;
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6) * 2;
                mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6);
            }
        } else
            mb_type = dep_dp->decode_decision(ctx->mb_type_contexts[2] + 6) + 1;
    } else
        mb_type = 0;

    if (mb_type <= 23)
        mb_type += 0;
    else if (!dep_dp->decode_terminate()) {
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 8) * 12;
        if (dep_dp->decode_decision(ctx->mb_type_contexts[1] + 9))
            mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 9) * 4 + 4;
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 10) * 2;
        mb_type += dep_dp->decode_decision(ctx->mb_type_contexts[1] + 10);
    } else
        mb_type = 48;

    se->value1 = mb_type;
}

static void readB8_typeInfo_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    uint8_t mb_type = 0;

    if (currSlice->slice_type == P_slice) {
        mb_type = 0;
        if (!dep_dp->decode_decision(ctx->b8_type_contexts[0] + 1)) {
            mb_type += 1;
            if (dep_dp->decode_decision(ctx->b8_type_contexts[0] + 3)) {
                mb_type += 1;
                mb_type += dep_dp->decode_decision(ctx->b8_type_contexts[0] + 4) ? 0 : 1;
            }
        }
    } else {
        mb_type = 0;
        if (dep_dp->decode_decision(ctx->b8_type_contexts[1])) {
            mb_type += 1;
            if (dep_dp->decode_decision(ctx->b8_type_contexts[1] + 1)) {
                mb_type += 2;
                if (dep_dp->decode_decision(ctx->b8_type_contexts[1] + 2)) {
                    mb_type += 4;
                    if (dep_dp->decode_decision(ctx->b8_type_contexts[1] + 3))
                        mb_type += 4;
                    else
                        mb_type += dep_dp->decode_decision(ctx->b8_type_contexts[1] + 3) * 2;
                } else
                    mb_type += dep_dp->decode_decision(ctx->b8_type_contexts[1] + 3) * 2;
            }
            mb_type += dep_dp->decode_decision(ctx->b8_type_contexts[1] + 3);
        }
    }

    se->value1 = mb_type;
}


static void readRefFrame_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    StorablePicture *dec_picture = currSlice->dec_picture;
    cabac_contexts_t *ctx = currSlice->mot_ctx;
    mb_t *neighborMB = NULL;

    int   addctx  = 0;
    int   a = 0, b = 0;
    int   act_ctx;
    int   act_sym;
    int   list = se->value2;

    PixelPos block_a, block_b;
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
    get4x4Neighbour(currMB, currMB->subblock_x - 1, currMB->subblock_y    , mb_size, &block_a);
    get4x4Neighbour(currMB, currMB->subblock_x,     currMB->subblock_y - 1, mb_size, &block_b);

#define IS_DIRECT(MB)   ((MB)->mb_type == 0 && (currSlice->slice_type == B_SLICE))

    if (block_b.available) {
        int b8b = ((block_b.x >> 1) & 0x01)+(block_b.y & 0x02);    
        neighborMB = &currSlice->mb_data[block_b.mb_addr];
        if (!( (neighborMB->mb_type==IPCM) || IS_DIRECT(neighborMB) || (neighborMB->b8mode[b8b]==0 && neighborMB->b8pdir[b8b]==2))) {
            if (currSlice->MbaffFrameFlag && (currMB->mb_field_decoding_flag == FALSE) && (neighborMB->mb_field_decoding_flag == TRUE))
                b = (dec_picture->mv_info[block_b.pos_y][block_b.pos_x].ref_idx[list] > 1 ? 2 : 0);
            else
                b = (dec_picture->mv_info[block_b.pos_y][block_b.pos_x].ref_idx[list] > 0 ? 2 : 0);
        }
    }

    if (block_a.available) {
        int b8a = ((block_a.x >> 1) & 0x01)+(block_a.y & 0x02);    
        neighborMB = &currSlice->mb_data[block_a.mb_addr];
        if (!((neighborMB->mb_type==IPCM) || IS_DIRECT(neighborMB) || (neighborMB->b8mode[b8a]==0 && neighborMB->b8pdir[b8a]==2))) {
            if (currSlice->MbaffFrameFlag && (currMB->mb_field_decoding_flag == FALSE) && (neighborMB->mb_field_decoding_flag == 1))
                a = (dec_picture->mv_info[block_a.pos_y][block_a.pos_x].ref_idx[list] > 1 ? 1 : 0);
            else
                a = (dec_picture->mv_info[block_a.pos_y][block_a.pos_x].ref_idx[list] > 0 ? 1 : 0);
        }
    }

#undef IS_DIRECT

    act_ctx = a + b;
    se->context = act_ctx; // store context

    //int binIdx[] = { act_ctx, 4, 5, 5, 5, 5, 5 };

    act_sym = dep_dp->decode_decision(ctx->ref_no_contexts[addctx] + act_ctx);
    if (act_sym != 0) {
        act_ctx = 4;
        act_sym = dep_dp->u(ctx->ref_no_contexts[addctx] + act_ctx, 1) + 1;
    }
    se->value1 = act_sym;
}

static unsigned int exp_golomb_decode_eq_prob(cabac_engine_t *dep_dp, int k)
{
    unsigned int l;
    int symbol = 0;
    int binary_symbol = 0;

    do {
        l = dep_dp->decode_bypass();
        if (l == 1) {
            symbol += (1<<k);
            ++k;
        }
    } while (l != 0);

    while (k--) //next binary part
        if (dep_dp->decode_bypass())
            binary_symbol |= (1<<k);

    return (unsigned int) (symbol + binary_symbol);
}

static unsigned int unary_exp_golomb_mv_decode(cabac_engine_t *dep_dp, cabac_context_t* ctx, unsigned int max_bin)
{
    unsigned int symbol = dep_dp->decode_decision(ctx);

    if (symbol == 0)
        return 0;
    else {
        unsigned int exp_start = 8;
        unsigned int l,k = 1;
        unsigned int bin = 1;

        symbol = 0;

        ++ctx;
        do {
            l = dep_dp->decode_decision(ctx);
            if ((++bin) == 2) ctx++;
            if (bin == max_bin) ++ctx;
            ++symbol;
            ++k;
        } while (l != 0 && k != exp_start);

        if (l != 0)
            symbol += exp_golomb_decode_eq_prob(dep_dp,3) + 1;
        return symbol;
    }
}

static void read_mvd_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;
    int i = currMB->subblock_x;
    int j = currMB->subblock_y;
    int a = 0, b = 0;
    int act_ctx;
    int act_sym;  
    int list_idx = se->value2 & 0x01;
    int k = (se->value2 >> 1); // MVD component

    PixelPos block_a, block_b;
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
    get4x4NeighbourBase(currMB, i - 1, j    , mb_size, &block_a);
    get4x4NeighbourBase(currMB, i    , j - 1, mb_size, &block_b);

    if (block_a.available) {
        a = iabs(currSlice->mb_data[block_a.mb_addr].mvd[list_idx][block_a.y][block_a.x][k]);
        if (currSlice->MbaffFrameFlag && (k==1)) {
            if ((currMB->mb_field_decoding_flag==0) && (currSlice->mb_data[block_a.mb_addr].mb_field_decoding_flag==1))
                a *= 2;
            else if ((currMB->mb_field_decoding_flag==1) && (currSlice->mb_data[block_a.mb_addr].mb_field_decoding_flag==0))
                a /= 2;
        }
    }
    if (block_b.available) {
        b = iabs(currSlice->mb_data[block_b.mb_addr].mvd[list_idx][block_b.y][block_b.x][k]);
        if (currSlice->MbaffFrameFlag && (k==1)) {
            if ((currMB->mb_field_decoding_flag==0) && (currSlice->mb_data[block_b.mb_addr].mb_field_decoding_flag==1))
                b *= 2;
            else if ((currMB->mb_field_decoding_flag==1) && (currSlice->mb_data[block_b.mb_addr].mb_field_decoding_flag==0))
                b /= 2;
        }
    }
    a += b;

    if (a < 3)
        act_ctx = 5 * k;
    else if (a > 32)
        act_ctx = 5 * k + 3;
    else
        act_ctx = 5 * k + 2;

    se->context = act_ctx;

    act_sym = dep_dp->decode_decision(&ctx->mv_res_contexts[0][act_ctx]);

    if (act_sym != 0) {
        act_ctx = 5 * k;
        act_sym = unary_exp_golomb_mv_decode(dep_dp, ctx->mv_res_contexts[1] + act_ctx, 3) + 1;

        if (dep_dp->decode_bypass())
            act_sym = -act_sym;
    }
    se->value1 = act_sym;
}

static void read_CBP_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int cbp = 0;

    //  coding of luma part (bit by bit)
    for (int mb_y = 0; mb_y < 4; mb_y += 2) {
        for (int mb_x = 0; mb_x < 4; mb_x += 2) {
            int condTermFlagA = 0;
            int condTermFlagB = 0;

            if (mb_x == 0) {
                PixelPos block_a;
                int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
                get4x4Neighbour(currMB, (mb_x << 2) - 1, (mb_y << 2), mb_size, &block_a);
                if (block_a.available && currSlice->mb_data[block_a.mb_addr].mb_type != IPCM) {
                    int a_cbp = currSlice->mb_data[block_a.mb_addr].cbp;
                    condTermFlagA = (a_cbp & (1 << (2 * (block_a.y >> 1) + 1))) == 0 ? 1 : 0;
                }
            } else
                condTermFlagA = (cbp & (1 << mb_y)) == 0 ? 1 : 0;
            if (mb_y == 0) {
                if (currMB->mb_up && currMB->mb_up->mb_type != IPCM) {
                    int b_cbp = currMB->mb_up->cbp;
                    condTermFlagB = (b_cbp & (1 << (2 + (mb_x >> 1)))) == 0 ? 2 : 0;
                }
            } else
                condTermFlagB = (cbp & (1 << (mb_x >> 1))) == 0 ? 2 : 0;

            int ctxIdxInc = condTermFlagA + condTermFlagB;
            if (dep_dp->decode_decision(ctx->cbp_contexts[0] + ctxIdxInc))
                cbp += (1 << (mb_y + (mb_x >> 1)));
        }
    }

    if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
        // coding of chroma part
        // CABAC decoding for BinIdx 0
        int condTermFlagA = currMB->mb_left && (currMB->mb_left->mb_type == IPCM || currMB->mb_left->cbp > 15) ? 1 : 0;
        int condTermFlagB = currMB->mb_up   && (currMB->mb_up  ->mb_type == IPCM || currMB->mb_up  ->cbp > 15) ? 2 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        // CABAC decoding for BinIdx 1
        if (dep_dp->decode_decision(ctx->cbp_contexts[1] + ctxIdxInc)) { // set the chroma bits
            condTermFlagA = currMB->mb_left && (currMB->mb_left->mb_type == IPCM || (currMB->mb_left->cbp >> 4) == 2) ? 1 : 0;
            condTermFlagB = currMB->mb_up   && (currMB->mb_up  ->mb_type == IPCM || (currMB->mb_up  ->cbp >> 4) == 2) ? 2 : 0;
            ctxIdxInc = condTermFlagA + condTermFlagB;

            cbp += dep_dp->decode_decision(ctx->cbp_contexts[2] + ctxIdxInc) ? 32 : 16;
        }
    }

    se->value1 = cbp;

    if (!cbp)
        currSlice->last_dquant = 0;
}



uint32_t parse_mb_skip_run(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;

    uint32_t mb_skip_run = 0;

    if (!pps->entropy_coding_mode_flag)
        mb_skip_run = bitstream->ue();

    return mb_skip_run;
}

bool parse_mb_skip_flag(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;
    cabac_contexts_t *ctx = slice->mot_ctx;

    bool mb_skip_flag = 0;

    if (pps->entropy_coding_mode_flag) {
        int tabIdx = slice->slice_type == B_slice ? 2 : 1;
        int ctxIdx = slice->slice_type == B_slice ? 7 : 0;

        int condTermFlagA = mb->mb_left && !mb->mb_left->mb_skip_flag ? 1 : 0;
        int condTermFlagB = mb->mb_up   && !mb->mb_up  ->mb_skip_flag ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        mb_skip_flag = dep_dp->decode_decision(&ctx->mb_type_contexts[tabIdx][ctxIdx + ctxIdxInc]);

        if (mb_skip_flag)
            slice->last_dquant = 0;
    }

    mb->ei_flag = 0;

    return mb_skip_flag;
}

bool parse_mb_field_decoding_flag(mb_t *mb)
{
	slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;
    cabac_contexts_t *ctx = slice->mot_ctx;

    bool mb_field_decoding_flag;

    if (!pps->entropy_coding_mode_flag)
        mb_field_decoding_flag = bitstream->f(1);
    else {
        int condTermFlagA = mb->mbAvailA && slice->mb_data[mb->mbAddrA].mb_field_decoding_flag ? 1 : 0;
        int condTermFlagB = mb->mbAvailB && slice->mb_data[mb->mbAddrB].mb_field_decoding_flag ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        mb_field_decoding_flag = dep_dp->decode_decision(ctx->mb_aff_contexts + ctxIdxInc);
    }

    return mb_field_decoding_flag;
}

uint32_t parse_mb_type(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;

    uint32_t mb_type;

    if (!pps->entropy_coding_mode_flag) {
        mb_type = bitstream->ue();
        if (slice->slice_type == P_slice || slice->slice_type == SP_slice)
            mb_type++;
    } else {
        syntax_element_t currSE;
        void (*reading)(macroblock_t*, syntax_element_t*, cabac_engine_t*) =
            slice->slice_type == I_slice || slice->slice_type == SI_slice ? readMB_typeInfo_CABAC_i_slice :
            slice->slice_type == P_slice || slice->slice_type == SP_slice ? readMB_typeInfo_CABAC_p_slice :
                                                                            readMB_typeInfo_CABAC_b_slice;
        reading(mb, &currSE, dep_dp);
        mb_type = currSE.value1;
    }

    mb->ei_flag = 0;

    return mb_type;
}

uint8_t parse_sub_mb_type(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;

    syntax_element_t currSE;
    if (!pps->entropy_coding_mode_flag) 
        currSE.value1 = bitstream->ue();
    else
        readB8_typeInfo_CABAC(mb, &currSE, dep_dp);

    return currSE.value1;
}


bool parse_transform_size_8x8_flag(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;
    cabac_contexts_t *ctx = slice->mot_ctx;

    bool transform_size_8x8_flag;

    if (!pps->entropy_coding_mode_flag)
        transform_size_8x8_flag = bitstream->f(1);
    else {
        int condTermFlagA = mb->mb_left && mb->mb_left->transform_size_8x8_flag ? 1 : 0;
        int condTermFlagB = mb->mb_up   && mb->mb_up  ->transform_size_8x8_flag ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        transform_size_8x8_flag = dep_dp->decode_decision(ctx->transform_size_contexts + ctxIdxInc);
    }

    return transform_size_8x8_flag;
}

int8_t parse_intra_pred_mode(mb_t *mb, uint8_t block4x4Idx)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;
    cabac_contexts_t *ctx = slice->mot_ctx;

    int8_t intra_pred_mode;

    if (!pps->entropy_coding_mode_flag) {
        if (bitstream->f(1))
            intra_pred_mode = -1;
        else
            intra_pred_mode = bitstream->f(3);
    } else {
        if (dep_dp->decode_decision(ctx->ipr_contexts))
            intra_pred_mode = -1;
        else {
            intra_pred_mode  = (dep_dp->decode_decision(ctx->ipr_contexts + 1)     );
            intra_pred_mode |= (dep_dp->decode_decision(ctx->ipr_contexts + 1) << 1);
            intra_pred_mode |= (dep_dp->decode_decision(ctx->ipr_contexts + 1) << 2);
        }
    }

    return intra_pred_mode;
}

uint8_t parse_intra_chroma_pred_mode(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;
    cabac_contexts_t *ctx = slice->mot_ctx;

    uint8_t intra_chroma_pred_mode;

    if (!pps->entropy_coding_mode_flag)
        intra_chroma_pred_mode = bitstream->ue();
    else {
        int condTermFlagA = mb->mb_left && mb->mb_left->intra_chroma_pred_mode != 0 && mb->mb_left->mb_type != IPCM ? 1 : 0;
        int condTermFlagB = mb->mb_up   && mb->mb_up  ->intra_chroma_pred_mode != 0 && mb->mb_up  ->mb_type != IPCM ? 1 : 0;
        int ctxIdxInc = condTermFlagA + condTermFlagB;

        // binIdx[] = { ctxIdxInc, 3, 3 };

        intra_chroma_pred_mode = dep_dp->decode_decision(ctx->cipr_contexts + ctxIdxInc);
        if (intra_chroma_pred_mode)
            intra_chroma_pred_mode = dep_dp->tu(ctx->cipr_contexts + 3, 0, 1) + 1;
    }

    return intra_chroma_pred_mode;
}

uint8_t parse_ref_idx(mb_t *mb, uint8_t b8mode, uint8_t list)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    bool refidx_present =
        slice->slice_type == B_slice || !slice->allrefzero || mb->mb_type != P8x8;
    int num_ref_idx_active = list == LIST_0 ?
        slice->num_ref_idx_l0_active_minus1 + 1 :
        slice->num_ref_idx_l1_active_minus1 + 1;

    if (!refidx_present || num_ref_idx_active <= 1)
        return 0;

    DataPartition *dP = &slice->partArr[0];
    syntax_element_t currSE;

    if (!pps->entropy_coding_mode_flag) {
        if (num_ref_idx_active == 2)
            currSE.value1 = 1 - dP->bitstream->f(1);
        else
            currSE.value1 = dP->bitstream->ue();
    } else {
        currSE.context = (b8mode >= 4);
        currSE.value2  = list;
        readRefFrame_CABAC(mb, &currSE, &dP->bitstream->de_cabac);
    }

    return currSE.value1;
}

int16_t parse_mvd(mb_t *mb, uint8_t xy, uint8_t list)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    DataPartition *dP = &slice->partArr[0];
    syntax_element_t currSE;

    if (!pps->entropy_coding_mode_flag) 
        currSE.value1 = dP->bitstream->se();
    else {
        currSE.value2 = xy * 2 + list;
        read_mvd_CABAC(mb, &currSE, &dP->bitstream->de_cabac);
    }

    return currSE.value1;
}

uint8_t parse_coded_block_pattern(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    sps_t *sps = slice->active_sps;
    pps_t *pps = slice->active_pps;

    DataPartition *dP = &slice->partArr[0];
    syntax_element_t currSE;

    if (!pps->entropy_coding_mode_flag) {
        //! gives CBP value from codeword number, both for intra and inter
        static const uint8_t NCBP[2][48][2] = {
            { { 15,  0 } , {  0,  1 } , {  7,  2 } , { 11,  4 } , { 13,  8 } , { 14,  3 },
              {  3,  5 } , {  5, 10 } , { 10, 12 } , { 12, 15 } , {  1,  7 } , {  2, 11 },
              {  4, 13 } , {  8, 14 } , {  6,  6 } , {  9,  9 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 },
              {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } , {  0,  0 } },
            { { 47,  0 } , { 31, 16 } , { 15,  1 } , {  0,  2 } , { 23,  4 } , { 27,  8 },
              { 29, 32 } , { 30,  3 } , {  7,  5 } , { 11, 10 } , { 13, 12 } , { 14, 15 },
              { 39, 47 } , { 43,  7 } , { 45, 11 } , { 46, 13 } , { 16, 14 } , {  3,  6 },
              {  5,  9 } , { 10, 31 } , { 12, 35 } , { 19, 37 } , { 21, 42 } , { 26, 44 },
              { 28, 33 } , { 35, 34 } , { 37, 36 } , { 42, 40 } , { 44, 39 } , {  1, 43 },
              {  2, 45 } , {  4, 46 } , {  8, 17 } , { 17, 18 } , { 18, 20 } , { 20, 24 },
              { 24, 19 } , {  6, 21 } , {  9, 26 } , { 22, 28 } , { 25, 23 } , { 32, 27 },
              { 33, 29 } , { 34, 30 } , { 36, 22 } , { 40, 25 } , { 38, 38 } , { 41, 41 } }
        };

        bool normal  = (sps->chroma_format_idc == 0 || sps->chroma_format_idc == 3 ? 0 : 1);
        bool inter   = (mb->is_intra_block ? 0 : 1);
        int  cbp_idx = dP->bitstream->ue();
        currSE.value1 = NCBP[normal][cbp_idx][inter];
    } else
        read_CBP_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    return currSE.value1;
}

int8_t parse_mb_qp_delta(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    Bitstream *bitstream = slice->partArr[0].bitstream;
    cabac_engine_t *dep_dp = &bitstream->de_cabac;
    cabac_contexts_t *ctx = slice->mot_ctx;

    int8_t mb_qp_delta;

    if (!pps->entropy_coding_mode_flag)
        mb_qp_delta = bitstream->se();
    else {
        int ctxIdxInc = slice->last_dquant != 0 ? 1 : 0;

        mb_qp_delta = dep_dp->decode_decision(ctx->delta_qp_contexts + ctxIdxInc);

        //int binIdx[] = { ctxIdxInc, 2, 3, 3, 3, 3, 3 };

        if (mb_qp_delta) {
            ctxIdxInc = 2;
            int act_sym = dep_dp->u(ctx->delta_qp_contexts + ctxIdxInc, 1) + 1;

            mb_qp_delta = (act_sym + 1) >> 1;
            if ((act_sym & 1) == 0) // lsb is signed bit
                mb_qp_delta = -mb_qp_delta;
        }

        slice->last_dquant = mb_qp_delta;
    }

    return mb_qp_delta;
}
