
#include <math.h>

#include "global.h"
#include "slice.h"
#include "dpb.h"
#include "bitstream_elements.h"
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



static void read_skip_flag_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;
    int tabIdx = currSlice->slice_type == B_SLICE ? 2 : 1;
    int ctxIdx = currSlice->slice_type == B_SLICE ? 7 : 0;

    int condTermFlagA = (currMB->mb_left != NULL) ? (currMB->mb_left->mb_skip_flag == 0) : 0;
    int condTermFlagB = (currMB->mb_up   != NULL) ? (currMB->mb_up  ->mb_skip_flag == 0) : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    se->value1 = dep_dp->decode_decision(&ctx->mb_type_contexts[tabIdx][ctxIdx + ctxIdxInc]) != 1;

    if (!se->value1)
        currMB->p_Slice->last_dquant = 0;
}

static void readFieldModeInfo_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int condTermFlagA = currMB->mbAvailA ? currSlice->mb_data[currMB->mbAddrA].mb_field_decoding_flag : 0;
    int condTermFlagB = currMB->mbAvailB ? currSlice->mb_data[currMB->mbAddrB].mb_field_decoding_flag : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    se->value1 = dep_dp->decode_decision(&ctx->mb_aff_contexts[ctxIdxInc]);
}

static void readMB_typeInfo_CABAC_i_slice(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int a = 0, b = 0;
    int act_ctx;
    int act_sym;
    int mode_sym;
    int curr_mb_type = 0;

    if (currSlice->slice_type == I_SLICE) { // INTRA-frame
        if (currMB->mb_up != NULL)
            b = (currMB->mb_up->mb_type != I4MB && currMB->mb_up->mb_type != I8MB ? 1 : 0 );
        if (currMB->mb_left != NULL)
            a = (currMB->mb_left->mb_type != I4MB && currMB->mb_left->mb_type != I8MB ? 1 : 0);

        act_ctx = a + b;
        act_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
        se->context = act_ctx; // store context

        if (act_sym == 0) // 4x4 Intra
            curr_mb_type = act_sym;
        else { // 16x16 Intra
            mode_sym = dep_dp->decode_terminate();
            if (mode_sym == 1)
                curr_mb_type = 25;
            else {
                act_sym = 1;
                act_ctx = 4;
                mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx); // decoding of AC/no AC
                act_sym += mode_sym * 12;
                act_ctx = 5;
                // decoding of cbp: 0,1,2
                mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                if (mode_sym != 0) {
                    act_ctx = 6;
                    mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += 4;
                    if (mode_sym != 0)
                        act_sym += 4;
                }
                // decoding of I pred-mode: 0,1,2,3
                act_ctx = 7;
                mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                act_sym += mode_sym * 2;
                act_ctx = 8;
                mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                act_sym += mode_sym;
                curr_mb_type = act_sym;
            }
        }
    } else if (currSlice->slice_type == SI_SLICE) { // SI-frame
        // special ctx's for SI4MB
        if (currMB->mb_up != NULL)
            b = (( (currMB->mb_up)->mb_type != SI4MB) ? 1 : 0 );
        if (currMB->mb_left != NULL)
            a = (( (currMB->mb_left)->mb_type != SI4MB) ? 1 : 0 );

        act_ctx = a + b;
        act_sym = dep_dp->decode_decision(ctx->mb_type_contexts[1] + act_ctx);
        se->context = act_ctx; // store context

        if (act_sym == 0) //  SI 4x4 Intra
            curr_mb_type = 0;
        else { // analog INTRA_IMG
            if (currMB->mb_up != NULL)
                b = (( (currMB->mb_up)->mb_type != I4MB) ? 1 : 0 );
            if (currMB->mb_left != NULL)
                a = (( (currMB->mb_left)->mb_type != I4MB) ? 1 : 0 );

            act_ctx = a + b;
            act_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
            se->context = act_ctx; // store context

            if (act_sym == 0) // 4x4 Intra
                curr_mb_type = 1;
            else { // 16x16 Intra
                mode_sym = dep_dp->decode_terminate();
                if (mode_sym == 1)
                    curr_mb_type = 26;
                else {
                    act_sym = 2;
                    act_ctx = 4;
                    mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx); // decoding of AC/no AC
                    act_sym += mode_sym * 12;
                    act_ctx = 5;
                    // decoding of cbp: 0,1,2
                    mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                    if (mode_sym != 0) {
                        act_ctx = 6;
                        mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                        act_sym += 4;
                        if (mode_sym != 0)
                          act_sym += 4;
                    }
                    // decoding of I pred-mode: 0,1,2,3
                    act_ctx = 7;
                    mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += mode_sym * 2;
                    act_ctx = 8;
                    mode_sym = dep_dp->decode_decision(ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += mode_sym;
                    curr_mb_type = act_sym;
                }
            }
        }
    }

    se->value1 = curr_mb_type;
}

static void readMB_typeInfo_CABAC_p_slice(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int act_ctx;
    int act_sym;
    int mode_sym;
    int curr_mb_type;
    cabac_context_t *mb_type_contexts = ctx->mb_type_contexts[1];

    if (dep_dp->decode_decision(&mb_type_contexts[4])) {
        if (dep_dp->decode_decision(&mb_type_contexts[7]))
            act_sym = 7;
        else
            act_sym = 6;
    } else {
        if (dep_dp->decode_decision(&mb_type_contexts[5])) {
            if (dep_dp->decode_decision(&mb_type_contexts[7])) 
                act_sym = 2;
            else
                act_sym = 3;
        } else {
            if (dep_dp->decode_decision(&mb_type_contexts[6]))
                act_sym = 4;
            else
                act_sym = 1;
        }
    }

    if (act_sym <= 6)
        curr_mb_type = act_sym;
    else { // additional info for 16x16 Intra-mode
        mode_sym = dep_dp->decode_terminate();
        if (mode_sym == 1)
            curr_mb_type = 31;
        else {
            act_ctx = 8;
            mode_sym =  dep_dp->decode_decision(mb_type_contexts + act_ctx); // decoding of AC/no AC
            act_sym += mode_sym*12;

            // decoding of cbp: 0,1,2
            act_ctx = 9;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
            if (mode_sym != 0) {
                act_sym += 4;
                mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
                if (mode_sym != 0)
                    act_sym += 4;
            }

            // decoding of I pred-mode: 0,1,2,3
            act_ctx = 10;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
            act_sym += mode_sym*2;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
            act_sym += mode_sym;
            curr_mb_type = act_sym;
        }
    }

    se->value1 = curr_mb_type;
}

static void readMB_typeInfo_CABAC_b_slice(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int a = 0, b = 0;
    int act_ctx;
    int act_sym;
    int mode_sym;
    int curr_mb_type;
    cabac_context_t *mb_type_contexts = ctx->mb_type_contexts[2];

    if (currMB->mb_up != NULL)
        b = (( (currMB->mb_up)->mb_type != 0) ? 1 : 0 );

    if (currMB->mb_left != NULL)
        a = (( (currMB->mb_left)->mb_type != 0) ? 1 : 0 );

    act_ctx = a + b;

    if (dep_dp->decode_decision(&mb_type_contexts[act_ctx])) {
        if (dep_dp->decode_decision(&mb_type_contexts[4])) {
            if (dep_dp->decode_decision(&mb_type_contexts[5])) {
                act_sym = 12;
                if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                    act_sym += 8;
                if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                    act_sym += 4;
                if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                    act_sym += 2;

                if (act_sym == 24)  
                    act_sym = 11;
                else if (act_sym == 26)  
                    act_sym = 22;
                else {
                    if (act_sym == 22)
                        act_sym = 23;
                    if (dep_dp->decode_decision(&mb_type_contexts[6]))
                        act_sym += 1;
                }
            } else {
                act_sym = 3;
                if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                    act_sym += 4;
                if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                    act_sym += 2;
                if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                    act_sym += 1;
            }
        } else {
            if (dep_dp->decode_decision(&mb_type_contexts[6])) 
                act_sym = 2;
            else
                act_sym = 1;
        }
    } else
        act_sym = 0;

    if (act_sym <= 23)
        curr_mb_type = act_sym;
    else { // additional info for 16x16 Intra-mode
        mode_sym = dep_dp->decode_terminate();
        if (mode_sym == 1)
            curr_mb_type = 48;
        else {
            mb_type_contexts = ctx->mb_type_contexts[1];
            act_ctx = 8;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx); // decoding of AC/no AC
            act_sym += mode_sym*12;

            // decoding of cbp: 0,1,2
            act_ctx = 9;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
            if (mode_sym != 0) {
                act_sym += 4;
                mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
                if (mode_sym != 0)
                    act_sym += 4;
            }

            // decoding of I pred-mode: 0,1,2,3
            act_ctx = 10;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
            act_sym += mode_sym * 2;
            mode_sym = dep_dp->decode_decision(mb_type_contexts + act_ctx);
            act_sym += mode_sym;
            curr_mb_type = act_sym;
        }
    }

    se->value1 = curr_mb_type;
}

static void readB8_typeInfo_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;
    cabac_context_t *b8_type_contexts;
    int act_sym = 0;

    if (currSlice->slice_type == P_slice) {
        b8_type_contexts = &ctx->b8_type_contexts[0][1];

        if (dep_dp->decode_decision(b8_type_contexts++))
            act_sym = 0;
        else {
            if (dep_dp->decode_decision(++b8_type_contexts))
                act_sym = (dep_dp->decode_decision(++b8_type_contexts)) ? 2: 3;
            else
                act_sym = 1;
        }
    } else {
        b8_type_contexts = ctx->b8_type_contexts[1];

        if (dep_dp->decode_decision(b8_type_contexts++)) {
            if (dep_dp->decode_decision(b8_type_contexts++)) {
                if (dep_dp->decode_decision(b8_type_contexts++)) {
                    if (dep_dp->decode_decision(b8_type_contexts)) {
                        act_sym = 10;
                        if (dep_dp->decode_decision(b8_type_contexts)) 
                            act_sym++;
                    } else {
                        act_sym = 6;
                        if (dep_dp->decode_decision(b8_type_contexts)) 
                            act_sym += 2;
                        if (dep_dp->decode_decision(b8_type_contexts)) 
                            act_sym++;
                    }
                } else {
                    act_sym = 2;
                    if (dep_dp->decode_decision(b8_type_contexts)) 
                        act_sym += 2;
                    if (dep_dp->decode_decision(b8_type_contexts)) 
                        act_sym++;
                }
            } else
                act_sym = dep_dp->decode_decision(++b8_type_contexts) ? 1 : 0;
            act_sym++;
        } else
            act_sym = 0;
    }

    se->value1 = act_sym;
}


static void readMB_transform_size_flag_CABAC( mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int b = (currMB->mb_up   == NULL) ? 0 : currMB->mb_up->transform_size_8x8_flag;
    int a = (currMB->mb_left == NULL) ? 0 : currMB->mb_left->transform_size_8x8_flag;

    int act_sym = dep_dp->decode_decision(ctx->transform_size_contexts + a + b);

    se->value1 = act_sym;
}

static void readIntraPredMode_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int act_sym = dep_dp->decode_decision(ctx->ipr_contexts);

    // remaining_mode_selector
    if (act_sym == 1)
        se->value1 = -1;
    else {
        se->value1  = (dep_dp->decode_decision(ctx->ipr_contexts + 1)     );
        se->value1 |= (dep_dp->decode_decision(ctx->ipr_contexts + 1) << 1);
        se->value1 |= (dep_dp->decode_decision(ctx->ipr_contexts + 1) << 2);
    }
}

static void readCIPredMode_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    mb_t *MbUp   = currMB->mb_up;
    mb_t *MbLeft = currMB->mb_left;

    int b = (MbUp != NULL)   ? (((MbUp->intra_chroma_pred_mode   != 0) && (MbUp->mb_type != IPCM)) ? 1 : 0) : 0;
    int a = (MbLeft != NULL) ? (((MbLeft->intra_chroma_pred_mode != 0) && (MbLeft->mb_type != IPCM)) ? 1 : 0) : 0;
    int act_ctx = a + b;
    int act_sym;

    // binIdx[] = { act_ctx, 3, 3 };

    act_sym = dep_dp->decode_decision(ctx->cipr_contexts + act_ctx);
    if (act_sym != 0)
        act_sym = dep_dp->tu(ctx->cipr_contexts + 3, 0, 1) + 1;

    se->value1 = act_sym;
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
    StorablePicture *dec_picture = currMB->p_Slice->dec_picture;
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;
    mb_t *neighborMB = NULL;

    int mb_x, mb_y;
    int a = 0, b = 0;
    int curr_cbp_ctx;
    int cbp = 0;
    int cbp_bit;
    int mask;
    PixelPos block_a;

    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    //  coding of luma part (bit by bit)
    for (mb_y = 0; mb_y < 4; mb_y += 2) {
        if (mb_y == 0) {
            neighborMB = currMB->mb_up;
            b = 0;
        }

        for (mb_x = 0; mb_x < 4; mb_x += 2) {
            if (mb_y == 0) {
                if (neighborMB != NULL) {
                    if (neighborMB->mb_type!=IPCM)
                        b = (( (neighborMB->cbp & (1<<(2 + (mb_x>>1)))) == 0) ? 2 : 0);
                }
            } else
                b = ( ((cbp & (1<<(mb_x/2))) == 0) ? 2: 0);

            if (mb_x == 0) {
                get4x4Neighbour(currMB, (mb_x<<2) - 1, (mb_y << 2), mb_size, &block_a);
                if (block_a.available) {
                    if (currSlice->mb_data[block_a.mb_addr].mb_type==IPCM)
                        a = 0;
                    else
                        a = (( (currSlice->mb_data[block_a.mb_addr].cbp & (1<<(2*(block_a.y/2)+1))) == 0) ? 1 : 0);
                } else
                    a = 0;
            } else
                a = ( ((cbp & (1<<mb_y)) == 0) ? 1: 0);

            curr_cbp_ctx = a + b;
            mask = (1 << (mb_y + (mb_x >> 1)));
            cbp_bit = dep_dp->decode_decision(ctx->cbp_contexts[0] + curr_cbp_ctx);
            if (cbp_bit) 
                cbp += mask;
        }
    }

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
        // coding of chroma part
        // CABAC decoding for BinIdx 0
        b = 0;
        neighborMB = currMB->mb_up;
        if (neighborMB != NULL) {
            if (neighborMB->mb_type==IPCM || (neighborMB->cbp > 15))
                b = 2;
        }

        a = 0;
        neighborMB = currMB->mb_left;
        if (neighborMB != NULL) {
            if (neighborMB->mb_type==IPCM || (neighborMB->cbp > 15))
                a = 1;
        }

        curr_cbp_ctx = a + b;
        cbp_bit = dep_dp->decode_decision(ctx->cbp_contexts[1] + curr_cbp_ctx);

        // CABAC decoding for BinIdx 1
        if (cbp_bit) { // set the chroma bits
            b = 0;
            neighborMB = currMB->mb_up;
            if (neighborMB != NULL) {
                if ((neighborMB->mb_type == IPCM) || ((neighborMB->cbp >> 4) == 2))
                    b = 2;
            }

            a = 0;
            neighborMB = currMB->mb_left;
            if (neighborMB != NULL) {
                if ((neighborMB->mb_type == IPCM) || ((neighborMB->cbp >> 4) == 2))
                    a = 1;
            }

            curr_cbp_ctx = a + b;
            cbp_bit = dep_dp->decode_decision(ctx->cbp_contexts[2] + curr_cbp_ctx);
            cbp += (cbp_bit == 1) ? 32 : 16;
        }
    }

    se->value1 = cbp;

    if (!cbp)
        currSlice->last_dquant = 0;
}

static void read_dQuant_CABAC(mb_t *currMB, syntax_element_t *se, cabac_engine_t *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    cabac_contexts_t *ctx = currSlice->mot_ctx;

    int dquant = 0;
    int act_ctx = (currSlice->last_dquant != 0 ? 1 : 0);
    int act_sym = dep_dp->decode_decision(ctx->delta_qp_contexts + act_ctx);

    //int binIdx[] = { act_ctx, 2, 3, 3, 3, 3, 3 };

    if (act_sym != 0) {
        act_ctx = 2;
        act_sym = dep_dp->u(ctx->delta_qp_contexts + act_ctx, 1) + 1;

        dquant = (act_sym + 1) >> 1;
        if ((act_sym & 1) == 0) // lsb is signed bit
            dquant = -dquant;
    }

    se->value1 = dquant;

    currSlice->last_dquant = dquant;
}


int check_next_mb_and_get_field_mode_CABAC(slice_t *slice)
{
    VideoParameters     *p_Vid = slice->p_Vid;
    cabac_context_t*       mb_type_ctx_copy[3];
    cabac_context_t*       mb_aff_ctx_copy;
    cabac_engine_t *dep_dp_copy;
    cabac_contexts_t  *mot_ctx  = slice->mot_ctx;  

    int length;

    int skip   = 0;
    int field  = 0;
    int i;

    mb_t *currMB;

    syntax_element_t currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    cabac_engine_t *dep_dp = &dP->bitstream->de_cabac;

    //get next MB
    ++slice->current_mb_nr;
  
    currMB = &slice->mb_data[slice->current_mb_nr];
    currMB->p_Vid    = p_Vid;
    currMB->p_Slice  = slice; 
    currMB->slice_nr = slice->current_slice_nr;
    currMB->mb_field_decoding_flag = slice->mb_data[slice->current_mb_nr-1].mb_field_decoding_flag;
    currMB->mbAddrX  = slice->current_mb_nr;

    CheckAvailabilityOfNeighborsMBAFF(currMB);
    CheckAvailabilityOfNeighborsCABAC(currMB);

    //create
    dep_dp_copy = (cabac_engine_t *)calloc(1, sizeof(cabac_engine_t) );
    for (i=0;i<3;++i)
        mb_type_ctx_copy[i] = (cabac_context_t*) calloc(NUM_MB_TYPE_CTX, sizeof(cabac_context_t) );
    mb_aff_ctx_copy = (cabac_context_t*) calloc(NUM_MB_AFF_CTX, sizeof(cabac_context_t) );

    //copy
    memcpy(dep_dp_copy,dep_dp,sizeof(cabac_engine_t));
    length = *(dep_dp_copy->Dcodestrm_len) = *(dep_dp->Dcodestrm_len);
    for (i=0;i<3;++i)
        memcpy(mb_type_ctx_copy[i], mot_ctx->mb_type_contexts[i],NUM_MB_TYPE_CTX*sizeof(cabac_context_t) );
    memcpy(mb_aff_ctx_copy, mot_ctx->mb_aff_contexts,NUM_MB_AFF_CTX*sizeof(cabac_context_t) );

    //check_next_mb
    slice->last_dquant = 0;
    read_skip_flag_CABAC(currMB, &currSE, dep_dp);

    skip = (currSE.value1 == 0);

    if (!skip) {
        readFieldModeInfo_CABAC(currMB, &currSE, dep_dp);
        field = currSE.value1;
        slice->mb_data[slice->current_mb_nr-1].mb_field_decoding_flag = (Boolean)field;
    }

    //reset
    slice->current_mb_nr--;

    memcpy(dep_dp,dep_dp_copy,sizeof(cabac_engine_t));
    *(dep_dp->Dcodestrm_len) = length;
    for (i=0;i<3;++i)
        memcpy(mot_ctx->mb_type_contexts[i],mb_type_ctx_copy[i], NUM_MB_TYPE_CTX*sizeof(cabac_context_t) );
    memcpy( mot_ctx->mb_aff_contexts,mb_aff_ctx_copy,NUM_MB_AFF_CTX*sizeof(cabac_context_t) );

    CheckAvailabilityOfNeighborsCABAC(currMB);

    //delete
    free(dep_dp_copy);
    for (i=0;i<3;++i)
        free(mb_type_ctx_copy[i]);
    free(mb_aff_ctx_copy);

    return skip;
}




uint32_t parse_mb_skip_run(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)
        currSE.value1 = dP->bitstream->ue();

    return currSE.value1;
}

bool parse_mb_skip_flag(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)
        currSE.value1 = dP->bitstream->ue();
    else
        read_skip_flag_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    mb->ei_flag = 0;

    return !currSE.value1;
}

bool parse_mb_field_decoding_flag(mb_t *mb)
{
	slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)
        currSE.value1 = dP->bitstream->f(1);
    else
        readFieldModeInfo_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    return currSE.value1;
}

uint32_t parse_mb_type(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)   
        currSE.value1 = dP->bitstream->ue();
    else {
        void (*reading)(macroblock_t*, syntax_element_t*, cabac_engine_t*) =
            slice->slice_type == I_slice || slice->slice_type == SI_slice ? readMB_typeInfo_CABAC_i_slice :
            slice->slice_type == P_slice || slice->slice_type == SP_slice ? readMB_typeInfo_CABAC_p_slice :
                                                                            readMB_typeInfo_CABAC_b_slice;
        reading(mb, &currSE, &dP->bitstream->de_cabac);
    }

    if (!pps->entropy_coding_mode_flag &&
        (slice->slice_type == P_slice || slice->slice_type == SP_slice))
        (currSE.value1)++;

    mb->ei_flag = 0;

    return currSE.value1;
}

uint8_t parse_sub_mb_type(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag) 
        currSE.value1 = dP->bitstream->ue();
    else
        readB8_typeInfo_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    return currSE.value1;
}


bool parse_transform_size_8x8_flag(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_HEADER;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)
        currSE.value1 = dP->bitstream->f(1);
    else
        readMB_transform_size_flag_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    return currSE.value1;
}

int8_t parse_intra_pred_mode(mb_t *mb, uint8_t block4x4Idx)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_INTRAPREDMODE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag) {
        if (dP->bitstream->f(1))
            currSE.value1 = -1;
        else
            currSE.value1 = dP->bitstream->f(3);
    } else {
        currSE.context = block4x4Idx;
        readIntraPredMode_CABAC(mb, &currSE, &dP->bitstream->de_cabac);
    }

    return currSE.value1;
}

uint8_t parse_intra_chroma_pred_mode(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = SE_INTRAPREDMODE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)
        currSE.value1 = dP->bitstream->ue();
    else
        readCIPredMode_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    return currSE.value1;
}

uint8_t parse_ref_idx(mb_t *mb, uint8_t b8mode, uint8_t list)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    bool refidx_present = (slice->slice_type == B_slice) ||
                          (!slice->allrefzero) ||
                          (mb->mb_type != P8x8);
    int num_ref_idx_active = list == LIST_0 ?
        slice->num_ref_idx_l0_active_minus1 + 1 :
        slice->num_ref_idx_l1_active_minus1 + 1;

    if (!refidx_present || num_ref_idx_active <= 1)
        return 0;

    syntax_element_t currSE;
    currSE.type = SE_REFFRAME;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

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

    syntax_element_t currSE;
    currSE.type = SE_MVD;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag) 
        currSE.value1 = dP->bitstream->se();
    else {
        currSE.value2 = xy * 2 + list;
        read_mvd_CABAC(mb, &currSE, &dP->bitstream->de_cabac);
    }

    return currSE.value1;
}

//! gives CBP value from codeword number, both for intra and inter
static const byte NCBP[2][48][2] = {
      // 0      1        2       3       4       5       6       7       8       9      10      11
    {{15, 0},{ 0, 1},{ 7, 2},{11, 4},{13, 8},{14, 3},{ 3, 5},{ 5,10},{10,12},{12,15},{ 1, 7},{ 2,11},
     { 4,13},{ 8,14},{ 6, 6},{ 9, 9},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},
     { 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},
     { 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0},{ 0, 0}},
    {{47, 0},{31,16},{15, 1},{ 0, 2},{23, 4},{27, 8},{29,32},{30, 3},{ 7, 5},{11,10},{13,12},{14,15},
     {39,47},{43, 7},{45,11},{46,13},{16,14},{ 3, 6},{ 5, 9},{10,31},{12,35},{19,37},{21,42},{26,44},
     {28,33},{35,34},{37,36},{42,40},{44,39},{ 1,43},{ 2,45},{ 4,46},{ 8,17},{17,18},{18,20},{20,24},
     {24,19},{ 6,21},{ 9,26},{22,28},{25,23},{32,27},{33,29},{34,30},{36,22},{40,25},{38,38},{41,41}}
};

uint8_t parse_coded_block_pattern(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    sps_t *sps = slice->active_sps;
    pps_t *pps = slice->active_pps;

    syntax_element_t currSE;
    currSE.type = mb->is_intra_block ? SE_CBP_INTRA : SE_CBP_INTER;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag) {
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

    syntax_element_t currSE;
    currSE.type = !mb->is_intra_block ? SE_DELTA_QUANT_INTER : SE_DELTA_QUANT_INTRA;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag)
        currSE.value1 = dP->bitstream->se();
    else
        read_dQuant_CABAC(mb, &currSE, &dP->bitstream->de_cabac);

    return currSE.value1;
}
