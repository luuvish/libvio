
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
#include "biaridecod.h"
#include "transform.h"
#include "mv_prediction.h"
#include "intra_prediction.h"
#include "inter_prediction.h"

#include "mb_read_syntax.h"



static void read_skip_flag_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx  = currSlice->mot_ctx;
    int tabIdx = currSlice->slice_type == B_SLICE ? 2 : 1;
    int ctxIdx = currSlice->slice_type == B_SLICE ? 7 : 0;
    int a = (currMB->mb_left != NULL) ? (currMB->mb_left->mb_skip_flag == 0) : 0;
    int b = (currMB->mb_up   != NULL) ? (currMB->mb_up  ->mb_skip_flag == 0) : 0;
    int act_ctx = ctxIdx + a + b;

    se->value1 = biari_decode_symbol(dep_dp, &ctx->mb_type_contexts[tabIdx][act_ctx]) != 1;

    if (!se->value1)
        currMB->p_Slice->last_dquant = 0;
}

static void readFieldModeInfo_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx  = currSlice->mot_ctx;
    int a = currMB->mbAvailA ? currSlice->mb_data[currMB->mbAddrA].mb_field_decoding_flag : 0;
    int b = currMB->mbAvailB ? currSlice->mb_data[currMB->mbAddrB].mb_field_decoding_flag : 0;
    int act_ctx = a + b;

    se->value1 = biari_decode_symbol(dep_dp, &ctx->mb_aff_contexts[act_ctx]);
}

static void readMB_typeInfo_CABAC_i_slice(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;

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
        act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
        se->context = act_ctx; // store context

        if (act_sym == 0) // 4x4 Intra
            curr_mb_type = act_sym;
        else { // 16x16 Intra
            mode_sym = biari_decode_final(dep_dp);
            if (mode_sym == 1)
                curr_mb_type = 25;
            else {
                act_sym = 1;
                act_ctx = 4;
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx); // decoding of AC/no AC
                act_sym += mode_sym * 12;
                act_ctx = 5;
                // decoding of cbp: 0,1,2
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                if (mode_sym != 0) {
                    act_ctx = 6;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += 4;
                    if (mode_sym != 0)
                        act_sym += 4;
                }
                // decoding of I pred-mode: 0,1,2,3
                act_ctx = 7;
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                act_sym += mode_sym * 2;
                act_ctx = 8;
                mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
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
        act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[1] + act_ctx);
        se->context = act_ctx; // store context

        if (act_sym == 0) //  SI 4x4 Intra
            curr_mb_type = 0;
        else { // analog INTRA_IMG
            if (currMB->mb_up != NULL)
                b = (( (currMB->mb_up)->mb_type != I4MB) ? 1 : 0 );
            if (currMB->mb_left != NULL)
                a = (( (currMB->mb_left)->mb_type != I4MB) ? 1 : 0 );

            act_ctx = a + b;
            act_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
            se->context = act_ctx; // store context

            if (act_sym == 0) // 4x4 Intra
                curr_mb_type = 1;
            else { // 16x16 Intra
                mode_sym = biari_decode_final(dep_dp);
                if (mode_sym == 1)
                    curr_mb_type = 26;
                else {
                    act_sym = 2;
                    act_ctx = 4;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx); // decoding of AC/no AC
                    act_sym += mode_sym * 12;
                    act_ctx = 5;
                    // decoding of cbp: 0,1,2
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    if (mode_sym != 0) {
                        act_ctx = 6;
                        mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                        act_sym += 4;
                        if (mode_sym != 0)
                          act_sym += 4;
                    }
                    // decoding of I pred-mode: 0,1,2,3
                    act_ctx = 7;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += mode_sym * 2;
                    act_ctx = 8;
                    mode_sym = biari_decode_symbol(dep_dp, ctx->mb_type_contexts[0] + act_ctx);
                    act_sym += mode_sym;
                    curr_mb_type = act_sym;
                }
            }
        }
    }

    se->value1 = curr_mb_type;
}

static void readMB_typeInfo_CABAC_p_slice(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;

    int act_ctx;
    int act_sym;
    int mode_sym;
    int curr_mb_type;
    BiContextType *mb_type_contexts = ctx->mb_type_contexts[1];

    if (biari_decode_symbol(dep_dp, &mb_type_contexts[4])) {
        if (biari_decode_symbol(dep_dp, &mb_type_contexts[7]))
            act_sym = 7;
        else
            act_sym = 6;
    } else {
        if (biari_decode_symbol(dep_dp, &mb_type_contexts[5])) {
            if (biari_decode_symbol(dep_dp, &mb_type_contexts[7])) 
                act_sym = 2;
            else
                act_sym = 3;
        } else {
            if (biari_decode_symbol(dep_dp, &mb_type_contexts[6]))
                act_sym = 4;
            else
                act_sym = 1;
        }
    }

    if (act_sym <= 6)
        curr_mb_type = act_sym;
    else { // additional info for 16x16 Intra-mode
        mode_sym = biari_decode_final(dep_dp);
        if (mode_sym == 1)
            curr_mb_type = 31;
        else {
            act_ctx = 8;
            mode_sym =  biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx ); // decoding of AC/no AC
            act_sym += mode_sym*12;

            // decoding of cbp: 0,1,2
            act_ctx = 9;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
            if (mode_sym != 0) {
                act_sym += 4;
                mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
                if (mode_sym != 0)
                    act_sym += 4;
            }

            // decoding of I pred-mode: 0,1,2,3
            act_ctx = 10;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
            act_sym += mode_sym*2;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx );
            act_sym += mode_sym;
            curr_mb_type = act_sym;
        }
    }

    se->value1 = curr_mb_type;
}

static void readMB_typeInfo_CABAC_b_slice(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;

    int a = 0, b = 0;
    int act_ctx;
    int act_sym;
    int mode_sym;
    int curr_mb_type;
    BiContextType *mb_type_contexts = ctx->mb_type_contexts[2];

    if (currMB->mb_up != NULL)
        b = (( (currMB->mb_up)->mb_type != 0) ? 1 : 0 );

    if (currMB->mb_left != NULL)
        a = (( (currMB->mb_left)->mb_type != 0) ? 1 : 0 );

    act_ctx = a + b;

    if (biari_decode_symbol(dep_dp, &mb_type_contexts[act_ctx])) {
        if (biari_decode_symbol(dep_dp, &mb_type_contexts[4])) {
            if (biari_decode_symbol(dep_dp, &mb_type_contexts[5])) {
                act_sym = 12;
                if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                    act_sym += 8;
                if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                    act_sym += 4;
                if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                    act_sym += 2;

                if (act_sym == 24)  
                    act_sym = 11;
                else if (act_sym == 26)  
                    act_sym = 22;
                else {
                    if (act_sym == 22)
                        act_sym = 23;
                    if (biari_decode_symbol(dep_dp, &mb_type_contexts[6]))
                        act_sym += 1;
                }
            } else {
                act_sym = 3;
                if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                    act_sym += 4;
                if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                    act_sym += 2;
                if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                    act_sym += 1;
            }
        } else {
            if (biari_decode_symbol(dep_dp, &mb_type_contexts[6])) 
                act_sym = 2;
            else
                act_sym = 1;
        }
    } else
        act_sym = 0;

    if (act_sym <= 23)
        curr_mb_type = act_sym;
    else { // additional info for 16x16 Intra-mode
        mode_sym = biari_decode_final(dep_dp);
        if (mode_sym == 1)
            curr_mb_type = 48;
        else {
            mb_type_contexts = ctx->mb_type_contexts[1];
            act_ctx = 8;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx); // decoding of AC/no AC
            act_sym += mode_sym*12;

            // decoding of cbp: 0,1,2
            act_ctx = 9;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx);
            if (mode_sym != 0) {
                act_sym += 4;
                mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx);
                if (mode_sym != 0)
                    act_sym += 4;
            }

            // decoding of I pred-mode: 0,1,2,3
            act_ctx = 10;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx);
            act_sym += mode_sym * 2;
            mode_sym = biari_decode_symbol(dep_dp, mb_type_contexts + act_ctx);
            act_sym += mode_sym;
            curr_mb_type = act_sym;
        }
    }

    se->value1 = curr_mb_type;
}

static void readB8_typeInfo_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;
    BiContextType *b8_type_contexts;
    int act_sym = 0;

    if (currSlice->slice_type == P_slice) {
        b8_type_contexts = &ctx->b8_type_contexts[0][1];

        if (biari_decode_symbol(dep_dp, b8_type_contexts++))
            act_sym = 0;
        else {
            if (biari_decode_symbol(dep_dp, ++b8_type_contexts))
                act_sym = (biari_decode_symbol(dep_dp, ++b8_type_contexts)) ? 2: 3;
            else
                act_sym = 1;
        }
    } else {
        b8_type_contexts = ctx->b8_type_contexts[1];

        if (biari_decode_symbol(dep_dp, b8_type_contexts++)) {
            if (biari_decode_symbol(dep_dp, b8_type_contexts++)) {
                if (biari_decode_symbol(dep_dp, b8_type_contexts++)) {
                    if (biari_decode_symbol(dep_dp, b8_type_contexts)) {
                        act_sym = 10;
                        if (biari_decode_symbol(dep_dp, b8_type_contexts)) 
                            act_sym++;
                    } else {
                        act_sym = 6;
                        if (biari_decode_symbol(dep_dp, b8_type_contexts)) 
                            act_sym += 2;
                        if (biari_decode_symbol(dep_dp, b8_type_contexts)) 
                            act_sym++;
                    }
                } else {
                    act_sym = 2;
                    if (biari_decode_symbol(dep_dp, b8_type_contexts)) 
                        act_sym += 2;
                    if (biari_decode_symbol(dep_dp, b8_type_contexts)) 
                        act_sym++;
                }
            } else
                act_sym = biari_decode_symbol(dep_dp, ++b8_type_contexts) ? 1 : 0;
            act_sym++;
        } else
            act_sym = 0;
    }

    se->value1 = act_sym;
}


static void readMB_transform_size_flag_CABAC( mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    TextureInfoContexts *ctx = currSlice->tex_ctx;

    int b = (currMB->mb_up   == NULL) ? 0 : currMB->mb_up->transform_size_8x8_flag;
    int a = (currMB->mb_left == NULL) ? 0 : currMB->mb_left->transform_size_8x8_flag;

    int act_sym = biari_decode_symbol(dep_dp, ctx->transform_size_contexts + a + b );

    se->value1 = act_sym;
}

static void readIntraPredMode_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    TextureInfoContexts *ctx = currSlice->tex_ctx;

    int act_sym = biari_decode_symbol(dep_dp, ctx->ipr_contexts);

    // remaining_mode_selector
    if (act_sym == 1)
        se->value1 = -1;
    else {
        se->value1  = (biari_decode_symbol(dep_dp, ctx->ipr_contexts + 1)     );
        se->value1 |= (biari_decode_symbol(dep_dp, ctx->ipr_contexts + 1) << 1);
        se->value1 |= (biari_decode_symbol(dep_dp, ctx->ipr_contexts + 1) << 2);
    }
}

static unsigned int unary_bin_max_decode(DecodingEnvironment *dep_dp, BiContextTypePtr ctx, int ctx_offset, unsigned int max_symbol)
{
    unsigned int symbol =  biari_decode_symbol(dep_dp, ctx );

    if (symbol == 0 || max_symbol == 0)
        return symbol;
    else {
        unsigned int l;
        ctx += ctx_offset;
        symbol = 0;
        do {
            l = biari_decode_symbol(dep_dp, ctx);
            ++symbol;
        } while (l != 0 && symbol < max_symbol);

        if (l != 0 && symbol == max_symbol)
            ++symbol;
        return symbol;
    }
}

static void readCIPredMode_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    TextureInfoContexts *ctx = currSlice->tex_ctx;
    int                 *act_sym  = &se->value1;

    mb_t          *MbUp   = currMB->mb_up;
    mb_t          *MbLeft = currMB->mb_left;

    int b = (MbUp != NULL)   ? (((MbUp->intra_chroma_pred_mode   != 0) && (MbUp->mb_type != IPCM)) ? 1 : 0) : 0;
    int a = (MbLeft != NULL) ? (((MbLeft->intra_chroma_pred_mode != 0) && (MbLeft->mb_type != IPCM)) ? 1 : 0) : 0;
    int act_ctx = a + b;

    *act_sym = biari_decode_symbol(dep_dp, ctx->cipr_contexts + act_ctx );

    if (*act_sym != 0)
        *act_sym = unary_bin_max_decode(dep_dp, ctx->cipr_contexts + 3, 0, 1) + 1;
}

static unsigned int unary_bin_decode(DecodingEnvironment *dep_dp, BiContextTypePtr ctx, int ctx_offset)
{
    unsigned int symbol = biari_decode_symbol(dep_dp, ctx);

    if (symbol == 0)
        return 0;
    else {
        unsigned int l;
        ctx += ctx_offset;;
        symbol = 0;
        do {
            l = biari_decode_symbol(dep_dp, ctx);
            ++symbol;
        }
        while (l != 0);
        return symbol;
    }
}

static void readRefFrame_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    StorablePicture *dec_picture = currSlice->dec_picture;
    MotionInfoContexts *ctx = currSlice->mot_ctx;
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

    act_sym = biari_decode_symbol(dep_dp,ctx->ref_no_contexts[addctx] + act_ctx );

    if (act_sym != 0) {
        act_ctx = 4;
        act_sym = unary_bin_decode(dep_dp,ctx->ref_no_contexts[addctx] + act_ctx,1);
        ++act_sym;
    }
    se->value1 = act_sym;
}

static unsigned int exp_golomb_decode_eq_prob(DecodingEnvironment *dep_dp, int k)
{
    unsigned int l;
    int symbol = 0;
    int binary_symbol = 0;

    do {
        l = biari_decode_symbol_eq_prob(dep_dp);
        if (l == 1) {
            symbol += (1<<k);
            ++k;
        }
    } while (l != 0);

    while (k--)                             //next binary part
        if (biari_decode_symbol_eq_prob(dep_dp) == 1)
            binary_symbol |= (1<<k);

    return (unsigned int) (symbol + binary_symbol);
}

static unsigned int unary_exp_golomb_mv_decode(DecodingEnvironment *dep_dp, BiContextTypePtr ctx, unsigned int max_bin)
{
    unsigned int symbol = biari_decode_symbol(dep_dp, ctx );

    if (symbol == 0)
        return 0;
    else {
        unsigned int exp_start = 8;
        unsigned int l,k = 1;
        unsigned int bin = 1;

        symbol = 0;

        ++ctx;
        do {
            l = biari_decode_symbol(dep_dp, ctx);
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

static void read_MVD_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{  
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;
    int i = currMB->subblock_x;
    int j = currMB->subblock_y;
    int a = 0;
    //int act_ctx;
    int act_sym;  
    int list_idx = se->value2 & 0x01;
    int k = (se->value2 >> 1); // MVD component

    PixelPos block_a, block_b;
    get4x4NeighbourBase(currMB, i - 1, j    , mb_size, &block_a);
    get4x4NeighbourBase(currMB, i    , j - 1, mb_size, &block_b);

    if (block_a.available)
        a = iabs(currSlice->mb_data[block_a.mb_addr].mvd[list_idx][block_a.y][block_a.x][k]);
    if (block_b.available)
        a += iabs(currSlice->mb_data[block_b.mb_addr].mvd[list_idx][block_b.y][block_b.x][k]);

    if (a < 3)
        a = 5 * k;
    else if (a > 32)
        a = 5 * k + 3;
    else
        a = 5 * k + 2;

    se->context = a;

    act_sym = biari_decode_symbol(dep_dp, ctx->mv_res_contexts[0] + a );

    if (act_sym != 0) {
        a = 5 * k;
        act_sym = unary_exp_golomb_mv_decode(dep_dp, ctx->mv_res_contexts[1] + a, 3) + 1;

        if (biari_decode_symbol_eq_prob(dep_dp))
            act_sym = -act_sym;
    }
    se->value1 = act_sym;
}

static void read_mvd_CABAC_mbaff(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    MotionInfoContexts *ctx = currSlice->mot_ctx;
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

    if (block_a.available) {
        a = iabs(currSlice->mb_data[block_a.mb_addr].mvd[list_idx][block_a.y][block_a.x][k]);
        if (currSlice->MbaffFrameFlag && (k==1)) {
            if ((currMB->mb_field_decoding_flag==0) && (currSlice->mb_data[block_a.mb_addr].mb_field_decoding_flag==1))
                a *= 2;
            else if ((currMB->mb_field_decoding_flag==1) && (currSlice->mb_data[block_a.mb_addr].mb_field_decoding_flag==0))
                a /= 2;
        }
    }

    get4x4NeighbourBase(currMB, i    , j - 1, mb_size, &block_b);

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

    act_sym = biari_decode_symbol(dep_dp,&ctx->mv_res_contexts[0][act_ctx] );

    if (act_sym != 0) {
        act_ctx = 5 * k;
        act_sym = unary_exp_golomb_mv_decode(dep_dp, ctx->mv_res_contexts[1] + act_ctx, 3) + 1;

        if (biari_decode_symbol_eq_prob(dep_dp))
            act_sym = -act_sym;
    }
    se->value1 = act_sym;
}


int check_next_mb_and_get_field_mode_CABAC(slice_t *slice)
{
    VideoParameters     *p_Vid = slice->p_Vid;
    BiContextTypePtr     mb_type_ctx_copy[3];
    BiContextTypePtr     mb_aff_ctx_copy;
    DecodingEnvironment *dep_dp_copy;
    MotionInfoContexts  *mot_ctx  = slice->mot_ctx;  

    int length;

    int skip   = 0;
    int field  = 0;
    int i;

    mb_t *currMB;

    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];
    currSE.reading = read_skip_flag_CABAC;

    DecodingEnvironment *dep_dp = &dP->bitstream->de_cabac;

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
    dep_dp_copy = (DecodingEnvironment *)calloc(1, sizeof(DecodingEnvironment) );
    for (i=0;i<3;++i)
        mb_type_ctx_copy[i] = (BiContextTypePtr) calloc(NUM_MB_TYPE_CTX, sizeof(BiContextType) );
    mb_aff_ctx_copy = (BiContextTypePtr) calloc(NUM_MB_AFF_CTX, sizeof(BiContextType) );

    //copy
    memcpy(dep_dp_copy,dep_dp,sizeof(DecodingEnvironment));
    length = *(dep_dp_copy->Dcodestrm_len) = *(dep_dp->Dcodestrm_len);
    for (i=0;i<3;++i)
        memcpy(mb_type_ctx_copy[i], mot_ctx->mb_type_contexts[i],NUM_MB_TYPE_CTX*sizeof(BiContextType) );
    memcpy(mb_aff_ctx_copy, mot_ctx->mb_aff_contexts,NUM_MB_AFF_CTX*sizeof(BiContextType) );

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

    memcpy(dep_dp,dep_dp_copy,sizeof(DecodingEnvironment));
    *(dep_dp->Dcodestrm_len) = length;
    for (i=0;i<3;++i)
        memcpy(mot_ctx->mb_type_contexts[i],mb_type_ctx_copy[i], NUM_MB_TYPE_CTX*sizeof(BiContextType) );
    memcpy( mot_ctx->mb_aff_contexts,mb_aff_ctx_copy,NUM_MB_AFF_CTX*sizeof(BiContextType) );

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

    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
        currSE.mapping = linfo_ue;
        currSE.value1 = dP->bitstream->ue();
    }

    return currSE.value1;
}

bool parse_mb_skip_flag(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
        currSE.mapping = linfo_ue;
    else
        currSE.reading = read_skip_flag_CABAC;
    dP->readSyntaxElement(mb, &currSE, dP);

    if (!dP->bitstream->ei_flag)
        mb->ei_flag = 0;

    return !currSE.value1;
}

bool parse_mb_field_decoding_flag(mb_t *mb)
{
	slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
        currSE.mapping = linfo_ue;
        currSE.value1 = dP->bitstream->f(1);
    } else {
        currSE.reading = readFieldModeInfo_CABAC;
        dP->readSyntaxElement(mb, &currSE, dP);
    }

    return currSE.value1;
}

uint32_t parse_mb_type(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)   
        currSE.mapping = linfo_ue;
    else
        currSE.reading =
            slice->slice_type == I_slice || slice->slice_type == SI_slice ? readMB_typeInfo_CABAC_i_slice :
            slice->slice_type == P_slice || slice->slice_type == SP_slice ? readMB_typeInfo_CABAC_p_slice :
                                                                            readMB_typeInfo_CABAC_b_slice;
    dP->readSyntaxElement(mb, &currSE, dP);

    if (!pps->entropy_coding_mode_flag &&
        (slice->slice_type == P_slice || slice->slice_type == SP_slice))
        (currSE.value1)++;

    if (!dP->bitstream->ei_flag)
        mb->ei_flag = 0;

    return currSE.value1;
}

uint8_t parse_sub_mb_type(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) 
        currSE.mapping = linfo_ue;
    else
        currSE.reading = readB8_typeInfo_CABAC;

    dP->readSyntaxElement(mb, &currSE, dP);

    return currSE.value1;
}


bool parse_transform_size_8x8_flag(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_HEADER;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
        currSE.value1 = dP->bitstream->f(1);
    else {
        currSE.reading = readMB_transform_size_flag_CABAC;
        dP->readSyntaxElement(mb, &currSE, dP);
    }

    return currSE.value1;
}

int8_t parse_intra_pred_mode(mb_t *mb, uint8_t block4x4Idx)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_INTRAPREDMODE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
        if (dP->bitstream->f(1))
            currSE.value1 = -1;
        else
            currSE.value1 = dP->bitstream->f(3);
    } else {
        currSE.context = block4x4Idx;
        currSE.reading = readIntraPredMode_CABAC;
        dP->readSyntaxElement(mb, &currSE, dP);
    }

    return currSE.value1;
}

uint8_t parse_intra_chroma_pred_mode(mb_t *mb)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_INTRAPREDMODE;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
        currSE.mapping = linfo_ue;
    else
        currSE.reading = readCIPredMode_CABAC;

    dP->readSyntaxElement(mb, &currSE, dP);

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

    SyntaxElement currSE;
    currSE.type = SE_REFFRAME;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
        currSE.mapping = linfo_ue;
    else
        currSE.reading = readRefFrame_CABAC;

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
        if (num_ref_idx_active == 2) {
            currSE.context = (b8mode >= 4);
            currSE.value1 = 1 - dP->bitstream->f(1);
            return currSE.value1;
        }
    }

    currSE.context = (b8mode >= 4);
    currSE.value2 = list;
    dP->readSyntaxElement(mb, &currSE, dP);

    return currSE.value1;
}

int16_t parse_mvd(mb_t *mb, uint8_t xy, uint8_t list)
{
    slice_t *slice = mb->p_Slice;
    pps_t *pps = slice->active_pps;

    SyntaxElement currSE;
    currSE.type = SE_MVD;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];
    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) 
        currSE.mapping = linfo_se;
    else
        currSE.reading = slice->MbaffFrameFlag ? read_mvd_CABAC_mbaff : read_MVD_CABAC;

    currSE.value2 = xy * 2 + list;
    dP->readSyntaxElement(mb, &currSE, dP);

    return currSE.value1;
}
