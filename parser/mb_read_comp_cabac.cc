/*!
 ***********************************************************************
 * \file read_comp_cabac.c
 *
 * \brief
 *     Read Coefficient Components
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

#include "global.h"
#include "slice.h"
#include "bitstream_elements.h"
#include "bitstream_cabac.h"
#include "bitstream.h"
#include "macroblock.h"
#include "mb_read.h"
#include "neighbour.h"
#include "transform.h"

#include "mb_read_syntax.h"


#define IS_I16MB(MB)    ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)
#define IS_DIRECT(MB)   ((MB)->mb_type == 0 && (currSlice->slice_type == B_SLICE))

static inline int rshift_rnd_sf(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}



//! for the linfo_levrun_inter routine
static const byte NTAB1[4][8][2] = {
    { {1, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0} },
    { {1, 1}, {1, 2}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0} },
    { {2, 0}, {1, 3}, {1, 4}, {1, 5}, {0, 0}, {0, 0}, {0, 0}, {0, 0} },
    { {3, 0}, {2, 1}, {2, 2}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {4, 0} }
};

static const byte LEVRUN1[16] = {
    4, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0
};

//! for the linfo_levrun__c2x2 routine
static const byte LEVRUN3[4] = { 2, 1, 0, 0 };

static const byte NTAB3[2][2][2] = {
    { {1, 0}, {0, 0} },
    { {2, 0}, {1, 1} }
};


static void linfo_levrun_inter(int len, int info, int *level, int *irun)
{
    if (len <= 9) {
        int l2  = imax(0,(len >> 1)-1);
        int inf = info >> 1;
        *level = NTAB1[l2][inf][0];
        *irun  = NTAB1[l2][inf][1];
        if (info % 2 == 1)
            *level = -*level;
    } else {
        *irun  = (info & 0x1e) >> 1;
        *level = LEVRUN1[*irun] + (info >> 5) + (1 << ((len >> 1) - 5));
        if (info % 2 == 1)
            *level = -*level;
    }

    if (len == 1)
        *level = 0;
}

static void linfo_levrun_c2x2(int len, int info, int *level, int *irun)
{
    if (len <= 5) {
        int l2  = imax(0, (len >> 1) - 1);
        int inf = info >> 1;
        *level = NTAB3[l2][inf][0];
        *irun  = NTAB3[l2][inf][1];
        if (info % 2 == 1)
            *level = -*level;
    } else {
        *irun  = (info & 0x06) >> 1;
        *level = LEVRUN3[*irun] + (info >> 3) + (1 << ((len >> 1) - 3));
        if (info % 2 == 1)
            *level = -*level;
    }

    if (len == 1)
        *level = 0;
}

static const short maxpos       [] = {15, 14, 63, 31, 31, 15,  3, 14,  7, 15, 15, 14, 63, 31, 31, 15, 15, 14, 63, 31, 31, 15};
static const short c1isdc       [] = { 1,  0,  1,  1,  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1};
static const short type2ctx_bcbp[] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 10, 11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20};
static const short type2ctx_map [] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}; // 8
static const short type2ctx_last[] = { 0,  1,  2,  3,  4,  5,  6,  7,  6,  6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}; // 8
static const short type2ctx_one [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 10, 11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20}; // 7
static const short type2ctx_abs [] = { 0,  1,  2,  3,  3,  4,  5,  6,  5,  5, 10, 11, 12, 13, 13, 14, 16, 17, 18, 19, 19, 20}; // 7
static const short max_c2       [] = { 4,  4,  4,  4,  4,  4,  3,  4,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4}; // 9

//===== position -> ctx for MAP =====
//--- zig-zag scan ----
static const byte  pos2ctx_map8x8 [] = { 0,  1,  2,  3,  4,  5,  5,  4,  4,  3,  3,  4,  4,  4,  5,  5,
                                         4,  4,  4,  4,  3,  3,  6,  7,  7,  7,  8,  9, 10,  9,  8,  7,
                                         7,  6, 11, 12, 13, 11,  6,  7,  8,  9, 14, 10,  9,  8,  6, 11,
                                        12, 13, 11,  6,  9, 14, 10,  9, 11, 12, 13, 11 ,14, 10, 12, 14}; // 15 CTX
static const byte  pos2ctx_map8x4 [] = { 0,  1,  2,  3,  4,  5,  7,  8,  9, 10, 11,  9,  8,  6,  7,  8,
                                         9, 10, 11,  9,  8,  6, 12,  8,  9, 10, 11,  9, 13, 13, 14, 14}; // 15 CTX
static const byte  pos2ctx_map4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14}; // 15 CTX
static const byte  pos2ctx_map2x4c[] = { 0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
static const byte  pos2ctx_map4x4c[] = { 0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
static const byte* pos2ctx_map    [] = {pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8, pos2ctx_map8x4,
                                        pos2ctx_map8x4, pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
                                        pos2ctx_map2x4c, pos2ctx_map4x4c, 
                                        pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8,pos2ctx_map8x4,
                                        pos2ctx_map8x4, pos2ctx_map4x4,
                                        pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8,pos2ctx_map8x4,
                                        pos2ctx_map8x4,pos2ctx_map4x4};
//--- interlace scan ----
//taken from ABT
static const byte  pos2ctx_map8x8i[] = { 0,  1,  1,  2,  2,  3,  3,  4,  5,  6,  7,  7,  7,  8,  4,  5,
                                         6,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 11, 12, 11,
                                         9,  9, 10, 10,  8, 11, 12, 11,  9,  9, 10, 10,  8, 13, 13,  9,
                                         9, 10, 10,  8, 13, 13,  9,  9, 10, 10, 14, 14, 14, 14, 14, 14}; // 15 CTX
static const byte  pos2ctx_map8x4i[] = { 0,  1,  2,  3,  4,  5,  6,  3,  4,  5,  6,  3,  4,  7,  6,  8,
                                         9,  7,  6,  8,  9, 10, 11, 12, 12, 10, 11, 13, 13, 14, 14, 14}; // 15 CTX
static const byte  pos2ctx_map4x8i[] = { 0,  1,  1,  1,  2,  3,  3,  4,  4,  4,  5,  6,  2,  7,  7,  8,
                                         8,  8,  5,  6,  9, 10, 10, 11, 11, 11, 12, 13, 13, 14, 14, 14}; // 15 CTX
static const byte* pos2ctx_map_int[] = {pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i,pos2ctx_map8x4i,
                                        pos2ctx_map4x8i,pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map4x4,
                                        pos2ctx_map2x4c, pos2ctx_map4x4c,
                                        pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i,pos2ctx_map8x4i,
                                        pos2ctx_map8x4i,pos2ctx_map4x4,
                                        pos2ctx_map4x4, pos2ctx_map4x4, pos2ctx_map8x8i,pos2ctx_map8x4i,
                                        pos2ctx_map8x4i,pos2ctx_map4x4};

//===== position -> ctx for LAST =====
static const byte  pos2ctx_last8x8 [] = { 0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                                          2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
                                          3,  3,  3,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,
                                          5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8}; //  9 CTX
static const byte  pos2ctx_last8x4 [] = { 0,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,
                                          3,  3,  3,  3,  4,  4,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8}; //  9 CTX

static const byte  pos2ctx_last4x4 [] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15}; // 15 CTX
static const byte  pos2ctx_last2x4c[] = { 0,  0,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
static const byte  pos2ctx_last4x4c[] = { 0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2}; // 15 CTX
static const byte* pos2ctx_last    [] = {pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8, pos2ctx_last8x4,
                                         pos2ctx_last8x4, pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last4x4,
                                         pos2ctx_last2x4c, pos2ctx_last4x4c,
                                         pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8,pos2ctx_last8x4,
                                         pos2ctx_last8x4, pos2ctx_last4x4,
                                         pos2ctx_last4x4, pos2ctx_last4x4, pos2ctx_last8x8,pos2ctx_last8x4,
                                         pos2ctx_last8x4, pos2ctx_last4x4};

static inline int get_bit(int64 x, int n)
{
    return (int)(((x >> n) & 1));
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

static unsigned int unary_exp_golomb_level_decode(DecodingEnvironment *dep_dp, BiContextTypePtr ctx)
{
    unsigned int symbol = biari_decode_symbol(dep_dp, ctx);

    if (symbol == 0)
        return 0;
    else {
        unsigned int l, k = 1;
        unsigned int exp_start = 13;

        symbol = 0;

        do {
            l = biari_decode_symbol(dep_dp, ctx);
            ++symbol;
            ++k;
        } while (l != 0 && k != exp_start);
        if (l != 0)
            symbol += exp_golomb_decode_eq_prob(dep_dp,0)+1;
        return symbol;
    }
}

static int read_significance_map(mb_t *currMB, DecodingEnvironment *dep_dp, int type, int coeff[])
{
    slice_t *currSlice = currMB->p_Slice;
    int               fld    = (currSlice->field_pic_flag || currMB->mb_field_decoding_flag);
    const byte *pos2ctx_Map = (fld) ? pos2ctx_map_int[type] : pos2ctx_map[type];
    const byte *pos2ctx_Last = pos2ctx_last[type];

    BiContextTypePtr  map_ctx  = currSlice->tex_ctx->map_contexts [fld][type2ctx_map [type]];
    BiContextTypePtr  last_ctx = currSlice->tex_ctx->last_contexts[fld][type2ctx_last[type]];

    int i;
    int coeff_ctr = 0;
    int i0        = 0;
    int i1        = maxpos[type];

    if (!c1isdc[type]) {
        ++i0;
        ++i1;
    }

    for (i = i0; i < i1; ++i) { // if last coeff is reached, it has to be significant
        //--- read significance symbol ---
        if (biari_decode_symbol(dep_dp, map_ctx + pos2ctx_Map[i])) {
            *(coeff++) = 1;
            ++coeff_ctr;
            //--- read last coefficient symbol ---
            if (biari_decode_symbol(dep_dp, last_ctx + pos2ctx_Last[i])) {
                memset(coeff, 0, (i1 - i) * sizeof(int));
                return coeff_ctr;
            }
        } else
            *(coeff++) = 0;
    }
    //--- last coefficient must be significant if no last symbol was received ---
    if (i < i1 + 1) {
        *coeff = 1;
        ++coeff_ctr;
    }

    return coeff_ctr;
}

static void read_significant_coefficients(DecodingEnvironment *dep_dp, TextureInfoContexts *tex_ctx, int type, int *coeff)
{
    BiContextType *one_contexts = tex_ctx->one_contexts[type2ctx_one[type]];
    BiContextType *abs_contexts = tex_ctx->abs_contexts[type2ctx_abs[type]];
    const short max_type = max_c2[type];
    int i = maxpos[type];
    int *cof = coeff + i;
    int   c1 = 1;
    int   c2 = 0;

    for (; i >= 0; i--) {
        if (*cof != 0) {
            *cof += biari_decode_symbol(dep_dp, one_contexts + c1);

            if (*cof == 2) {
                *cof += unary_exp_golomb_level_decode(dep_dp, abs_contexts + c2);
                c2 = imin (++c2, max_type);
                c1 = 0;
            } else if (c1)
                c1 = imin (++c1, 4);

            if (biari_decode_symbol_eq_prob(dep_dp))
                *cof = - *cof;
        }
        cof--;
    }
}

static int read_and_store_CBP_block_bit_444(mb_t *currMB, DecodingEnvironment *dep_dp, int type)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    TextureInfoContexts *tex_ctx = currSlice->tex_ctx;
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

    if (sps->separate_colour_plane_flag && type != LUMA_8x8) {
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
    } else {
        if (block_b.available) {
            if (!size_8x8_flag || mb_data[block_b.mb_addr].transform_size_8x8_flag) {
                if (mb_data[block_b.mb_addr].mb_type == IPCM)
                    upper_bit = 1;
                else if (size_8x8_flag)
                    upper_bit = get_bit(mb_data[block_b.mb_addr].s_cbp[pl].bits_8x8, bit + bit_pos_b);
                else
                    upper_bit = get_bit(mb_data[block_b.mb_addr].s_cbp[pl].bits, bit + bit_pos_b);
            }
        }
        if (block_a.available) {
            if (!size_8x8_flag || mb_data[block_a.mb_addr].transform_size_8x8_flag) {
                if (mb_data[block_a.mb_addr].mb_type == IPCM)
                    left_bit = 1;
                else if (size_8x8_flag)
                    left_bit = get_bit(mb_data[block_a.mb_addr].s_cbp[pl].bits_8x8, bit + bit_pos_a);
                else
                    left_bit = get_bit(mb_data[block_a.mb_addr].s_cbp[pl].bits, bit + bit_pos_a);
            }
        }
    }

    ctx = 2 * upper_bit + left_bit;
    cbp_bit = biari_decode_symbol(dep_dp, tex_ctx->bcbp_contexts[type2ctx_bcbp[type]] + ctx);

    bit = y_dc ?  0 :
          y_ac ?  1 + j + (i >> 2) :
          u_dc ? 17 :
          v_dc ? 18 :
          u_ac ? 19 + j + (i >> 2) :
                 35 + j + (i >> 2);

    if (cbp_bit) {
        if (size_8x8_flag) {
            currMB->s_cbp[pl].bits     |= (int64)(0x33 << bit);
            currMB->s_cbp[pl].bits_8x8 |= (int64)(0x33 << bit);
        } else if (type == LUMA_8x4 || type == CB_8x4 || type == CR_8x4)
            currMB->s_cbp[pl].bits |= (int64)(0x03 << bit);
        else if (type == LUMA_4x8 || type == CB_4x8 || type == CR_4x8)
            currMB->s_cbp[pl].bits |= (int64)(0x11 << bit);
        else
            currMB->s_cbp[pl].bits |= (int64)(0x01 << bit);
    }
    return cbp_bit;
}

static int read_and_store_CBP_block_bit_normal(mb_t *currMB, DecodingEnvironment *dep_dp, int type)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    TextureInfoContexts *tex_ctx = currSlice->tex_ctx;
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

    if (type != LUMA_8x8) {
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
    }

    if (type != LUMA_8x8) {
        ctx = 2 * upper_bit + left_bit;     
        cbp_bit = biari_decode_symbol(dep_dp, tex_ctx->bcbp_contexts[type2ctx_bcbp[type]] + ctx);
    }

    bit = y_dc ?  0 :
          y_ac ?  1 + j + (i >> 2) :
          u_dc ? 17 :
          v_dc ? 18 :
          u_ac ? 19 + j + (i >> 2) :
                 35 + j + (i >> 2);

    if (cbp_bit) {
        if (type == LUMA_8x8)
            currMB->s_cbp[0].bits |= (int64)(0x33 << bit);
        else if (type == LUMA_8x4)
            currMB->s_cbp[0].bits |= (int64)(0x03 << bit);
        else if (type == LUMA_4x8)
            currMB->s_cbp[0].bits |= (int64)(0x11 << bit);
        else
            currMB->s_cbp[0].bits |= i64_power2(bit);
    }
    return cbp_bit;
}

static int read_and_store_CBP_block_bit(mb_t *currMB, DecodingEnvironment *dep_dp, int type)
{
    if (currMB->p_Slice->active_sps->chroma_format_idc == YUV444)
        return read_and_store_CBP_block_bit_444(currMB, dep_dp, type);
    else
        return read_and_store_CBP_block_bit_normal(currMB, dep_dp, type);
}

static void readRunLevel_CABAC(mb_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    int  *coeff_ctr = &currSlice->coeff_ctr;
    int  *coeff = currSlice->coeff;

    if (*coeff_ctr < 0) {
        if ((*coeff_ctr = read_and_store_CBP_block_bit(currMB, dep_dp, se->context) ) != 0) {
            *coeff_ctr = read_significance_map(currMB, dep_dp, se->context, coeff);
            read_significant_coefficients(dep_dp, currSlice->tex_ctx, se->context, coeff);
        }
    }

    if (*coeff_ctr) {
        for (se->value2 = 0; coeff[currSlice->pos] == 0; ++currSlice->pos, ++se->value2);
        se->value1 = coeff[currSlice->pos++];
    } else {
        se->value1 = se->value2 = 0;
    }
    if ((*coeff_ctr)-- == 0) 
        currSlice->pos = 0;
}


static void read_tc_luma(mb_t *currMB, ColorPlane pl)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    const byte (*pos_scan4x4)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    const byte (*pos_scan8x8)[2] = !currSlice->field_pic_flag && !currMB->mb_field_decoding_flag ? SNGL_SCAN8x8 : FIELD_SCAN8x8;

    if (IS_I16MB(currMB) && !currMB->dpl_flag) {
        SyntaxElement currSE;
        currSE.type    = SE_LUM_DC_INTRA;
        currSE.context = LUMA_16DC;
        DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][currSE.type]];
        if (dP->bitstream->ei_flag)
            currSE.mapping = linfo_levrun_inter;
        else
            currSE.reading = readRunLevel_CABAC;

        int coef_ctr = -1;
        int level = 1;
        for (int k = 0; k < 17 && level != 0; ++k) {
            dP->readSyntaxElement(currMB, &currSE, dP);
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

        SyntaxElement currSE;
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
                        currSE.type = currMB->is_intra_block ?
                            (k == 0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) :
                            (k == 0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER);
                        DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][currSE.type]];
                        if (dP->bitstream->ei_flag)
                            currSE.mapping = linfo_levrun_inter;
                        else
                            currSE.reading = readRunLevel_CABAC;

                        dP->readSyntaxElement(currMB, &currSE, dP);
                        level = currSE.value1;
                        if (level != 0) {
                            coef_ctr += currSE.value2 + 1;
                            currMB->s_cbp[pl].blk |= (int64)(0x01 << (j * 4 + i));
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

        SyntaxElement currSE;
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
                    currSE.type = currMB->is_intra_block ?
                        (k == 0 ? SE_LUM_DC_INTRA : SE_LUM_AC_INTRA) :
                        (k == 0 ? SE_LUM_DC_INTER : SE_LUM_AC_INTER);
                    DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][currSE.type]];
                    if (dP->bitstream->ei_flag)
                        currSE.mapping = linfo_levrun_inter;
                    else
                        currSE.reading = readRunLevel_CABAC;

                    dP->readSyntaxElement(currMB, &currSE, dP);
                    level = currSE.value1;
                    if (level != 0) {
                        coef_ctr += currSE.value2 + 1;
                        currMB->s_cbp[pl].blk |= (int64)(0x33 << (block_y * 4 + block_x));
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
        SyntaxElement currSE;
        currSE.context = sps->ChromaArrayType == 1 ? CHROMA_DC : CHROMA_DC_2x4;
        currSE.type    = currMB->is_intra_block ? SE_CHR_DC_INTRA : SE_CHR_DC_INTER;
        DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][currSE.type]];
        if (dP->bitstream->ei_flag)
            currSE.mapping = linfo_levrun_c2x2;
        else
            currSE.reading = readRunLevel_CABAC;

        for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
            currMB->is_v_block = iCbCr * 2;

            int coef_ctr = -1;
            int level = 1;
            for (int k = 0; k < NumC8x8 * 4 + 1 && level != 0; ++k) {
                dP->readSyntaxElement(currMB, &currSE, dP);
                level = currSE.value1;

                int i0, j0;
                if (level != 0) {
                    coef_ctr += currSE.value2 + 1;
                    assert(coef_ctr < NumC8x8 * 4);
                    if (sps->ChromaArrayType == 1) {
                        i0 = coef_ctr % 2;
                        j0 = coef_ctr / 2;
                        currMB->s_cbp[0].blk |= (int64)(0xf << (iCbCr * 4 + 16));
                    }
                    if (sps->ChromaArrayType == 2) {
                        i0 = FIELD_SCAN[coef_ctr][0];
                        j0 = FIELD_SCAN[coef_ctr][1];
                        currMB->s_cbp[0].blk |= (int64)(0xff << (iCbCr * 8 + 16));
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
        SyntaxElement currSE;
        currSE.context = CHROMA_AC;
        currSE.type    = currMB->is_intra_block ? SE_CHR_AC_INTRA : SE_CHR_AC_INTER;
        DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][currSE.type]];
        if (dP->bitstream->ei_flag)
            currSE.mapping = linfo_levrun_inter;
        else
            currSE.reading = readRunLevel_CABAC;

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
                        dP->readSyntaxElement(currMB, &currSE, dP);
                        level = currSE.value1;

                        if (level != 0) {
                            currMB->s_cbp[0].blk |= (int64)(0x1 << (i8x8 * 4 + i4x4 + 16));
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

    SyntaxElement currSE;
    currSE.type    = SE_LUM_DC_INTRA;
    currSE.context = LUMA_16DC;
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][currSE.type]];
    if (dP->bitstream->ei_flag)
        currSE.mapping = linfo_levrun_inter;
    else
        currSE.reading = readRunLevel_CABAC;

    const byte (*pos_scan4x4)[2] = !slice->field_pic_flag && !this->mb_field_decoding_flag ? SNGL_SCAN : FIELD_SCAN;
    int coef_ctr = -1;
    int level = 1;
    for (int k = 0; k < 17 && level != 0; ++k) {
        dP->readSyntaxElement(this, &currSE, dP);
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
