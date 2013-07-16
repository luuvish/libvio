/*!
 ***********************************************************************
 * \file macroblock.c
 *
 * \brief
 *     Decode a mb_t
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Inge Lille-Langøy               <inge.lille-langoy@telenor.com>
 *    - Rickard Sjoberg                 <rickard.sjoberg@era.ericsson.se>
 *    - Jani Lainema                    <jani.lainema@nokia.com>
 *    - Sebastian Purreiter             <sebastian.purreiter@mch.siemens.de>
 *    - Thomas Wedi                     <wedi@tnt.uni-hannover.de>
 *    - Detlev Marpe
 *    - Gabi Blaettermann
 *    - Ye-Kui Wang                     <wyk@ieee.org>
 *    - Lowell Winger                   <lwinger@lsil.com>
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 ***********************************************************************
*/

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
#include "bitstream_ctx.h"


// Table 7-11 mb_t types for I slices
const uint8_t mb_types_i_slice[26][5] = {
    { I_NxN        , Intra_8x8  , NA,  0,  0 },
    { I_16x16_0_0_0, Intra_16x16,  0,  0,  0 },
    { I_16x16_1_0_0, Intra_16x16,  1,  0,  0 },
    { I_16x16_2_0_0, Intra_16x16,  2,  0,  0 },
    { I_16x16_3_0_0, Intra_16x16,  3,  0,  0 },
    { I_16x16_0_1_0, Intra_16x16,  0,  1,  0 },
    { I_16x16_1_1_0, Intra_16x16,  1,  1,  0 },
    { I_16x16_2_1_0, Intra_16x16,  2,  1,  0 },
    { I_16x16_3_1_0, Intra_16x16,  3,  1,  0 },
    { I_16x16_0_2_0, Intra_16x16,  0,  2,  0 },
    { I_16x16_1_2_0, Intra_16x16,  1,  2,  0 },
    { I_16x16_2_2_0, Intra_16x16,  2,  2,  0 },
    { I_16x16_3_2_0, Intra_16x16,  3,  2,  0 },
    { I_16x16_0_0_1, Intra_16x16,  0,  0, 15 },
    { I_16x16_1_0_1, Intra_16x16,  1,  0, 15 },
    { I_16x16_2_0_1, Intra_16x16,  2,  0, 15 },
    { I_16x16_3_0_1, Intra_16x16,  3,  0, 15 },
    { I_16x16_0_1_1, Intra_16x16,  0,  1, 15 },
    { I_16x16_1_1_1, Intra_16x16,  1,  1, 15 },
    { I_16x16_2_1_1, Intra_16x16,  2,  1, 15 },
    { I_16x16_3_1_1, Intra_16x16,  3,  1, 15 },
    { I_16x16_0_2_1, Intra_16x16,  0,  2, 15 },
    { I_16x16_1_2_1, Intra_16x16,  1,  2, 15 },
    { I_16x16_2_2_1, Intra_16x16,  2,  2, 15 },
    { I_16x16_3_2_1, Intra_16x16,  3,  2, 15 },
    { I_PCM,                  NA, NA, NA, NA }
};

// Table 7-13 mb_t type value 0 to 4 for P and SP slices
const uint8_t mb_types_p_slice[6][6] = {
    { P_L0_16x16  , 1, Pred_L0,      NA, 16, 16 },
    { P_L0_L0_16x8, 2, Pred_L0, Pred_L0, 16,  8 },
    { P_L0_L0_8x16, 2, Pred_L0, Pred_L0,  8, 16 },
    { P_8x8       , 4,      NA,      NA,  8,  8 },
    { P_8x8ref0   , 4,      NA,      NA,  8,  8 },
    { P_Skip      , 1, Pred_L0,      NA, 16, 16 }
};

// Table 7-14 mb_t type value 0 to 22 for B slices
const uint8_t mb_types_b_slice[24][6] = {
    { B_Direct_16x16, NA, Direct ,      NA,  8,  8 },
    { B_L0_16x16    ,  1, Pred_L0,      NA, 16, 16 },
    { B_L1_16x16    ,  1, Pred_L1,      NA, 16, 16 },
    { B_Bi_16x16    ,  1, BiPred ,      NA, 16, 16 },
    { B_L0_L0_16x8  ,  2, Pred_L0, Pred_L0, 16,  8 },
    { B_L0_L0_8x16  ,  2, Pred_L0, Pred_L0,  8, 16 },
    { B_L1_L1_16x8  ,  2, Pred_L1, Pred_L1, 16,  8 },
    { B_L1_L1_8x16  ,  2, Pred_L1, Pred_L1,  8, 16 },
    { B_L0_L1_16x8  ,  2, Pred_L0, Pred_L1, 16,  8 },
    { B_L0_L1_8x16  ,  2, Pred_L0, Pred_L1,  8, 16 },
    { B_L1_L0_16x8  ,  2, Pred_L1, Pred_L0, 16,  8 },
    { B_L1_L0_8x16  ,  2, Pred_L1, Pred_L0,  8, 16 },
    { B_L0_Bi_16x8  ,  2, Pred_L0, BiPred , 16,  8 },
    { B_L0_Bi_8x16  ,  2, Pred_L0, BiPred ,  8, 16 },
    { B_L1_Bi_16x8  ,  2, Pred_L1, BiPred , 16,  8 },
    { B_L1_Bi_8x16  ,  2, Pred_L1, BiPred ,  8, 16 },
    { B_Bi_L0_16x8  ,  2, BiPred , Pred_L0, 16,  8 },
    { B_Bi_L0_8x16  ,  2, BiPred , Pred_L0,  8, 16 },
    { B_Bi_L1_16x8  ,  2, BiPred , Pred_L1, 16,  8 },
    { B_Bi_L1_8x16  ,  2, BiPred , Pred_L1,  8, 16 },
    { B_Bi_Bi_16x8  ,  2, BiPred , BiPred , 16,  8 },
    { B_Bi_Bi_8x16  ,  2, BiPred , BiPred ,  8, 16 },
    { B_8x8         ,  4,      NA,      NA,  8,  8 },
    { B_Skip        , NA, Direct ,      NA,  8,  8 }
};

// Table 7-17 Sub-macroblock types in P macroblocks
const uint8_t sub_mb_types_p_slice[4][5] = {
    { P_L0_8x8, 1, Pred_L0, 8, 8 },
    { P_L0_8x4, 2, Pred_L0, 8, 4 },
    { P_L0_4x8, 2, Pred_L0, 4, 8 },
    { P_L0_4x4, 4, Pred_L0, 4, 4 }
};

// Table 7-18 Sub-macroblock types in B macroblocks
const uint8_t sub_mb_types_b_slice[14][5] = {
    { B_Direct_8x8, 4, Direct , 4, 4 },
    { B_L0_8x8    , 1, Pred_L0, 8, 8 },
    { B_L1_8x8    , 1, Pred_L1, 8, 8 },
    { B_Bi_8x8    , 1, BiPred , 8, 8 },
    { B_L0_8x4    , 2, Pred_L0, 8, 4 },
    { B_L0_4x8    , 2, Pred_L0, 4, 8 },
    { B_L1_8x4    , 2, Pred_L1, 8, 4 },
    { B_L1_4x8    , 2, Pred_L1, 4, 8 },
    { B_Bi_8x4    , 2, BiPred , 8, 4 },
    { B_Bi_4x8    , 2, BiPred , 4, 8 },
    { B_L0_4x4    , 4, Pred_L0, 4, 4 },
    { B_L1_4x4    , 4, Pred_L1, 4, 4 },
    { B_Bi_4x4    , 4, BiPred , 4, 4 },
    {           NA, 4, Direct , 4, 4 }
};


static void readB8_typeInfo_CABAC(mb_t *currMB, 
                                  SyntaxElement *se,
                                  DecodingEnvironment *dep_dp)
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

static void readIntraPredMode_CABAC(mb_t *currMB, 
                                    SyntaxElement *se,
                                    DecodingEnvironment *dep_dp)
{
    slice_t *currSlice = currMB->p_Slice;
    TextureInfoContexts *ctx = currSlice->tex_ctx;
    // use_most_probable_mode
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

static char readRefPictureIdx(mb_t *currMB, char b8mode, int list)
{
    slice_t *currSlice = currMB->p_Slice;
    pps_t *pps = currSlice->active_pps;

    bool refidx_present = (currSlice->slice_type == B_slice) ||
                          (!currSlice->allrefzero) ||
                          (currMB->mb_type != P8x8);
    int num_ref_idx_active = list == LIST_0 ?
        currSlice->num_ref_idx_l0_active_minus1 + 1 :
        currSlice->num_ref_idx_l1_active_minus1 + 1;

    if (!refidx_present || num_ref_idx_active <= 1)
        return 0;

    DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][SE_REFFRAME]];
    SyntaxElement currSE;
    currSE.type = SE_REFFRAME;
    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
        currSE.mapping = linfo_ue;
    else
        currSE.reading = readRefFrame_CABAC;

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
        if (num_ref_idx_active == 2) {
            currSE.context = (b8mode >= 4);
            currSE.value1 = 1 - dP->bitstream->f(1);
            return (char) currSE.value1;
        }
    }

    currSE.context = (b8mode >= 4);
    currSE.value2 = list;
    dP->readSyntaxElement(currMB, &currSE, dP);
    return (char) currSE.value1;
}




static inline void reset_mv_info(PicMotionParams *mv_info, int slice_no)
{
    mv_info->ref_pic[LIST_0] = NULL;
    mv_info->ref_pic[LIST_1] = NULL;
    mv_info->mv     [LIST_0] = zero_mv;
    mv_info->mv     [LIST_1] = zero_mv;
    mv_info->ref_idx[LIST_0] = -1;
    mv_info->ref_idx[LIST_1] = -1;
    mv_info->slice_no = slice_no;
}

static inline void reset_mv_info_list(PicMotionParams *mv_info, int list, int slice_no)
{
    mv_info->ref_pic[list] = NULL;
    mv_info->mv     [list] = zero_mv;
    mv_info->ref_idx[list] = -1;
    mv_info->slice_no = slice_no;
}


static void init_macroblock_basic(mb_t *currMB)
{
    int j, i;
    PicMotionParams **mv_info = &currMB->p_Slice->dec_picture->mv_info[currMB->block_y];
    int slice_no =  currMB->p_Slice->current_slice_nr;
    // reset vectors and pred. modes
    for (j = 0; j < BLOCK_SIZE; ++j) {
        i = currMB->block_x;
        reset_mv_info_list(*mv_info + (i++), LIST_1, slice_no);
        reset_mv_info_list(*mv_info + (i++), LIST_1, slice_no);
        reset_mv_info_list(*mv_info + (i++), LIST_1, slice_no);
        reset_mv_info_list(*(mv_info++) + i, LIST_1, slice_no);
    }
}

static void init_macroblock_direct(mb_t *currMB)
{
    int slice_no = currMB->p_Slice->current_slice_nr;
    PicMotionParams **mv_info = &currMB->p_Slice->dec_picture->mv_info[currMB->block_y]; 
    int i, j;

    i = currMB->block_x;
    for (j = 0; j < BLOCK_SIZE; ++j) {
        (*mv_info+i)->slice_no = slice_no;
        (*mv_info+i+1)->slice_no = slice_no;
        (*mv_info+i+2)->slice_no = slice_no;
        (*(mv_info++)+i+3)->slice_no = slice_no;
    }
}

static void init_macroblock(mb_t *currMB)
{
    int j, i;
    slice_t *currSlice = currMB->p_Slice;
    PicMotionParams **mv_info = &currSlice->dec_picture->mv_info[currMB->block_y]; 
    int slice_no = currSlice->current_slice_nr;
    // reset vectors and pred. modes

    for (j = 0; j < BLOCK_SIZE; ++j) {
        i = currMB->block_x;
        reset_mv_info(*mv_info + (i++), slice_no);
        reset_mv_info(*mv_info + (i++), slice_no);
        reset_mv_info(*mv_info + (i++), slice_no);
        reset_mv_info(*(mv_info++) + i, slice_no);
    }
}

static void concealIPCMcoeffs(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int i, j, k;

    for (i = 0; i < 16; ++i) {
        for (j = 0; j < 16; ++j)
        currSlice->cof[0][i][j] = (1 << (sps->BitDepthY - 1));
    }

    if (sps->chroma_format_idc != YUV400 && !sps->separate_colour_plane_flag) {
        for (k = 0; k < 2; ++k)
            for (i = 0; i < sps->MbHeightC; ++i)
                for (j = 0; j < sps->MbWidthC; ++j)
                    currSlice->cof[k][i][j] = (1 << (sps->BitDepthC - 1));
    }
}


void readIPCM_CABAC(slice_t *currSlice, struct datapartition_dec *dP)
{
    sps_t *sps = currSlice->active_sps;
    StorablePicture *dec_picture = currSlice->dec_picture;
    Bitstream* currStream = dP->bitstream;
    DecodingEnvironment *dep = &dP->bitstream->de_cabac;
    byte *buf = currStream->streamBuffer;
    int BitstreamLengthInBits = (dP->bitstream->bitstream_length << 3) + 7;

    int val = 0;
    int bits_read = 0;
    int bitoffset, bitdepth;
    int uv, i, j;

    while (dep->DbitsLeft >= 8) {
        dep->Dvalue   >>= 8;
        dep->DbitsLeft -= 8;
        (*dep->Dcodestrm_len)--;
    }

    bitoffset = (*dep->Dcodestrm_len) << 3;

    // read luma values
    bitdepth = sps->BitDepthY;
    for (i = 0; i < MB_BLOCK_SIZE; i++) {
        for (j = 0; j < MB_BLOCK_SIZE; j++) {
            bits_read += GetBits(buf, bitoffset, &val, BitstreamLengthInBits, bitdepth);
            currSlice->cof[0][i][j] = val;

            bitoffset += bitdepth;
        }
    }

    // read chroma values
    bitdepth = sps->BitDepthC;
    if (dec_picture->chroma_format_idc != YUV400 && !sps->separate_colour_plane_flag) {
        for (uv = 1; uv < 3; ++uv) {
            for (i = 0; i < sps->MbHeightC; ++i) {
                for (j = 0; j < sps->MbWidthC; ++j) {
                    bits_read += GetBits(buf, bitoffset, &val, BitstreamLengthInBits, bitdepth);
                    currSlice->cof[uv][i][j] = val;

                    bitoffset += bitdepth;
                }
            }
        }
    }

    (*dep->Dcodestrm_len) += ( bits_read >> 3);
    if (bits_read & 7)
        ++(*dep->Dcodestrm_len);

    arideco_start_decoding(&currStream->de_cabac, currStream->streamBuffer, currStream->read_len, &currStream->read_len);
}

static void read_IPCM_coeffs_from_NAL(slice_t *currSlice, struct datapartition_dec *dP)
{
    sps_t *sps = currSlice->active_sps;
    pps_t *pps = currSlice->active_pps;

    Bitstream *currStream = dP->bitstream;
    int i, j;

    //For CABAC, we don't need to read bits to let stream byte aligned
    //  because we have variable for integer bytes position
    if (pps->entropy_coding_mode_flag) {
        readIPCM_CABAC(currSlice, dP);
        return;
    }

    //read bits to let stream byte aligned
    if ((dP->bitstream->frame_bitoffset & 0x07) != 0)
        currStream->f(8 - (currStream->frame_bitoffset & 0x07));

    //read luma and chroma IPCM coefficients
    for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++)
            currSlice->cof[0][i][j] = currStream->f(sps->BitDepthY);
    }

    if (sps->chroma_format_idc != YUV400 && !sps->separate_colour_plane_flag) {
        for (i = 0; i < sps->MbHeightC; i++) {
            for (j = 0; j < sps->MbWidthC; j++)
                currSlice->cof[1][i][j] = currStream->f(sps->BitDepthC);
        }
        for (i = 0; i < sps->MbHeightC; i++) {
            for (j = 0; j < sps->MbWidthC; j++)
                currSlice->cof[2][i][j] = currStream->f(sps->BitDepthC);
        }
    }
}


static inline void reset_coeffs(mb_t *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    if (!p_Vid->active_pps->entropy_coding_mode_flag)
        memset(p_Vid->nz_coeff[currMB->mbAddrX][0][0], 0, 3 * BLOCK_PIXELS * sizeof(byte));
}

static inline void field_flag_inference(mb_t *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    if (currMB->mbAvailA)
        currMB->mb_field_decoding_flag = p_Vid->mb_data[currMB->mbAddrA].mb_field_decoding_flag;
    else
        // check top macroblock pair
        currMB->mb_field_decoding_flag = currMB->mbAvailB ? p_Vid->mb_data[currMB->mbAddrB].mb_field_decoding_flag : FALSE;
}


static void skip_macroblock(mb_t *currMB)
{
    MotionVector pred_mv;
    int zeroMotionAbove;
    int zeroMotionLeft;
    PixelPos mb[4];    // neighbor blocks
    int   i, j;
    int   a_mv_y = 0;
    int   a_ref_idx = 0;
    int   b_mv_y = 0;
    int   b_ref_idx = 0;
    int   img_block_y   = currMB->block_y;
    VideoParameters *p_Vid = currMB->p_Vid;
    slice_t *currSlice = currMB->p_Slice;
    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    StorablePicture *dec_picture = currSlice->dec_picture;
    MotionVector *a_mv = NULL;
    MotionVector *b_mv = NULL;

    get_neighbors(currMB, mb, 0, 0, MB_BLOCK_SIZE);
    if (!currSlice->MbaffFrameFlag) {
        if (mb[0].available) {
            a_mv      = &dec_picture->mv_info[mb[0].pos_y][mb[0].pos_x].mv[LIST_0];
            a_mv_y    = a_mv->mv_y;    
            a_ref_idx = dec_picture->mv_info[mb[0].pos_y][mb[0].pos_x].ref_idx[LIST_0];
        }
        if (mb[1].available) {
            b_mv      = &dec_picture->mv_info[mb[1].pos_y][mb[1].pos_x].mv[LIST_0];
            b_mv_y    = b_mv->mv_y;
            b_ref_idx = dec_picture->mv_info[mb[1].pos_y][mb[1].pos_x].ref_idx[LIST_0];
        }
    } else {
        if (mb[0].available) {
            a_mv      = &dec_picture->mv_info[mb[0].pos_y][mb[0].pos_x].mv[LIST_0];
            a_mv_y    = a_mv->mv_y;    
            a_ref_idx = dec_picture->mv_info[mb[0].pos_y][mb[0].pos_x].ref_idx[LIST_0];

            if (currMB->mb_field_decoding_flag && !p_Vid->mb_data[mb[0].mb_addr].mb_field_decoding_flag) {
                a_mv_y    /=2;
                a_ref_idx *=2;
            }
            if (!currMB->mb_field_decoding_flag && p_Vid->mb_data[mb[0].mb_addr].mb_field_decoding_flag) {
                a_mv_y    *=2;
                a_ref_idx >>=1;
            }
        }

        if (mb[1].available) {
            b_mv      = &dec_picture->mv_info[mb[1].pos_y][mb[1].pos_x].mv[LIST_0];
            b_mv_y    = b_mv->mv_y;
            b_ref_idx = dec_picture->mv_info[mb[1].pos_y][mb[1].pos_x].ref_idx[LIST_0];

            if (currMB->mb_field_decoding_flag && !p_Vid->mb_data[mb[1].mb_addr].mb_field_decoding_flag) {
                b_mv_y    /=2;
                b_ref_idx *=2;
            }
            if (!currMB->mb_field_decoding_flag && p_Vid->mb_data[mb[1].mb_addr].mb_field_decoding_flag) {
                b_mv_y    *=2;
                b_ref_idx >>=1;
            }
        }
    }

    zeroMotionLeft  = !mb[0].available ? 1 : a_ref_idx==0 && a_mv->mv_x == 0 && a_mv_y==0 ? 1 : 0;
    zeroMotionAbove = !mb[1].available ? 1 : b_ref_idx==0 && b_mv->mv_x == 0 && b_mv_y==0 ? 1 : 0;

    currMB->cbp = 0;
    reset_coeffs(currMB);

    if (zeroMotionAbove || zeroMotionLeft) {
        PicMotionParams **dec_mv_info = &dec_picture->mv_info[img_block_y];
        StorablePicture *cur_pic = currSlice->listX[list_offset][0];
        PicMotionParams *mv_info = NULL;
    
        for (j = 0; j < BLOCK_SIZE; ++j) {
            for (i = currMB->block_x; i < currMB->block_x + BLOCK_SIZE; ++i) {
                mv_info = &dec_mv_info[j][i];
                mv_info->ref_pic[LIST_0] = cur_pic;
                mv_info->mv     [LIST_0] = zero_mv;
                mv_info->ref_idx[LIST_0] = 0;
            }
        }
    } else {
        PicMotionParams **dec_mv_info = &dec_picture->mv_info[img_block_y];
        PicMotionParams *mv_info = NULL;
        StorablePicture *cur_pic = currSlice->listX[list_offset][0];
        GetMVPredictor(currMB, mb, &pred_mv, 0, dec_picture->mv_info, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

        // Set first block line (position img_block_y)
        for (j = 0; j < BLOCK_SIZE; ++j) {
            for (i = currMB->block_x; i < currMB->block_x + BLOCK_SIZE; ++i) {
                mv_info = &dec_mv_info[j][i];
                mv_info->ref_pic[LIST_0] = cur_pic;
                mv_info->mv     [LIST_0] = pred_mv;
                mv_info->ref_idx[LIST_0] = 0;
            }
        }
    }
}





static bool check_mb_skip_cavlc(mb_t *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int mb_nr = currMB->mbAddrX; 

    if (currSlice->MbaffFrameFlag && (mb_nr % 2) != 0)
        currMB->mb_field_decoding_flag = p_Vid->mb_data[mb_nr - 1].mb_field_decoding_flag;
    else
        currMB->mb_field_decoding_flag = 0;

    currMB->update_qp(currSlice->SliceQpY);

    //  read MB mode *****************************************************************
    DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][SE_MBTYPE]];
    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    currSE.mapping = linfo_ue;

    // VLC Non-Intra  
    if (currSlice->cod_counter == -1)
        currSlice->cod_counter = dP->bitstream->ue();

    if (currSlice->cod_counter == 0) {
        currSlice->cod_counter--;
        if (currSlice->MbaffFrameFlag) {
            int prevMbSkipped = 0;
            if ((mb_nr % 2) != 0) {
                mb_t *topMB = &p_Vid->mb_data[mb_nr-1];
                prevMbSkipped = topMB->mb_skip_flag;
            } else
                prevMbSkipped = 0;

            // read MB aff
            if ((((mb_nr % 2) == 0) || ((mb_nr % 2) && prevMbSkipped))) {
                currMB->mb_field_decoding_flag = dP->bitstream->f(1);
            }
        }
        return 1;
    }

    currSlice->cod_counter--;
    currMB->ei_flag = 0;
    currMB->mb_skip_flag = 1;      
    currMB->mb_type      = 0;

    if (currSlice->MbaffFrameFlag && (mb_nr % 2) == 0) {
        // read field flag of bottom block
        if (currSlice->cod_counter == 0) {
            currMB->mb_field_decoding_flag = dP->bitstream->f(1);
            dP->bitstream->frame_bitoffset--;
        } else if (currSlice->cod_counter > 0) {
            // check left macroblock pair first
            if (mb_is_available(mb_nr - 2, currMB) && (mb_nr % (sps->PicWidthInMbs * 2)) != 0)
                currMB->mb_field_decoding_flag = p_Vid->mb_data[mb_nr - 2].mb_field_decoding_flag;
            // check top macroblock pair
            else if (mb_is_available(mb_nr - 2 * sps->PicWidthInMbs, currMB))
                currMB->mb_field_decoding_flag = p_Vid->mb_data[mb_nr - 2 * sps->PicWidthInMbs].mb_field_decoding_flag;
            else
                currMB->mb_field_decoding_flag = FALSE;
        }
    }

    return 0;
}

static bool check_mb_skip_cabac(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;  
    VideoParameters *p_Vid = currMB->p_Vid;
    int mb_nr = currMB->mbAddrX;

    if (currSlice->MbaffFrameFlag && (mb_nr % 2) != 0)
        currMB->mb_field_decoding_flag = p_Vid->mb_data[mb_nr - 1].mb_field_decoding_flag;
    else
        currMB->mb_field_decoding_flag = 0;

    currMB->update_qp(currSlice->SliceQpY);

    //  read MB mode *****************************************************************
    DataPartition *dP = &currSlice->partArr[assignSE2partition[currSlice->dp_mode][SE_MBTYPE]];
    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    if (dP->bitstream->ei_flag)   
        currSE.mapping = linfo_ue;

    if (currSlice->MbaffFrameFlag) {
        // read MB skip_flag
        int prevMbSkipped = 0;
        if ((mb_nr % 2) != 0) {
            mb_t *topMB = &p_Vid->mb_data[mb_nr - 1];
            prevMbSkipped = topMB->mb_skip_flag;
        } else
            prevMbSkipped = 0;
        if ((mb_nr % 2) == 0 || prevMbSkipped)
            field_flag_inference(currMB);
    }

    CheckAvailabilityOfNeighborsCABAC(currMB);
    currSE.reading = read_skip_flag_CABAC;
    dP->readSyntaxElement(currMB, &currSE, dP);

    currMB->mb_type      =  currSE.value1;
    currMB->mb_skip_flag = !currSE.value1;
    currMB->cbp          =  currSE.value1;
    if (!dP->bitstream->ei_flag)
        currMB->ei_flag = 0;

    if (currSE.value1 == 0)
        currSlice->cod_counter = 0;

    if (currSlice->MbaffFrameFlag) {
        // read MB AFF
        int check_bottom, read_bottom, read_top;  
        check_bottom = read_bottom = read_top = 0;
        if ((mb_nr % 2) == 0) {
            check_bottom = currMB->mb_skip_flag;
            read_top = !check_bottom;
        } else {
            mb_t *topMB = &p_Vid->mb_data[mb_nr - 1];
            read_bottom = topMB->mb_skip_flag && !currMB->mb_skip_flag;
        }

        if (read_bottom || read_top) {
            currSE.reading = readFieldModeInfo_CABAC;
            dP->readSyntaxElement(currMB, &currSE, dP);
            currMB->mb_field_decoding_flag = currSE.value1;
        }

        if (check_bottom)
            check_next_mb_and_get_field_mode_CABAC(currSlice, &currSE, dP);
        CheckAvailabilityOfNeighborsCABAC(currMB);    
    }

    return currMB->mb_type != 0;
}


void macroblock_t::parse()
{
    slice_t *slice = this->p_Slice;

    switch (slice->slice_type) {
    case I_SLICE:
    case SI_SLICE:
        this->parse_i_slice();
        return;
    case P_SLICE:
    case SP_SLICE:
    case B_SLICE:
        this->parse_pb_slice();
        return;
    }
}

void macroblock_t::parse_i_slice()
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;
    int mb_nr = this->mbAddrX; 

    if ((mb_nr % 2) != 0)
        this->mb_field_decoding_flag = slice->mb_data[mb_nr - 1].mb_field_decoding_flag;
    else
        this->mb_field_decoding_flag = 0;

    this->update_qp(slice->SliceQpY);

    //  read MB mode *****************************************************************
    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_MBTYPE]];
    SyntaxElement currSE;
    currSE.type = SE_MBTYPE;
    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)   
        currSE.mapping = linfo_ue;

    // read MB aff
    if (slice->MbaffFrameFlag && (mb_nr & 0x01) == 0) {
        if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
            currSE.value1 = dP->bitstream->f(1);
        } else {
            currSE.reading = readFieldModeInfo_CABAC;
            dP->readSyntaxElement(this, &currSE, dP);
        }
        this->mb_field_decoding_flag = currSE.value1;
    }

    if (pps->entropy_coding_mode_flag)
        CheckAvailabilityOfNeighborsCABAC(this);

    //  read MB type
    //if (pps->entropy_coding_mode_flag)
    //    currSE.reading = readMB_typeInfo_CABAC_i_slice;
    //dP->readSyntaxElement(this, &currSE, dP);
    currSE.value1 = getSE(this, SE_MBTYPE);
    this->mb_type = currSE.value1;
    //if (!dP->bitstream->ei_flag)
    //    this->ei_flag = 0;

    PicMotionParamsOld *motion = &slice->dec_picture->motion;
    motion->mb_field_decoding_flag[mb_nr] = this->mb_field_decoding_flag;

    this->interpret_mb_mode();

    //init NoMbPartLessThan8x8Flag
    this->NoMbPartLessThan8x8Flag = TRUE;

    if (this->mb_type == IPCM)
        this->parse_i_pcm();
    else
        this->parse_intra();
}

void macroblock_t::parse_pb_slice()
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;
    int mb_nr = this->mbAddrX;

    bool mb_skip;
    if (pps->entropy_coding_mode_flag)
        mb_skip = check_mb_skip_cabac(this);
    else
        mb_skip = check_mb_skip_cavlc(this);

    // read MB type
    if (mb_skip) {
        //DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_MBTYPE]];
        //currSE.type = SE_MBTYPE;
        //if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)   
        //    currSE.mapping = linfo_ue;
        //else
        //    currSE.reading = slice->slice_type != B_SLICE ? readMB_typeInfo_CABAC_p_slice :
        //                                                        readMB_typeInfo_CABAC_b_slice;
        //dP->readSyntaxElement(this, &currSE, dP);
        SyntaxElement currSE;
        currSE.value1 = getSE(this, SE_MBTYPE);
        if (!pps->entropy_coding_mode_flag && slice->slice_type != B_SLICE)
            ++(currSE.value1);
        this->mb_type = (short) currSE.value1;
        //if (!dP->bitstream->ei_flag)
        //    this->ei_flag = 0;
        if (!pps->entropy_coding_mode_flag)
            this->mb_skip_flag = 0;
    }

    PicMotionParamsOld *motion = &slice->dec_picture->motion;
    if (!slice->MbaffFrameFlag)
        motion->mb_field_decoding_flag[mb_nr] = 0;
    else
        motion->mb_field_decoding_flag[mb_nr] = this->mb_field_decoding_flag;

    this->interpret_mb_mode();    

    if (slice->MbaffFrameFlag && this->mb_field_decoding_flag) {
        slice->num_ref_idx_l0_active_minus1 = ((slice->num_ref_idx_l0_active_minus1 + 1) << 1) - 1;
        slice->num_ref_idx_l1_active_minus1 = ((slice->num_ref_idx_l1_active_minus1 + 1) << 1) - 1;
    }

    if (slice->slice_type != B_SLICE)
        //init NoMbPartLessThan8x8Flag
        this->NoMbPartLessThan8x8Flag = TRUE;

    if (this->mb_type == IPCM) // I_PCM mode
        this->parse_i_pcm();
    else if (this->mb_type == PSKIP || this->mb_type == BSKIP_DIRECT)
        this->parse_skip();
    else if (this->is_intra_block) // all other intra modes
        this->parse_intra();
    else // all other remaining modes
        this->parse_inter();
}

void macroblock_t::parse_skip()
{
    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;
    pps_t *pps = slice->active_pps;

    if (slice->slice_type == B_slice) {
        //init NoMbPartLessThan8x8Flag
        this->NoMbPartLessThan8x8Flag = sps->direct_8x8_inference_flag;
        this->transform_size_8x8_flag = 0;
        if (pps->constrained_intra_pred_flag)
            slice->intra_block[this->mbAddrX] = 0;

        //--- init macroblock data ---
        init_macroblock_direct(this);

        if (slice->cod_counter >= 0) {
            this->cbp = 0;
            if (pps->entropy_coding_mode_flag) {
                slice->is_reset_coeff = TRUE;
                slice->cod_counter = -1;
            } else
                reset_coeffs(this);
        } else {
            // read CBP and Coeffs  ***************************************************************
            if (pps->entropy_coding_mode_flag)
                this->read_CBP_and_coeffs_from_NAL_CABAC();
            else
                this->read_CBP_and_coeffs_from_NAL_CAVLC();
        }
    } else {
        this->transform_size_8x8_flag = 0;
        if (pps->constrained_intra_pred_flag)
            slice->intra_block[this->mbAddrX] = 0;

        //--- init macroblock data ---
        init_macroblock_basic(this);
        skip_macroblock(this);
    }
}

void macroblock_t::parse_intra()
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;

    if (this->mb_type != I4MB)
        this->NoMbPartLessThan8x8Flag = 1;

    this->transform_size_8x8_flag = 0;
    if (this->mb_type == I4MB && pps->transform_8x8_mode_flag) {
        DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_HEADER]];
        SyntaxElement currSE;
        currSE.type = SE_HEADER;
        if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
            currSE.value1 = dP->bitstream->f(1);
        } else {
            currSE.reading = readMB_transform_size_flag_CABAC;
            dP->readSyntaxElement(this, &currSE, dP);
        }

        this->transform_size_8x8_flag = currSE.value1;
        if (this->transform_size_8x8_flag) {
            this->mb_type = I8MB;
            for (int i = 0; i < 4; i++) {
                this->b8mode[i] = I8MB;
                this->b8pdir[i] = -1;
            }
        }
    }

    init_macroblock(this);
    this->parse_ipred_modes();

    if (pps->entropy_coding_mode_flag)
        this->read_CBP_and_coeffs_from_NAL_CABAC();
    else
        this->read_CBP_and_coeffs_from_NAL_CAVLC();
}

void macroblock_t::parse_inter()
{
    static const char p_v2b8 [ 5] = { 4, 5, 6, 7, IBLOCK};
    static const char p_v2pd [ 5] = { 0, 0, 0, 0, -1};
    static const char b_v2b8 [14] = { 0, 4, 4, 4, 5, 6, 5, 6, 5, 6, 7, 7, 7, IBLOCK};
    static const char b_v2pd [14] = { 2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1, 2, -1};

    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;
    pps_t *pps = slice->active_pps;

    //init NoMbPartLessThan8x8Flag
    this->NoMbPartLessThan8x8Flag = 1;
    this->transform_size_8x8_flag = 0;
    if (this->mb_type == P8x8) {
        DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_MBTYPE]];
        SyntaxElement currSE;
        currSE.type = SE_MBTYPE;      
        if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) 
            currSE.mapping = linfo_ue;
        else
            currSE.reading = readB8_typeInfo_CABAC;

        for (int mbPartIdx = 0; mbPartIdx < 4; mbPartIdx++) {
            dP->readSyntaxElement(this, &currSE, dP);
            this->sub_mb_type[mbPartIdx] = currSE.value1;

            const uint8_t (*sub_mb_types)[5] = slice->slice_type == B_slice ?
                sub_mb_types_b_slice : sub_mb_types_p_slice;

            this->NumSubMbPart   [mbPartIdx] = sub_mb_types[this->sub_mb_type[mbPartIdx]][1];
            this->SubMbPredMode  [mbPartIdx] = sub_mb_types[this->sub_mb_type[mbPartIdx]][2];
            this->SubMbPartWidth [mbPartIdx] = sub_mb_types[this->sub_mb_type[mbPartIdx]][3];
            this->SubMbPartHeight[mbPartIdx] = sub_mb_types[this->sub_mb_type[mbPartIdx]][4];

            if (slice->slice_type == B_slice) {
                this->b8mode[mbPartIdx] = b_v2b8[currSE.value1];
                this->b8pdir[mbPartIdx] = b_v2pd[currSE.value1];
            } else {
                this->b8mode[mbPartIdx] = p_v2b8[currSE.value1];
                this->b8pdir[mbPartIdx] = p_v2pd[currSE.value1];
            }

            //set NoMbPartLessThan8x8Flag for P8x8 mode
            this->NoMbPartLessThan8x8Flag &= 
                (this->b8mode[mbPartIdx] == 0 && sps->direct_8x8_inference_flag) ||
                (this->b8mode[mbPartIdx] == 4);
        }
    }

    if (pps->constrained_intra_pred_flag)
        slice->intra_block[this->mbAddrX] = 0;

    init_macroblock(this);
    this->parse_motion_info();

    if (pps->entropy_coding_mode_flag)
        this->read_CBP_and_coeffs_from_NAL_CABAC();
    else
        this->read_CBP_and_coeffs_from_NAL_CAVLC();
}

void macroblock_t::parse_i_pcm()
{
    slice_t *slice = this->p_Slice;

    this->NoMbPartLessThan8x8Flag = 1;
    this->transform_size_8x8_flag = 0;

    //--- init macroblock data ---
    init_macroblock(this);

    //read pcm_alignment_zero_bit and pcm_byte[i]

    // here dP is assigned with the same dP as SE_MBTYPE, because IPCM syntax is in the
    // same category as MBTYPE
    if (slice->dp_mode && slice->dpB_NotPresent)
        concealIPCMcoeffs(this);
    else {
        DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_LUM_DC_INTRA]];
        read_IPCM_coeffs_from_NAL(slice, dP);
    }
}


void macroblock_t::parse_ipred_modes()
{
    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;
    pps_t *pps = slice->active_pps;

    if (this->mb_type == I8MB)
        this->parse_ipred_8x8_modes();
    else if (this->mb_type == I4MB)
        this->parse_ipred_4x4_modes();

    if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
        DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_INTRAPREDMODE]];
        SyntaxElement currSE;
        currSE.type = SE_INTRAPREDMODE;
        if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
            currSE.mapping = linfo_ue;
        else
            currSE.reading = readCIPredMode_CABAC;

        dP->readSyntaxElement(this, &currSE, dP);
        this->intra_chroma_pred_mode = currSE.value1;

        if (this->intra_chroma_pred_mode < Intra_Chroma_DC || this->intra_chroma_pred_mode > Intra_Chroma_Plane)
            error("illegal chroma intra pred mode!\n", 600);
    }
}

void macroblock_t::parse_ipred_4x4_modes()
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;

    for (int luma4x4BlkIdx = 0; luma4x4BlkIdx < 16; luma4x4BlkIdx++) {
        int bx = ((luma4x4BlkIdx / 4) % 2) * 8 + ((luma4x4BlkIdx % 4) % 2) * 4;
        int by = ((luma4x4BlkIdx / 4) / 2) * 8 + ((luma4x4BlkIdx % 4) / 2) * 4;

        DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_INTRAPREDMODE]];
        SyntaxElement currSE;
        currSE.type = SE_INTRAPREDMODE;
        if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
            if (dP->bitstream->f(1))
                currSE.value1 = -1;
            else
                currSE.value1 = dP->bitstream->f(3);
        } else {
            currSE.context = luma4x4BlkIdx;
            currSE.reading = readIntraPredMode_CABAC;
            dP->readSyntaxElement(this, &currSE, dP);
        }

        this->prev_intra4x4_pred_mode_flag[luma4x4BlkIdx] = currSE.value1 == -1;
        this->rem_intra4x4_pred_mode      [luma4x4BlkIdx] = currSE.value1;

        int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
        PixelPos left_block, top_block;
        get4x4Neighbour(this, bx - 1, by    , mb_size, &left_block);
        get4x4Neighbour(this, bx    , by - 1, mb_size, &top_block );

        //get from array and decode
        if (pps->constrained_intra_pred_flag) {
            left_block.available &= slice->p_Vid->mb_data[left_block.mb_addr].is_intra_block;
            top_block.available  &= slice->p_Vid->mb_data[top_block.mb_addr ].is_intra_block;
        }
        // !! KS: not sure if the following is still correct...
        if (slice->slice_type == SI_slice) { // need support for MBINTLC1
            left_block.available &= slice->p_Vid->mb_data[left_block.mb_addr].mb_type == SI4MB;
            top_block.available  &= slice->p_Vid->mb_data[top_block.mb_addr ].mb_type == SI4MB;
        }

        bool dcPredModePredictedFlag = !left_block.available || !top_block.available;

        int scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };
        uint8_t intraMxMPredModeA = Intra_4x4_DC;
        uint8_t intraMxMPredModeB = Intra_4x4_DC;
        if (!dcPredModePredictedFlag) {
            uint8_t left_mb_type = slice->p_Vid->mb_data[left_block.mb_addr].mb_type;
            uint8_t top_mb_type  = slice->p_Vid->mb_data[top_block.mb_addr ].mb_type;
            if (left_mb_type == I8MB)
                intraMxMPredModeA = slice->p_Vid->mb_data[left_block.mb_addr].Intra8x8PredMode[scan[left_block.y * 4 + left_block.x] / 4];
            else if (left_mb_type == I4MB)
                intraMxMPredModeA = slice->p_Vid->mb_data[left_block.mb_addr].Intra4x4PredMode[scan[left_block.y * 4 + left_block.x]];
            if (top_mb_type == I8MB)
                intraMxMPredModeB = slice->p_Vid->mb_data[top_block.mb_addr].Intra8x8PredMode[scan[top_block.y * 4 + top_block.x] / 4];
            else if (top_mb_type == I4MB)
                intraMxMPredModeB = slice->p_Vid->mb_data[top_block.mb_addr].Intra4x4PredMode[scan[top_block.y * 4 + top_block.x]];
        }

        uint8_t predIntra4x4PredMode = imin(intraMxMPredModeA, intraMxMPredModeB);
        if (this->prev_intra4x4_pred_mode_flag[luma4x4BlkIdx])
            this->Intra4x4PredMode[luma4x4BlkIdx] = predIntra4x4PredMode;
        else if (this->rem_intra4x4_pred_mode[luma4x4BlkIdx] < predIntra4x4PredMode)
            this->Intra4x4PredMode[luma4x4BlkIdx] = this->rem_intra4x4_pred_mode[luma4x4BlkIdx];
        else
            this->Intra4x4PredMode[luma4x4BlkIdx] = this->rem_intra4x4_pred_mode[luma4x4BlkIdx] + 1;

        //if ((luma4x4BlkIdx % 4) == 0)
        //    this->Intra8x8PredMode[luma4x4BlkIdx / 4] = this->Intra4x4PredMode[luma4x4BlkIdx];
    }
}

void macroblock_t::parse_ipred_8x8_modes()
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;

    for (int luma8x8BlkIdx = 0; luma8x8BlkIdx < 4; luma8x8BlkIdx++) {
        int bx = (luma8x8BlkIdx % 2) * 8;
        int by = (luma8x8BlkIdx / 2) * 8;

        DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_INTRAPREDMODE]];
        SyntaxElement currSE;
        currSE.type = SE_INTRAPREDMODE;
        if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) {
            if (dP->bitstream->f(1))
                currSE.value1 = -1;
            else
                currSE.value1 = dP->bitstream->f(3);
        } else {
            currSE.context = luma8x8BlkIdx * 4;
            currSE.reading = readIntraPredMode_CABAC;
            dP->readSyntaxElement(this, &currSE, dP);
        }

        this->prev_intra8x8_pred_mode_flag[luma8x8BlkIdx] = currSE.value1 == -1;
        this->rem_intra8x8_pred_mode      [luma8x8BlkIdx] = currSE.value1;

        int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
        PixelPos left_block, top_block;
        get4x4Neighbour(this, bx - 1, by    , mb_size, &left_block);
        get4x4Neighbour(this, bx    , by - 1, mb_size, &top_block);

        //get from array and decode
        if (pps->constrained_intra_pred_flag) {
            left_block.available &= slice->p_Vid->mb_data[left_block.mb_addr].is_intra_block;
            top_block.available  &= slice->p_Vid->mb_data[top_block.mb_addr ].is_intra_block;
        }

        bool dcPredModePredictedFlag = !left_block.available || !top_block.available;

        int scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };
        uint8_t intraMxMPredModeA = Intra_8x8_DC;
        uint8_t intraMxMPredModeB = Intra_8x8_DC;
        if (!dcPredModePredictedFlag) {
            uint8_t left_mb_type = slice->p_Vid->mb_data[left_block.mb_addr].mb_type;
            uint8_t top_mb_type  = slice->p_Vid->mb_data[top_block.mb_addr ].mb_type;
            if (left_mb_type == I8MB)
                intraMxMPredModeA = slice->p_Vid->mb_data[left_block.mb_addr].Intra8x8PredMode[scan[left_block.y * 4 + left_block.x] / 4];
            else if (left_mb_type == I4MB)
                intraMxMPredModeA = slice->p_Vid->mb_data[left_block.mb_addr].Intra4x4PredMode[scan[left_block.y * 4 + left_block.x]];
            if (top_mb_type == I8MB)
                intraMxMPredModeB = slice->p_Vid->mb_data[top_block.mb_addr].Intra8x8PredMode[scan[top_block.y * 4 + top_block.x] / 4];
            else if (top_mb_type == I4MB)
                intraMxMPredModeB = slice->p_Vid->mb_data[top_block.mb_addr].Intra4x4PredMode[scan[top_block.y * 4 + top_block.x]];
        }

        uint8_t predIntra8x8PredMode = imin(intraMxMPredModeA, intraMxMPredModeB);
        if (this->prev_intra8x8_pred_mode_flag[luma8x8BlkIdx])
            this->Intra8x8PredMode[luma8x8BlkIdx] = predIntra8x8PredMode;
        else if (this->rem_intra8x8_pred_mode[luma8x8BlkIdx] < predIntra8x8PredMode)
            this->Intra8x8PredMode[luma8x8BlkIdx] = this->rem_intra8x8_pred_mode[luma8x8BlkIdx];
        else
            this->Intra8x8PredMode[luma8x8BlkIdx] = this->rem_intra8x8_pred_mode[luma8x8BlkIdx] + 1;

        //this->Intra4x4PredMode[luma8x8BlkIdx * 4    ] = this->Intra8x8PredMode[luma8x8BlkIdx];
        //this->Intra4x4PredMode[luma8x8BlkIdx * 4 + 1] = this->Intra8x8PredMode[luma8x8BlkIdx];
        //this->Intra4x4PredMode[luma8x8BlkIdx * 4 + 2] = this->Intra8x8PredMode[luma8x8BlkIdx];
        //this->Intra4x4PredMode[luma8x8BlkIdx * 4 + 3] = this->Intra8x8PredMode[luma8x8BlkIdx];
    }
}


void macroblock_t::parse_motion_info()
{
    slice_t *slice = this->p_Slice;

    if (slice->slice_type == B_slice && this->mb_type == P8x8)
        update_direct_mv_info(this);   

    this->parse_ref_pic_idx(LIST_0);
    if (slice->slice_type == B_slice)
        this->parse_ref_pic_idx(LIST_1);

    this->parse_motion_vectors(LIST_0);
    if (slice->slice_type == B_slice)
        this->parse_motion_vectors(LIST_1);

    int list_offset = slice->MbaffFrameFlag && this->mb_field_decoding_flag ?
                      this->mbAddrX % 2 ? 4 : 2 : 0;
    StorablePicture **list0 = slice->listX[LIST_0 + list_offset];
    StorablePicture **list1 = slice->listX[LIST_1 + list_offset];
    PicMotionParams **p_mv_info = &slice->dec_picture->mv_info[this->block_y];
    // record reference picture Ids for deblocking decisions
    PicMotionParams *mv_info;
    short ref_idx;
    for (int j4 = 0; j4 < 4; j4++) {
        for (int i4 = this->block_x; i4 < this->block_x + 4; ++i4) {
            mv_info = &p_mv_info[j4][i4];
            ref_idx = mv_info->ref_idx[LIST_0];
            mv_info->ref_pic[LIST_0] = (ref_idx >= 0) ? list0[ref_idx] : NULL;
            if (slice->slice_type == B_slice) {
                ref_idx = mv_info->ref_idx[LIST_1];
                mv_info->ref_pic[LIST_1] = (ref_idx >= 0) ? list1[ref_idx] : NULL;
            }
        }
    }
}

//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

void macroblock_t::parse_ref_pic_idx(int list)
{
    slice_t *slice = this->p_Slice;
    PicMotionParams **mv_info = slice->dec_picture->mv_info;

    int partmode = this->mb_type == P8x8 ? 4 : this->mb_type;
    int step_h0  = BLOCK_STEP [partmode][0];
    int step_v0  = BLOCK_STEP [partmode][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {      
            int k = 2 * (j0 >> 1) + (i0 >> 1);

            if ((this->b8pdir[k] == list || this->b8pdir[k] == BI_PRED) &&
                this->b8mode[k] != 0) {
                this->subblock_x = i0 << 2;
                this->subblock_y = j0 << 2;
                char refframe = readRefPictureIdx(this, this->b8mode[k], list);
                for (int j = 0; j < step_v0; j++) {
                    for (int i = 0; i < step_h0; i++)
                        mv_info[this->block_y + j0 + j][this->block_x + i0 + i].ref_idx[list] = refframe;
                }
            }
        }
    }
}

void macroblock_t::parse_motion_vectors(int list)
{
    slice_t *slice = this->p_Slice;
    PicMotionParams **mv_info = slice->dec_picture->mv_info;

    int partmode = (this->mb_type == P8x8 ? 4 : this->mb_type);
    int step_h0  = BLOCK_STEP [partmode][0];
    int step_v0  = BLOCK_STEP [partmode][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int kk = 2 * (j0 >> 1) + (i0 >> 1);
            int step_h4 = BLOCK_STEP[(int)this->b8mode[kk]][0];
            int step_v4 = BLOCK_STEP[(int)this->b8mode[kk]][1];
            char cur_ref_idx = mv_info[this->block_y + j0][this->block_x + i0].ref_idx[list];

            if ((this->b8pdir[kk] == list || this->b8pdir[kk] == BI_PRED) &&
                (this->b8mode[kk] != 0)) { //has forward vector
                for (int j = j0; j < j0 + step_v0; j += step_v4) {
                    for (int i = i0; i < i0 + step_h0; i += step_h4)
                        this->parse_motion_vector(list, step_h4 * 4, step_v4 * 4, i, j, cur_ref_idx);
                }
            }
        }
    }
}

void macroblock_t::parse_motion_vector(int list, int step_h4, int step_v4, int i, int j, char cur_ref_idx)
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;

    PicMotionParams **mv_info = slice->dec_picture->mv_info;
    //const uint8_t (*sub_mb_types)[5] = this->p_Slice->slice_type == B_slice ?
    //                                   sub_mb_types_b_slice : sub_mb_types_p_slice;
    //int step_h4 = sub_mb_types[this->sub_mb_type[kk]][3];
    //int step_v4 = sub_mb_types[this->sub_mb_type[kk]][4];

    this->subblock_x = i * 4;
    this->subblock_y = j * 4;

    PixelPos block[4]; // neighbor blocks
    MotionVector pred_mv;
    get_neighbors(this, block, BLOCK_SIZE * i, BLOCK_SIZE * j, step_h4);
    GetMVPredictor(this, block, &pred_mv, cur_ref_idx, mv_info, list, BLOCK_SIZE * i, BLOCK_SIZE * j, step_h4, step_v4);

    DataPartition *dP = &slice->partArr[assignSE2partition[slice->dp_mode][SE_MVD]];
    SyntaxElement currSE;
    currSE.type = SE_MVD;
    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag) 
        currSE.mapping = linfo_se;
    else
        currSE.reading = slice->MbaffFrameFlag ? read_mvd_CABAC_mbaff : read_MVD_CABAC;

    short curr_mvd[2];
    for (int k = 0; k < 2; ++k) {
        currSE.value2 = k * 2 + list;
        dP->readSyntaxElement(this, &currSE, dP);
        curr_mvd[k] = currSE.value1;              
    }

    MotionVector curr_mv;
    curr_mv.mv_x = curr_mvd[0] + pred_mv.mv_x;
    curr_mv.mv_y = curr_mvd[1] + pred_mv.mv_y;

    int i4 = this->block_x + i;
    int j4 = this->block_y + j;
    for (int jj = 0; jj < step_v4 / 4; ++jj) {
        for (int ii = 0; ii < step_h4 / 4; ++ii) {
            mv_info[jj + j4][ii + i4].mv[list] = curr_mv;
            this->mvd[list][jj + j][ii + i][0] = curr_mvd[0];
            this->mvd[list][jj + j][ii + i][1] = curr_mvd[1];
        }
    }
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

void linfo_cbp_intra_normal(int len, int info, int *cbp, int *dummy)
{
    int cbp_idx;
    linfo_ue(len, info, &cbp_idx, dummy);
    *cbp = NCBP[1][cbp_idx][0];
}

void linfo_cbp_intra_other(int len, int info, int *cbp, int *dummy)
{
    int cbp_idx;
    linfo_ue(len, info, &cbp_idx, dummy);
    *cbp = NCBP[0][cbp_idx][0];
}

void linfo_cbp_inter_normal(int len, int info, int *cbp, int *dummy)
{
    int cbp_idx;
    linfo_ue(len, info, &cbp_idx, dummy);
    *cbp = NCBP[1][cbp_idx][1];
}

void linfo_cbp_inter_other(int len, int info, int *cbp, int *dummy)
{
    int cbp_idx;
    linfo_ue(len, info, &cbp_idx, dummy);
    *cbp = NCBP[0][cbp_idx][1];
}

void macroblock_t::read_delta_quant(SyntaxElement *currSE, DataPartition *dP, const byte *partMap, int type)
{
    mb_t *currMB = this;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    pps_t *pps = currSlice->active_pps;
 
    currSE->type = type;
    dP = &currSlice->partArr[partMap[currSE->type]];

    if (!pps->entropy_coding_mode_flag || dP->bitstream->ei_flag)
        currSE->mapping = linfo_se;
    else
        currSE->reading = read_dQuant_CABAC;

    dP->readSyntaxElement(currMB, currSE, dP);
    currMB->delta_quant = (short) currSE->value1;
    if (currMB->delta_quant < -(26 + sps->QpBdOffsetY / 2) ||
        currMB->delta_quant >  (25 + sps->QpBdOffsetY / 2)) {
        printf("mb_qp_delta is out of range (%d)\n", currMB->delta_quant);
        currMB->delta_quant = iClip3(-(26 + sps->QpBdOffsetY/2), (25 + sps->QpBdOffsetY/2), currMB->delta_quant);
    }

    currSlice->SliceQpY = ((currSlice->SliceQpY + currMB->delta_quant + 52 + 2*sps->QpBdOffsetY)%(52+sps->QpBdOffsetY)) - sps->QpBdOffsetY;
    currMB->update_qp(currSlice->SliceQpY);
}
