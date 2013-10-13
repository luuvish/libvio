/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : parser_mb.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "parser.h"
#include "intra_prediction.h"
#include "transform.h"


namespace vio  {
namespace h264 {


// Table 7-11 mb_t types for I slices
static const uint8_t mb_types_i_slice[26][5] = {
    { I_NxN  , Intra_NxN  , NA,  0,  0 },
    { I_16x16, Intra_16x16,  0,  0,  0 },
    { I_16x16, Intra_16x16,  1,  0,  0 },
    { I_16x16, Intra_16x16,  2,  0,  0 },
    { I_16x16, Intra_16x16,  3,  0,  0 },
    { I_16x16, Intra_16x16,  0,  1,  0 },
    { I_16x16, Intra_16x16,  1,  1,  0 },
    { I_16x16, Intra_16x16,  2,  1,  0 },
    { I_16x16, Intra_16x16,  3,  1,  0 },
    { I_16x16, Intra_16x16,  0,  2,  0 },
    { I_16x16, Intra_16x16,  1,  2,  0 },
    { I_16x16, Intra_16x16,  2,  2,  0 },
    { I_16x16, Intra_16x16,  3,  2,  0 },
    { I_16x16, Intra_16x16,  0,  0, 15 },
    { I_16x16, Intra_16x16,  1,  0, 15 },
    { I_16x16, Intra_16x16,  2,  0, 15 },
    { I_16x16, Intra_16x16,  3,  0, 15 },
    { I_16x16, Intra_16x16,  0,  1, 15 },
    { I_16x16, Intra_16x16,  1,  1, 15 },
    { I_16x16, Intra_16x16,  2,  1, 15 },
    { I_16x16, Intra_16x16,  3,  1, 15 },
    { I_16x16, Intra_16x16,  0,  2, 15 },
    { I_16x16, Intra_16x16,  1,  2, 15 },
    { I_16x16, Intra_16x16,  2,  2, 15 },
    { I_16x16, Intra_16x16,  3,  2, 15 },
    { I_PCM  , NA         , NA, NA, NA }
};

// Table 7-12 Macroblock type with value 0 for SI slices
static const uint8_t mb_types_si_slice[1][5] = {
    { SI     , Intra_4x4  , NA, NA, NA }
    // [1..26] are mb_types_i_slice[0..25]
};

// Table 7-13 mb_t type value 0 to 4 for P and SP slices
static const uint8_t mb_types_p_slice[5][6] = {
//  { P_Skip        ,  1, Pred_L0, NA     , 16, 16 },
    { P_16x16       ,  1, Pred_L0, NA     , 16, 16 },
    { P_16x8        ,  2, Pred_L0, Pred_L0, 16,  8 },
    { P_8x16        ,  2, Pred_L0, Pred_L0,  8, 16 },
    { P_8x8         ,  4, NA     , NA     ,  8,  8 },
    { P_8x8ref0     ,  4, NA     , NA     ,  8,  8 }
    // [5..30] are mb_types_i_slice[0..25]
};

// Table 7-14 mb_t type value 0 to 22 for B slices
static const uint8_t mb_types_b_slice[23][6] = {
//  { B_Skip        , NA, Direct , NA     ,  8,  8 },
    { B_Direct_16x16, NA, Direct , NA     ,  8,  8 },
    { B_16x16       ,  1, Pred_L0, NA     , 16, 16 },
    { B_16x16       ,  1, Pred_L1, NA     , 16, 16 },
    { B_16x16       ,  1, BiPred , NA     , 16, 16 },
    { B_16x8        ,  2, Pred_L0, Pred_L0, 16,  8 },
    { B_8x16        ,  2, Pred_L0, Pred_L0,  8, 16 },
    { B_16x8        ,  2, Pred_L1, Pred_L1, 16,  8 },
    { B_8x16        ,  2, Pred_L1, Pred_L1,  8, 16 },
    { B_16x8        ,  2, Pred_L0, Pred_L1, 16,  8 },
    { B_8x16        ,  2, Pred_L0, Pred_L1,  8, 16 },
    { B_16x8        ,  2, Pred_L1, Pred_L0, 16,  8 },
    { B_8x16        ,  2, Pred_L1, Pred_L0,  8, 16 },
    { B_16x8        ,  2, Pred_L0, BiPred , 16,  8 },
    { B_8x16        ,  2, Pred_L0, BiPred ,  8, 16 },
    { B_16x8        ,  2, Pred_L1, BiPred , 16,  8 },
    { B_8x16        ,  2, Pred_L1, BiPred ,  8, 16 },
    { B_16x8        ,  2, BiPred , Pred_L0, 16,  8 },
    { B_8x16        ,  2, BiPred , Pred_L0,  8, 16 },
    { B_16x8        ,  2, BiPred , Pred_L1, 16,  8 },
    { B_8x16        ,  2, BiPred , Pred_L1,  8, 16 },
    { B_16x8        ,  2, BiPred , BiPred , 16,  8 },
    { B_8x16        ,  2, BiPred , BiPred ,  8, 16 },
    { B_8x8         ,  4, NA     , NA     ,  8,  8 }
    // [23..48] are mb_types_i_slice[0..25]
};

// Table 7-17 Sub-macroblock types in P macroblocks
static const uint8_t sub_mb_types_p_slice[4][5] = {
    { P_8x8         ,  1, Pred_L0, 8, 8 },
    { P_8x4         ,  2, Pred_L0, 8, 4 },
    { P_4x8         ,  2, Pred_L0, 4, 8 },
    { P_4x4         ,  4, Pred_L0, 4, 4 }
};

// Table 7-18 Sub-macroblock types in B macroblocks
static const uint8_t sub_mb_types_b_slice[13][5] = {
//  {           NA  ,  4, Direct , 4, 4 },
    { B_Direct_8x8  ,  4, Direct , 4, 4 },
    { B_8x8         ,  1, Pred_L0, 8, 8 },
    { B_8x8         ,  1, Pred_L1, 8, 8 },
    { B_8x8         ,  1, BiPred , 8, 8 },
    { B_8x4         ,  2, Pred_L0, 8, 4 },
    { B_4x8         ,  2, Pred_L0, 4, 8 },
    { B_8x4         ,  2, Pred_L1, 8, 4 },
    { B_4x8         ,  2, Pred_L1, 4, 8 },
    { B_8x4         ,  2, BiPred , 8, 4 },
    { B_4x8         ,  2, BiPred , 4, 8 },
    { B_4x4         ,  4, Pred_L0, 4, 4 },
    { B_4x4         ,  4, Pred_L1, 4, 4 },
    { B_4x4         ,  4, BiPred , 4, 4 }
};


void Parser::init(slice_t& slice)
{
    shr_t& shr = slice.header;

    this->mb_skip_run = -1;

    this->prescan_skip_read = false;
    this->prescan_mb_field_decoding_read = false;

    this->QpY = shr.SliceQpY;

    this->is_reset_coeff    = false;
    this->is_reset_coeff_cr = false;

    if (slice.active_pps->entropy_coding_mode_flag) {
        this->mot_ctx.init(shr.slice_type, shr.cabac_init_idc, shr.SliceQpY);
        this->last_dquant = 0;
    }
}


void Parser::parse(mb_t& mb)
{
    Macroblock mbp { mb };
    mbp.parse();
}


Parser::Macroblock::Macroblock(mb_t& _mb) :
    sps { *_mb.p_Slice->active_sps },
    pps { *_mb.p_Slice->active_pps },
    slice { *_mb.p_Slice },
    mb { _mb },
    se { _mb },
    re { _mb }
{
}

Parser::Macroblock::~Macroblock()
{
}

void Parser::Macroblock::parse()
{
    shr_t& shr = slice.header;

    int CurrMbAddr = mb.mbAddrX;
    bool moreDataFlag = 1;

    if (shr.slice_type != I_slice && shr.slice_type != SI_slice) {
        if (pps.entropy_coding_mode_flag) {
            if (slice.parser.prescan_skip_read) {
                slice.parser.prescan_skip_read = false;
                mb.mb_skip_flag = slice.parser.prescan_skip_flag;
            } else {
                mb.mb_skip_flag = se.mb_skip_flag();
                if (mb.mb_skip_flag)
                    slice.parser.last_dquant = 0;
                mb.ei_flag = 0;
            }
        } else {
            if (slice.parser.mb_skip_run == -1)
                slice.parser.mb_skip_run = se.mb_skip_run();
            mb.mb_skip_flag = (slice.parser.mb_skip_run > 0);
            if (mb.mb_skip_flag)
                mb.ei_flag = 0;
        }

        mb.mb_type                 = !mb.mb_skip_flag;
        mb.CodedBlockPatternLuma   = !mb.mb_skip_flag;
        mb.CodedBlockPatternChroma = 0;

        if (shr.MbaffFrameFlag && CurrMbAddr % 2 == 0) {
            if (pps.entropy_coding_mode_flag) {
                if (mb.mb_skip_flag) {
                    //get next MB
                    mb_t& nextMb = slice.neighbour.mb_data[mb.mbAddrX + 1];
                    nextMb.p_Slice  = &slice;
                    nextMb.slice_nr = mb.slice_nr;
                    nextMb.mbAddrX  = mb.mbAddrX + 1;
                    nextMb.mb_field_decoding_flag = mb.mb_field_decoding_flag;
                    SyntaxElement se { nextMb };

                    //check_next_mb
                    slice.parser.last_dquant = 0;
                    slice.parser.prescan_skip_read = true;
                    slice.parser.prescan_skip_flag = se.mb_skip_flag();
                    if (slice.parser.prescan_skip_flag)
                        slice.parser.last_dquant = 0;
                    mb.ei_flag = 0;
                    if (!slice.parser.prescan_skip_flag) {
                        slice.parser.prescan_mb_field_decoding_read = true;
                        slice.parser.prescan_mb_field_decoding_flag = se.mb_field_decoding_flag();
                        mb.mb_field_decoding_flag = slice.parser.prescan_mb_field_decoding_flag;
                    }
                }
            } else {
                if (slice.parser.mb_skip_run == 1)
                    mb.mb_field_decoding_flag = slice.parser.partArr[0].next_bits(1);
            }
        }

        if (pps.entropy_coding_mode_flag) {
            if (mb.mb_skip_flag)
                slice.parser.mb_skip_run = 0;
        } else
            --slice.parser.mb_skip_run;

        moreDataFlag = !mb.mb_skip_flag;
    }

    if (moreDataFlag) {
        bool prevMbSkipped = (CurrMbAddr % 2 == 1) ?
            slice.neighbour.mb_data[CurrMbAddr - 1].mb_skip_flag : 0;
        if (shr.MbaffFrameFlag &&
            (CurrMbAddr % 2 == 0 || (CurrMbAddr % 2 == 1 && prevMbSkipped))) {
            if (slice.parser.prescan_mb_field_decoding_read) {
                slice.parser.prescan_mb_field_decoding_read = false;
                mb.mb_field_decoding_flag = slice.parser.prescan_mb_field_decoding_flag;
            } else
                mb.mb_field_decoding_flag = se.mb_field_decoding_flag();
        }
    }

    if (moreDataFlag) {
        mb.mb_type = se.mb_type();
        mb.mb_type += (shr.slice_type == P_slice || shr.slice_type == SP_slice) ? 1 : 0;
        mb.ei_flag = 0;
    }

    this->mb_type(mb.mb_type);

    slice.dec_picture->motion.mb_field_decoding_flag[mb.mbAddrX] = mb.mb_field_decoding_flag;
    if (shr.slice_type != I_slice && shr.slice_type != SI_slice) {
        if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag) {
            shr.num_ref_idx_l0_active_minus1 = ((shr.num_ref_idx_l0_active_minus1 + 1) << 1) - 1;
            shr.num_ref_idx_l1_active_minus1 = ((shr.num_ref_idx_l1_active_minus1 + 1) << 1) - 1;
        }
    }

    if (mb.mb_type == I_PCM)
        return this->parse_i_pcm();

    if (!mb.is_intra_block)
        this->sub_mb_type();

    if (mb.is_intra_block)
        this->mb_pred_intra();
    else
        this->mb_pred_inter();

    this->update_qp(slice.parser.QpY);

    if (mb.mb_type == 0) {
        mb.transform_size_8x8_flag = 0;

        if (shr.slice_type == P_slice)
            return;
        if (slice.parser.mb_skip_run >= 0) {
            if (pps.entropy_coding_mode_flag) {
                slice.parser.is_reset_coeff = true;
                slice.parser.mb_skip_run = -1;
            } else
                memset(mb.nz_coeff, 0, 3 * 16 * sizeof(uint8_t));
            return;
        }
    }

    if (mb.mb_type != I_16x16)
        this->coded_block_pattern();
    if (mb.CodedBlockPatternLuma > 0 || mb.CodedBlockPatternChroma > 0 || mb.mb_type == I_16x16) {
        this->mb_qp_delta();
        this->erc_dpl();
    }
    this->update_qp(slice.parser.QpY);

    //if (mb.CodedBlockPatternLuma > 0 || mb.CodedBlockPatternChroma > 0 || mb.mb_type == I_16x16)
        re.residual();
}

void Parser::Macroblock::mb_type(uint8_t mb_type)
{
    shr_t& shr = slice.header;

    switch (shr.slice_type) {
    case I_slice:
        break;
    case SI_slice:
        if (mb_type < 1)
            return this->mb_type_si_slice(mb_type);
        mb_type -= 1;
        break;
    case P_slice:
    case SP_slice:
        if (mb_type < 6)
            return this->mb_type_p_slice(mb_type);
        mb_type -= 6;
        break;
    case B_slice:
        if (mb_type < 23)
            return this->mb_type_b_slice(mb_type);
        mb_type -= 23;
        break;
    default:
        return;
    }

    return this->mb_type_i_slice(mb_type);
}

void Parser::Macroblock::mb_type_i_slice(uint8_t mb_type)
{
    assert(mb_type <= 25);

    const uint8_t* mb_types = mb_types_i_slice[mb_type];

    mb.is_intra_block          = true;
    mb.mb_type                 = mb_types[0];
    mb.Intra16x16PredMode      = mb_types[2];
    mb.CodedBlockPatternLuma   = mb_types[4];
    mb.CodedBlockPatternChroma = mb_types[3];
}

void Parser::Macroblock::mb_type_si_slice(uint8_t mb_type)
{
    assert(mb_type == 0);

    const uint8_t* mb_types = mb_types_si_slice[mb_type];

    mb.is_intra_block = true;
    mb.mb_type        = mb_types[0];
}

void Parser::Macroblock::mb_type_p_slice(uint8_t mb_type)
{
    assert(mb_type <= 6);

    if (mb_type == 0) {
        mb.is_intra_block = false;
        mb.mb_type        = 0;
        memset(mb.SubMbType,     0, 4 * sizeof(char));
        memset(mb.SubMbPredMode, 0, 4 * sizeof(char));
        return;
    }

    const uint8_t* mb_types = mb_types_p_slice[mb_type - 1];

    mb.is_intra_block = false;
    mb.mb_type        = mb_types[0] == P_8x8ref0 ? P_8x8 : mb_types[0];
    mb.allrefzero     = mb_types[0] == P_8x8ref0;
    memset(mb.SubMbType,     mb_types[0], 4 * sizeof(char));
    memset(mb.SubMbPredMode, mb_types[2], 4 * sizeof(char));
}

void Parser::Macroblock::mb_type_b_slice(uint8_t mb_type)
{
    assert(mb_type <= 23);

    const uint8_t* mb_types = mb_types_b_slice[mb_type];

    mb.is_intra_block = false;
    mb.mb_type        = mb_types[0];
    memset(mb.SubMbType, mb_types[0], 4 * sizeof(char));
    for (int i = 0; i < 4; ++i) {
        mb.SubMbPredMode[i] = mb_types[0] == B_16x16 ? mb_types[2] :
                              mb_types[0] == B_16x8  ? mb_types[2 + (i / 2)] :
                                                       mb_types[2 + (i % 2)];
    }
}

void Parser::Macroblock::parse_i_pcm()
{
    mb.transform_size_8x8_flag = 0;

    // for deblocking filter
    this->update_qp(0);

    memset(mb.nz_coeff, 16, 3 * 16 * sizeof(byte));

    // for CABAC decoding of MB skip flag
    mb.mb_skip_flag = 0;

    //for deblocking filter CABAC
    mb.cbp_blks[0] = 0xFFFF;

    //For CABAC decoding of Dquant
    slice.parser.last_dquant = 0;
    slice.parser.is_reset_coeff = false;
    slice.parser.is_reset_coeff_cr = false;

    auto mv_info = slice.dec_picture->mv_info; 
    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            auto& mv = mv_info[mb.mb.y * 4 + y][mb.mb.x * 4 + x];
            mv.slice_no = mb.slice_nr;
            for (int list = 0; list < 2; list++) {
                mv.ref_pic[list] = NULL;
                mv.ref_idx[list] = -1;
                mv.mv     [list] = {0, 0};
            }
        }
    }

    if (slice.parser.dp_mode && slice.dpB_NotPresent) {
        for (int y = 0; y < 16; y++) {
            for (int x = 0; x < 16; x++)
                slice.decoder.transform->cof[0][y][x] = 1 << (sps.BitDepthY - 1);
        }

        if (sps.chroma_format_idc != YUV400 && !sps.separate_colour_plane_flag) {
            for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
                for (int y = 0; y < sps.MbHeightC; y++) {
                    for (int x = 0; x < sps.MbWidthC; x++)
                        slice.decoder.transform->cof[iCbCr + 1][y][x] = 1 << (sps.BitDepthC - 1);
                }
            }
        }
    } else {
        data_partition_t* dp = &slice.parser.partArr[slice.parser.dp_mode ? 1 : 0];

        if (dp->frame_bitoffset & 7)
            dp->f(8 - (dp->frame_bitoffset & 7));

        for (int y = 0; y < 16; y++) {
            for (int x = 0; x < 16; x++)
                slice.decoder.transform->cof[0][y][x] = dp->f(sps.BitDepthY);
        }

        if (sps.chroma_format_idc != YUV400 && !sps.separate_colour_plane_flag) {
            for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
                for (int y = 0; y < sps.MbHeightC; y++) {
                    for (int x = 0; x < sps.MbWidthC; x++)
                        slice.decoder.transform->cof[iCbCr + 1][y][x] = dp->f(sps.BitDepthC);
                }
            }
        }

        if (pps.entropy_coding_mode_flag)
            dp->de_cabac.init(dp);
    }
}

void Parser::Macroblock::sub_mb_type()
{
    shr_t& shr = slice.header;

    mb.noSubMbPartSizeLessThan8x8Flag = 1;
    mb.transform_size_8x8_flag = 0;

    if (mb.mb_type == P_8x8 || mb.mb_type == B_8x8) {
        for (int mbPartIdx = 0; mbPartIdx < 4; mbPartIdx++) {
            uint8_t sub_mb_type = se.sub_mb_type();

            assert(shr.slice_type != B_slice ? sub_mb_type < 4 : sub_mb_type < 13);

            const uint8_t (*sub_mb_types)[5] = (shr.slice_type != B_slice) ?
                sub_mb_types_p_slice : sub_mb_types_b_slice;
            mb.SubMbType    [mbPartIdx] = sub_mb_types[sub_mb_type][0];
            mb.SubMbPredMode[mbPartIdx] = sub_mb_types[sub_mb_type][2];

            mb.noSubMbPartSizeLessThan8x8Flag &= 
                (mb.SubMbType[mbPartIdx] == P_8x8 || mb.SubMbType[mbPartIdx] == B_8x8) ||
                (mb.SubMbType[mbPartIdx] == B_Direct_8x8 && sps.direct_8x8_inference_flag);
        }
    }
}


void Parser::Macroblock::mb_pred_intra()
{
    auto mv_info = slice.dec_picture->mv_info; 
    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            auto& mv = mv_info[mb.mb.y * 4 + y][mb.mb.x * 4 + x];
            mv.slice_no = mb.slice_nr;
            for (int list = 0; list < 2; list++) {
                mv.ref_pic[list] = NULL;
                mv.ref_idx[list] = -1;
                mv.mv     [list] = {0, 0};
            }
        }
    }

    mb.transform_size_8x8_flag = 0;
    if (pps.transform_8x8_mode_flag && mb.mb_type == I_NxN) {
        mb.transform_size_8x8_flag = se.transform_size_8x8_flag();
        mb.mb_type = mb.transform_size_8x8_flag ? I_8x8 : I_4x4;
    }

    if (mb.mb_type == I_4x4) {
        for (int luma4x4BlkIdx = 0; luma4x4BlkIdx < 16; ++luma4x4BlkIdx) {
            int bx = ((luma4x4BlkIdx / 4) % 2) * 8 + ((luma4x4BlkIdx % 4) % 2) * 4;
            int by = ((luma4x4BlkIdx / 4) / 2) * 8 + ((luma4x4BlkIdx % 4) / 2) * 4;
            int val = se.intra_pred_mode();

            bool    prev_intra4x4_pred_mode_flag = val == -1;
            uint8_t rem_intra4x4_pred_mode       = val;

            uint8_t predIntra4x4PredMode = intra_4x4_pred_mode(mb, bx, by);
            if (prev_intra4x4_pred_mode_flag)
                mb.Intra4x4PredMode[luma4x4BlkIdx] = predIntra4x4PredMode;
            else if (rem_intra4x4_pred_mode < predIntra4x4PredMode)
                mb.Intra4x4PredMode[luma4x4BlkIdx] = rem_intra4x4_pred_mode;
            else
                mb.Intra4x4PredMode[luma4x4BlkIdx] = rem_intra4x4_pred_mode + 1;
        }
    } else if (mb.mb_type == I_8x8) {
        for (int luma8x8BlkIdx = 0; luma8x8BlkIdx < 4; ++luma8x8BlkIdx) {
            int bx = (luma8x8BlkIdx % 2) * 8;
            int by = (luma8x8BlkIdx / 2) * 8;
            int val = se.intra_pred_mode();

            bool    prev_intra8x8_pred_mode_flag = val == -1;
            uint8_t rem_intra8x8_pred_mode       = val;

            uint8_t predIntra8x8PredMode = intra_8x8_pred_mode(mb, bx, by);
            if (prev_intra8x8_pred_mode_flag)
                mb.Intra8x8PredMode[luma8x8BlkIdx] = predIntra8x8PredMode;
            else if (rem_intra8x8_pred_mode < predIntra8x8PredMode)
                mb.Intra8x8PredMode[luma8x8BlkIdx] = rem_intra8x8_pred_mode;
            else
                mb.Intra8x8PredMode[luma8x8BlkIdx] = rem_intra8x8_pred_mode + 1;
        }
    }

    if (sps.ChromaArrayType == 1 || sps.ChromaArrayType == 2) {
        mb.intra_chroma_pred_mode = se.intra_chroma_pred_mode();

        assert(mb.intra_chroma_pred_mode >= IntraPrediction::Intra_Chroma_DC ||
               mb.intra_chroma_pred_mode <= IntraPrediction::Intra_Chroma_Plane);
    }
}

void Parser::Macroblock::mb_pred_inter()
{
    shr_t& shr = slice.header;

    if (shr.slice_type != B_slice && mb.mb_type == P_Skip)
        return this->skip_macroblock();

    if (mb.mb_type != 0) {
        auto mv_info = slice.dec_picture->mv_info; 
        for (int y = 0; y < 4; ++y) {
            for (int x = 0; x < 4; ++x) {
                auto& mv = mv_info[mb.mb.y * 4 + y][mb.mb.x * 4 + x];
                mv.slice_no = mb.slice_nr;
                for (int list = 0; list < 2; list++) {
                    mv.ref_pic[list] = NULL;
                    mv.ref_idx[list] = -1;
                    mv.mv     [list] = {0, 0};
                }
            }
        }
    }

    if (shr.slice_type == B_slice && (mb.mb_type == B_Skip || mb.mb_type == B_8x8)) {
        if (shr.direct_spatial_mv_pred_flag)
            this->get_direct_spatial();
        else
            this->get_direct_temporal();
    }

    if (shr.slice_type == B_slice && mb.mb_type == B_Skip)
        return;

    this->ref_idx_l(LIST_0);
    if (shr.slice_type == B_slice)
        this->ref_idx_l(LIST_1);

    this->mvd_l(LIST_0);
    if (shr.slice_type == B_slice)
        this->mvd_l(LIST_1);

    pic_motion_params** p_mv_info = slice.dec_picture->mv_info;
    // record reference picture Ids for deblocking decisions
    for (int j4 = 0; j4 < 4; j4++) {
        for (int i4 = 0; i4 < 4; ++i4) {
            auto mv_info = &p_mv_info[mb.mb.y * 4 + j4][mb.mb.x * 4 + i4];
            int  ref_idx = mv_info->ref_idx[LIST_0];
            mv_info->ref_pic[LIST_0] = get_ref_pic(mb, slice.RefPicList[LIST_0], ref_idx);
            if (shr.slice_type == B_slice) {
                ref_idx = mv_info->ref_idx[LIST_1];
                mv_info->ref_pic[LIST_1] = get_ref_pic(mb, slice.RefPicList[LIST_1], ref_idx);
            }
        }
    }
}

//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

void Parser::Macroblock::ref_idx_l(int list)
{
    auto mv_info = slice.dec_picture->mv_info;

    int step_h0 = BLOCK_STEP[mb.mb_type][0];
    int step_v0 = BLOCK_STEP[mb.mb_type][1];

    auto ref_idx_l = (list == 0) ? mb.ref_idx_l0 : mb.ref_idx_l1;

    for (int y8 = 0; y8 < 4; y8 += step_v0) {
        for (int x8 = 0; x8 < 4; x8 += step_h0) {      
            int mbPartIdx = 2 * (y8 >> 1) + (x8 >> 1);

            if ((mb.SubMbPredMode[mbPartIdx] == list || mb.SubMbPredMode[mbPartIdx] == BiPred) &&
                (mb.SubMbType[mbPartIdx] != 0)) {
                ref_idx_l[mbPartIdx] = se.ref_idx_l(list, x8, y8);

                for (int y4 = 0; y4 < step_v0; ++y4) {
                    for (int x4 = 0; x4 < step_h0; ++x4) {
                        auto& mv = mv_info[mb.mb.y * 4 + y8 + y4][mb.mb.x * 4 + x8 + x4];
                        mv.ref_idx[list] = ref_idx_l[mbPartIdx];
                    }
                }
            }
        }
    }
}

void Parser::Macroblock::mvd_l(int list)
{
    auto mv_info = slice.dec_picture->mv_info;

    int step_h0 = BLOCK_STEP[mb.mb_type][0];
    int step_v0 = BLOCK_STEP[mb.mb_type][1];

    auto mvd_l = (list == 0) ? mb.mvd_l0 : mb.mvd_l1;

    for (int y8 = 0; y8 < 4; y8 += step_v0) {
        for (int x8 = 0; x8 < 4; x8 += step_h0) {
            int mbPartIdx = 2 * (y8 >> 1) + (x8 >> 1);

            if ((mb.SubMbPredMode[mbPartIdx] == list || mb.SubMbPredMode[mbPartIdx] == BiPred) &&
                (mb.SubMbType[mbPartIdx] != 0)) {
                int step_h4 = BLOCK_STEP[mb.SubMbType[mbPartIdx]][0];
                int step_v4 = BLOCK_STEP[mb.SubMbType[mbPartIdx]][1];
                char cur_ref_idx = mv_info[mb.mb.y * 4 + y8][mb.mb.x * 4 + x8].ref_idx[list];

                for (int y4 = 0; y4 < step_v0; y4 += step_v4) {
                    for (int x4 = 0; x4 < step_h0; x4 += step_h4) {
                        //int subMbPartIdx = 2 * y4 + x4;

                        int16_t curr_mvd[2];
                        for (int compIdx = 0; compIdx < 2; ++compIdx)
                            curr_mvd[compIdx] = se.mvd_l(list, x8 + x4, y8 + y4, compIdx);

                        mv_t pred_mv = this->GetMVPredictor(cur_ref_idx, list, x8 + x4, y8 + y4, step_h4 * 4, step_v4 * 4);
                        mv_t curr_mv;
                        curr_mv.mv_x = curr_mvd[0] + pred_mv.mv_x;
                        curr_mv.mv_y = curr_mvd[1] + pred_mv.mv_y;

                        for (int y2 = 0; y2 < step_v4; ++y2) {
                            for (int x2 = 0; x2 < step_h4; ++x2) {
                                auto& mv = mv_info[mb.mb.y * 4 + y8 + y4 + y2][mb.mb.x * 4 + x8 + x4 + x2];
                                mv.mv[list] = curr_mv;
                                for (int compIdx = 0; compIdx < 2; ++compIdx)
                                    mvd_l[y8 + y4 + y2][x8 + x4 + x2][compIdx] = curr_mvd[compIdx];
                            }
                        }
                    }
                }
            }
        }
    }
}

void Parser::Macroblock::coded_block_pattern()
{
    shr_t& shr = slice.header;

    uint8_t coded_block_pattern = se.coded_block_pattern();
    mb.CodedBlockPatternLuma   = coded_block_pattern % 16;
    mb.CodedBlockPatternChroma = coded_block_pattern / 16;
    if (pps.entropy_coding_mode_flag) {
        if (!mb.CodedBlockPatternLuma && !mb.CodedBlockPatternChroma)
            slice.parser.last_dquant = 0;
    }

#define IS_DIRECT(MB) ((MB)->mb_type == 0 && shr.slice_type == B_slice)
    int need_transform_size_flag =
        mb.CodedBlockPatternLuma > 0 &&
        pps.transform_8x8_mode_flag && !mb.is_intra_block &&
        mb.noSubMbPartSizeLessThan8x8Flag &&
        (!IS_DIRECT(&mb) || sps.direct_8x8_inference_flag);
#undef IS_DIRECT

    if (need_transform_size_flag)
        mb.transform_size_8x8_flag = se.transform_size_8x8_flag();
}

void Parser::Macroblock::mb_qp_delta()
{
    mb.mb_qp_delta = se.mb_qp_delta();

    if (pps.entropy_coding_mode_flag)        
        slice.parser.last_dquant = mb.mb_qp_delta;

    assert(mb.mb_qp_delta >= -(26 + sps.QpBdOffsetY / 2) &&
           mb.mb_qp_delta <=  (25 + sps.QpBdOffsetY / 2));

    slice.parser.QpY =
        ((slice.parser.QpY + mb.mb_qp_delta + 52 + 2 * sps.QpBdOffsetY)
            % (52 + sps.QpBdOffsetY)) - sps.QpBdOffsetY;
}

void Parser::Macroblock::erc_dpl()
{
    if (slice.parser.dp_mode) {
        if (!mb.is_intra_block && slice.dpC_NotPresent)
            mb.dpl_flag = 1;
        if (mb.is_intra_block && slice.dpB_NotPresent) {
            mb.ei_flag  = 1;
            mb.dpl_flag = 1;
        }

        mb_t* mbA = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {-1, 0});
        mb_t* mbB = slice.neighbour.get_mb(&slice, false, mb.mbAddrX, {0, -1});
        mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
        mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

        if (!(pps.constrained_intra_pred_flag && mb.is_intra_block)) {
            if (mbA)
                mb.dpl_flag |= mbA->dpl_flag;
            if (mbB)
                mb.dpl_flag |= mbB->dpl_flag;
        }

        if (mb.dpl_flag) {
            mb.CodedBlockPatternLuma   = 0;
            mb.CodedBlockPatternChroma = 0;
        }
    }
}

// Table 8-15 Specification of QPc as a function of qPi

static const uint8_t QP_SCALE_CR[52] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35,
    35, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 39
};

void Parser::Macroblock::update_qp(int qp)
{
    shr_t& shr = slice.header;

    int QpOffset[2] = { pps.chroma_qp_index_offset, pps.second_chroma_qp_index_offset };

    mb.QpY          = qp;
    mb.qp_scaled[0] = qp + sps.QpBdOffsetY;

    for (int i = 0; i < 2; i++) {
        int8_t QpI = clip3(-(sps.QpBdOffsetC), 51, mb.QpY + QpOffset[i]);
        int8_t QpC = QpI < 30 ? QpI : QP_SCALE_CR[QpI];
        mb.QpC[i]           = QpC;
        mb.qp_scaled[i + 1] = QpC + sps.QpBdOffsetC;

        int8_t QsI = clip3(-(sps.QpBdOffsetC), 51, shr.QsY + QpOffset[i]);
        int8_t QsC = QsI < 30 ? QsI : QP_SCALE_CR[QsI];
        mb.QsC[i]           = QsC;
    }

    mb.TransformBypassModeFlag = (sps.qpprime_y_zero_transform_bypass_flag && mb.qp_scaled[0] == 0);
}


}
}
