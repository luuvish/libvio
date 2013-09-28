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
    { I_PCM  ,          NA, NA, NA, NA }
};

// Table 7-12 Macroblock type with value 0 for SI slices
static const uint8_t mb_types_si_slice[1][5] = {
    { SI     , Intra_4x4  , NA, NA, NA }
    // [1..26] are mb_types_i_slice[0..25]
};

// Table 7-13 mb_t type value 0 to 4 for P and SP slices
static const uint8_t mb_types_p_slice[5][6] = {
//  { P_Skip   , 1, Pred_L0,      NA, 16, 16 },
    { P_16x16  , 1, Pred_L0,      NA, 16, 16 },
    { P_16x8   , 2, Pred_L0, Pred_L0, 16,  8 },
    { P_8x16   , 2, Pred_L0, Pred_L0,  8, 16 },
    { P_8x8    , 4,      NA,      NA,  8,  8 },
    { P_8x8ref0, 4,      NA,      NA,  8,  8 }
    // [5..30] are mb_types_i_slice[0..25]
};

// Table 7-14 mb_t type value 0 to 22 for B slices
static const uint8_t mb_types_b_slice[23][6] = {
//  { B_Skip        , NA, Direct ,      NA,  8,  8 },
    { B_Direct_16x16, NA, Direct ,      NA,  8,  8 },
    { B_16x16       ,  1, Pred_L0,      NA, 16, 16 },
    { B_16x16       ,  1, Pred_L1,      NA, 16, 16 },
    { B_16x16       ,  1, BiPred ,      NA, 16, 16 },
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
    { B_8x8         ,  4,      NA,      NA,  8,  8 }
    // [23..48] are mb_types_i_slice[0..25]
};

// Table 7-17 Sub-macroblock types in P macroblocks
static const uint8_t sub_mb_types_p_slice[4][5] = {
    { P_8x8, 1, Pred_L0, 8, 8 },
    { P_8x4, 2, Pred_L0, 8, 4 },
    { P_4x8, 2, Pred_L0, 4, 8 },
    { P_4x4, 4, Pred_L0, 4, 4 }
};

// Table 7-18 Sub-macroblock types in B macroblocks
static const uint8_t sub_mb_types_b_slice[14][5] = {
//  {           NA, 4, Direct , 4, 4 },
    { B_Direct_8x8, 4, Direct , 4, 4 },
    { B_8x8       , 1, Pred_L0, 8, 8 },
    { B_8x8       , 1, Pred_L1, 8, 8 },
    { B_8x8       , 1, BiPred , 8, 8 },
    { B_8x4       , 2, Pred_L0, 8, 4 },
    { B_4x8       , 2, Pred_L0, 4, 8 },
    { B_8x4       , 2, Pred_L1, 8, 4 },
    { B_4x8       , 2, Pred_L1, 4, 8 },
    { B_8x4       , 2, BiPred , 8, 4 },
    { B_4x8       , 2, BiPred , 4, 8 },
    { B_4x4       , 4, Pred_L0, 4, 4 },
    { B_4x4       , 4, Pred_L1, 4, 4 },
    { B_4x4       , 4, BiPred , 4, 4 }
};

static void interpret_mb_mode_P(mb_t* mb)
{
    static const short ICBPTAB[6] = {0,16,32,15,31,47};
    short mbmode = mb->mb_type;

    if (mbmode < 4) {
        mb->mb_type = mbmode;
        memset(mb->SubMbType, mbmode, 4 * sizeof(char));
        memset(mb->SubMbPredMode,  0, 4 * sizeof(char));
    } else if (mbmode == 4 || mbmode == 5) {
        mb->mb_type = P8x8;
        mb->p_Slice->allrefzero = (mbmode == 5);
    } else if (mbmode == 6) {
        mb->is_intra_block = true;
        mb->mb_type = I4MB;
        memset(mb->SubMbType,   I4MB, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else if (mbmode == 31) {
        mb->is_intra_block = true;
        mb->mb_type = IPCM;
        mb->CodedBlockPatternLuma   = 15;
        mb->CodedBlockPatternChroma =  3;
        mb->Intra16x16PredMode = 0;

        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else {
        mb->is_intra_block = true;
        mb->mb_type = I16MB;
        uint8_t coded_block_pattern = ICBPTAB[((mbmode-7))>>2];
        mb->CodedBlockPatternLuma   = coded_block_pattern % 16;
        mb->CodedBlockPatternChroma = coded_block_pattern / 16;
        mb->Intra16x16PredMode = ((mbmode-7)) & 0x03;
        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    }
}

static void interpret_mb_mode_I(mb_t* mb)
{
    static const short ICBPTAB[6] = {0,16,32,15,31,47};
    short mbmode = mb->mb_type;

    if (mbmode == 0) {
        mb->is_intra_block = true;
        mb->mb_type = I4MB;
        memset(mb->SubMbType,   I4MB, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else if (mbmode == 25) {
        mb->is_intra_block = true;
        mb->mb_type = IPCM;
        mb->CodedBlockPatternLuma   = 15;
        mb->CodedBlockPatternChroma =  3;
        mb->Intra16x16PredMode = 0;

        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else {
        mb->is_intra_block = true;
        mb->mb_type = I16MB;
        uint8_t coded_block_pattern = ICBPTAB[(mbmode-1)>>2];
        mb->CodedBlockPatternLuma   = coded_block_pattern % 16;
        mb->CodedBlockPatternChroma = coded_block_pattern / 16;
        mb->Intra16x16PredMode = (mbmode-1) & 0x03;
        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    }
}

static void interpret_mb_mode_B(mb_t* mb)
{
    static const char offset2pdir16x16[3] = { 0, 1, 2 };
    static const char offset2pdir16x8[9][2] = {
        { 0, 0 }, { 1, 1 }, { 0, 1 },
        { 1, 0 }, { 0, 2 }, { 1, 2 },
        { 2, 0 }, { 2, 1 }, { 2, 2 }
    };
    static const char ICBPTAB[6] = {0,16,32,15,31,47};
    int mbtype = mb->mb_type;

    //--- set mbtype, b8type, and SubMbPredMode ---
    if (mbtype == 0) { // direct
        mb->mb_type = BSKIP_DIRECT;
        memset(mb->SubMbType,     0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, 2, 4 * sizeof(char));
    } else if (mbtype < 4) { // 16x16
        mb->mb_type = P16x16;
        memset(mb->SubMbType,                      mb->mb_type, 4 * sizeof(char));
        memset(mb->SubMbPredMode, offset2pdir16x16[mbtype - 1], 4 * sizeof(char));
    } else if (mbtype < 22) {
        mb->mb_type = (mbtype % 2 == 0 ? P16x8 : P8x16);
        memset(mb->SubMbType, mb->mb_type, 4 * sizeof(char));
        for (int i = 0; i < 4; i++)
            mb->SubMbPredMode[i] = offset2pdir16x8[(mbtype - 4) / 2][mbtype % 2 == 0 ? i / 2 : i % 2];
    } else if (mbtype == 22) { // 8x8(+split)
        mb->mb_type = P8x8;       // SubMbType and pdir is transmitted in additional codewords
    } else if (mbtype == 23) { // intra4x4
        mb->is_intra_block = true;
        mb->mb_type = I4MB;
        memset(mb->SubMbType,   I4MB, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else if (mbtype > 23 && mbtype < 48) { // intra16x16
        mb->is_intra_block = true;
        mb->mb_type = I16MB;
        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));

        uint8_t coded_block_pattern = ICBPTAB[(mbtype-24)>>2];
        mb->CodedBlockPatternLuma   = coded_block_pattern % 16;
        mb->CodedBlockPatternChroma = coded_block_pattern / 16;
        mb->Intra16x16PredMode = (mbtype-24) & 0x03;
    } else if (mbtype == 48) {
        mb->is_intra_block = true;
        mb->mb_type = IPCM;
        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));

        mb->CodedBlockPatternLuma   = 15;
        mb->CodedBlockPatternChroma =  3;
        mb->Intra16x16PredMode = 0;
    }
}

static void interpret_mb_mode_SI(mb_t* mb)
{
    const int ICBPTAB[6] = {0,16,32,15,31,47};
    short mbmode = mb->mb_type;

    if (mbmode == 0) {
        mb->is_intra_block = true;
        mb->mb_type = SI4MB;
        memset(mb->SubMbType,   I4MB, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else if (mbmode == 1) {
        mb->is_intra_block = true;
        mb->mb_type = I4MB;
        memset(mb->SubMbType,   I4MB, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else if (mbmode == 26) {
        mb->is_intra_block = true;
        mb->mb_type = IPCM;
        mb->CodedBlockPatternLuma   = 15;
        mb->CodedBlockPatternChroma =  3;
        mb->Intra16x16PredMode = 0;
        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    } else {
        mb->is_intra_block = true;
        mb->mb_type = I16MB;
        uint8_t coded_block_pattern = ICBPTAB[(mbmode - 2) >> 2];
        mb->CodedBlockPatternLuma   = coded_block_pattern % 16;
        mb->CodedBlockPatternChroma = coded_block_pattern / 16;
        mb->Intra16x16PredMode = (mbmode - 2) & 0x03;
        memset(mb->SubMbType,      0, 4 * sizeof(char));
        memset(mb->SubMbPredMode, -1, 4 * sizeof(char));
    }
}

void interpret_mb_mode(mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    switch (slice->slice_type) {
    case P_SLICE: 
    case SP_SLICE:
        interpret_mb_mode_P(mb);
        break;
    case B_SLICE:
        interpret_mb_mode_B(mb);
        break;
    case I_SLICE: 
        interpret_mb_mode_I(mb);
        break;
    case SI_SLICE: 
        interpret_mb_mode_SI(mb);
        break;
    default:
        printf("Unsupported slice type\n");
        break;
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
    int CurrMbAddr = mb.mbAddrX;
    bool moreDataFlag = 1;

    data_partition_t* dp = &slice.parser.partArr[0];

    if (slice.slice_type != I_slice && slice.slice_type != SI_slice) {
        if (pps.entropy_coding_mode_flag) {
            if (slice.prescan_skip_read) {
                slice.prescan_skip_read = false;
                mb.mb_skip_flag = slice.prescan_skip_flag;
            } else {
                mb.mb_skip_flag = se.mb_skip_flag();
                if (mb.mb_skip_flag)
                    slice.last_dquant = 0;
                mb.ei_flag = 0;
            }
        } else {
            if (slice.mb_skip_run == -1)
                slice.mb_skip_run = se.mb_skip_run();
            mb.mb_skip_flag = (slice.mb_skip_run > 0);
            if (mb.mb_skip_flag)
                mb.ei_flag = 0;
        }

        mb.mb_type                 = !mb.mb_skip_flag;
        mb.CodedBlockPatternLuma   = !mb.mb_skip_flag;
        mb.CodedBlockPatternChroma = 0;

        if (slice.MbaffFrameFlag && CurrMbAddr % 2 == 0) {
            if (pps.entropy_coding_mode_flag) {
                if (mb.mb_skip_flag) {
                    //get next MB
                    ++slice.current_mb_nr;

                    mb_t *currMB = &slice.mb_data[slice.current_mb_nr];
                    currMB->p_Slice  = &slice;
                    currMB->slice_nr = slice.current_slice_nr;
                    currMB->mb_field_decoding_flag = slice.mb_data[slice.current_mb_nr-1].mb_field_decoding_flag;
                    currMB->mbAddrX  = slice.current_mb_nr;

                    SyntaxElement se { *currMB };

                    //check_next_mb
                    slice.last_dquant = 0;
                    slice.prescan_skip_read = true;
                    slice.prescan_skip_flag = se.mb_skip_flag();
                    if (slice.prescan_skip_flag)
                        slice.last_dquant = 0;
                    mb.ei_flag = 0;
                    if (!slice.prescan_skip_flag) {
                        slice.prescan_mb_field_decoding_read = true;
                        slice.prescan_mb_field_decoding_flag = se.mb_field_decoding_flag();
                        slice.mb_data[slice.current_mb_nr-1].mb_field_decoding_flag = slice.prescan_mb_field_decoding_flag;
                    }

                    //reset
                    slice.current_mb_nr--;
                }
            } else {
                if (slice.mb_skip_run == 1)
                    mb.mb_field_decoding_flag = dp->next_bits(1);
            }
        }

        if (pps.entropy_coding_mode_flag) {
            if (mb.mb_skip_flag)
                slice.mb_skip_run = 0;
        } else
            slice.mb_skip_run--;

        moreDataFlag = !mb.mb_skip_flag;
    }

    if (moreDataFlag) {
        bool prevMbSkipped = (CurrMbAddr % 2 == 1) ?
            slice.mb_data[CurrMbAddr - 1].mb_skip_flag : 0;
        if (slice.MbaffFrameFlag &&
            (CurrMbAddr % 2 == 0 || (CurrMbAddr % 2 == 1 && prevMbSkipped))) {
            if (slice.prescan_mb_field_decoding_read) {
                slice.prescan_mb_field_decoding_read = false;
                mb.mb_field_decoding_flag = slice.prescan_mb_field_decoding_flag;
            } else
                mb.mb_field_decoding_flag = se.mb_field_decoding_flag();
        }
    }

    if (moreDataFlag) {
        mb.mb_type = se.mb_type();
        mb.ei_flag = 0;
    }

    if (slice.slice_type != I_slice && slice.slice_type != SI_slice) {
        if (slice.MbaffFrameFlag && mb.mb_field_decoding_flag) {
            slice.num_ref_idx_l0_active_minus1 = ((slice.num_ref_idx_l0_active_minus1 + 1) << 1) - 1;
            slice.num_ref_idx_l1_active_minus1 = ((slice.num_ref_idx_l1_active_minus1 + 1) << 1) - 1;
        }
    }

    interpret_mb_mode(&mb);

    if (slice.slice_type != B_slice)
        mb.noSubMbPartSizeLessThan8x8Flag = 1;
    slice.dec_picture->motion.mb_field_decoding_flag[mb.mbAddrX] = mb.mb_field_decoding_flag;

    if (mb.mb_type == IPCM)
        this->parse_i_pcm();
    else if (mb.mb_type == PSKIP || mb.mb_type == BSKIP_DIRECT)
        this->parse_skip();
    else if (mb.is_intra_block)
        this->parse_intra();
    else
        this->parse_inter();
}


void Parser::Macroblock::parse_i_pcm()
{
    mb.noSubMbPartSizeLessThan8x8Flag = 1;
    mb.transform_size_8x8_flag = 0;

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

    if (slice.dp_mode && slice.dpB_NotPresent) {
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
        data_partition_t* dp = &slice.parser.partArr[slice.dp_mode ? 1 : 0];

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

void Parser::Macroblock::parse_skip()
{
    mb.transform_size_8x8_flag = 0;
    if (pps.constrained_intra_pred_flag)
        mb.is_intra_block = 0;

    auto mv_info = slice.dec_picture->mv_info;
    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            auto& mv = mv_info[mb.mb.y * 4 + y][mb.mb.x * 4 + x];
            mv.slice_no = mb.slice_nr;
            if (slice.slice_type != B_slice) {
                mv.ref_pic[LIST_1] = NULL;
                mv.ref_idx[LIST_1] = -1;
                mv.mv     [LIST_1] = {0, 0};
            }
        }
    }

    if (slice.slice_type == B_slice) {
        if (slice.direct_spatial_mv_pred_flag)
            this->get_direct_spatial();
        else
            this->get_direct_temporal();

        mb.noSubMbPartSizeLessThan8x8Flag = sps.direct_8x8_inference_flag;

        if (slice.mb_skip_run >= 0) {
            if (pps.entropy_coding_mode_flag) {
                slice.is_reset_coeff = true;
                slice.mb_skip_run = -1;
            } else
                memset(mb.nz_coeff, 0, 3 * 16 * sizeof(uint8_t));
        } else {
            this->parse_cbp_qp();
            re.residual();
        }
    } else
        this->skip_macroblock();
}

void Parser::Macroblock::parse_intra()
{
    if (mb.mb_type != I4MB)
        mb.noSubMbPartSizeLessThan8x8Flag = 1;

    mb.transform_size_8x8_flag = 0;
    if (mb.mb_type == I4MB && pps.transform_8x8_mode_flag) {
        mb.transform_size_8x8_flag = se.transform_size_8x8_flag();
        if (mb.transform_size_8x8_flag) {
            mb.mb_type = I8MB;
            for (int i = 0; i < 4; i++) {
                mb.SubMbType    [i] = I8MB;
                mb.SubMbPredMode[i] = -1;
            }
        }
    }

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

    this->parse_ipred_modes();
    this->parse_cbp_qp();

    re.residual();
}

void Parser::Macroblock::parse_inter()
{
    static const char p_v2b8 [ 4] = { 4, 5, 6, 7 };
    static const char p_v2pd [ 4] = { 0, 0, 0, 0 };
    static const char b_v2b8 [13] = { 0, 4, 4, 4, 5, 6, 5, 6, 5, 6, 7, 7, 7 };
    static const char b_v2pd [13] = { 2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1, 2 };

    mb.noSubMbPartSizeLessThan8x8Flag = 1;
    mb.transform_size_8x8_flag = 0;
    if (mb.mb_type == P8x8) {
        for (int mbPartIdx = 0; mbPartIdx < 4; mbPartIdx++) {
            int val = se.sub_mb_type();

            //const uint8_t (*sub_mb_types)[5] = slice->slice_type == B_slice ?
            //    sub_mb_types_b_slice : sub_mb_types_p_slice;

            if (slice.slice_type != B_slice) {
                assert(val < 4);
                mb.SubMbType    [mbPartIdx] = p_v2b8[val];
                mb.SubMbPredMode[mbPartIdx] = p_v2pd[val];
            } else {
                assert(val < 13);
                mb.SubMbType    [mbPartIdx] = b_v2b8[val];
                mb.SubMbPredMode[mbPartIdx] = b_v2pd[val];
            }

            //set noSubMbPartSizeLessThan8x8Flag for P8x8 mode
            mb.noSubMbPartSizeLessThan8x8Flag &= 
                (mb.SubMbType[mbPartIdx] == 0 && sps.direct_8x8_inference_flag) ||
                (mb.SubMbType[mbPartIdx] == P8x8);
        }
    }

    if (pps.constrained_intra_pred_flag)
        mb.is_intra_block = 0;

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

    this->parse_motion_info();
    this->parse_cbp_qp();

    re.residual();
}


void Parser::Macroblock::parse_ipred_modes()
{
    if (mb.mb_type == I8MB)
        this->parse_ipred_8x8_modes();
    else if (mb.mb_type == I4MB)
        this->parse_ipred_4x4_modes();

    if (sps.chroma_format_idc != YUV400 && sps.chroma_format_idc != YUV444) {
        mb.intra_chroma_pred_mode = se.intra_chroma_pred_mode();

        assert(mb.intra_chroma_pred_mode >= IntraPrediction::Intra_Chroma_DC ||
               mb.intra_chroma_pred_mode <= IntraPrediction::Intra_Chroma_Plane);
    }
}

void Parser::Macroblock::parse_ipred_4x4_modes()
{
    for (int luma4x4BlkIdx = 0; luma4x4BlkIdx < 16; luma4x4BlkIdx++) {
        int bx = ((luma4x4BlkIdx / 4) % 2) * 8 + ((luma4x4BlkIdx % 4) % 2) * 4;
        int by = ((luma4x4BlkIdx / 4) / 2) * 8 + ((luma4x4BlkIdx % 4) / 2) * 4;
        int val = se.intra_pred_mode();

        bool    prev_intra4x4_pred_mode_flag = val == -1;
        uint8_t rem_intra4x4_pred_mode       = val;

        nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx - 1, by});
        nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx, by - 1});

        nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
        nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

        //get from array and decode
        if (pps.constrained_intra_pred_flag) {
            nbA.mb = nbA.mb && nbA.mb->is_intra_block ? nbA.mb : nullptr;
            nbB.mb = nbB.mb && nbB.mb->is_intra_block ? nbB.mb : nullptr;
        }
        // !! KS: not sure if the following is still correct...
        if (slice.slice_type == SI_slice) { // need support for MBINTLC1
            nbA.mb = nbA.mb && nbA.mb->mb_type == SI4MB ? nbA.mb : nullptr;
            nbB.mb = nbB.mb && nbB.mb->mb_type == SI4MB ? nbB.mb : nullptr;
        }

        bool dcPredModePredictedFlag = !nbA.mb || !nbB.mb;

        int scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };
        uint8_t intraMxMPredModeA = IntraPrediction::Intra_4x4_DC;
        uint8_t intraMxMPredModeB = IntraPrediction::Intra_4x4_DC;
        if (!dcPredModePredictedFlag) {
            if (nbA.mb->mb_type == I8MB)
                intraMxMPredModeA = nbA.mb->Intra8x8PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4] / 4];
            else if (nbA.mb->mb_type == I4MB)
                intraMxMPredModeA = nbA.mb->Intra4x4PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4]];
            if (nbB.mb->mb_type == I8MB)
                intraMxMPredModeB = nbB.mb->Intra8x8PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4] / 4];
            else if (nbB.mb->mb_type == I4MB)
                intraMxMPredModeB = nbB.mb->Intra4x4PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4]];
        }

        uint8_t predIntra4x4PredMode = min(intraMxMPredModeA, intraMxMPredModeB);
        if (prev_intra4x4_pred_mode_flag)
            mb.Intra4x4PredMode[luma4x4BlkIdx] = predIntra4x4PredMode;
        else if (rem_intra4x4_pred_mode < predIntra4x4PredMode)
            mb.Intra4x4PredMode[luma4x4BlkIdx] = rem_intra4x4_pred_mode;
        else
            mb.Intra4x4PredMode[luma4x4BlkIdx] = rem_intra4x4_pred_mode + 1;

        //if ((luma4x4BlkIdx % 4) == 0)
        //    mb.Intra8x8PredMode[luma4x4BlkIdx / 4] = mb.Intra4x4PredMode[luma4x4BlkIdx];
    }
}

void Parser::Macroblock::parse_ipred_8x8_modes()
{
    for (int luma8x8BlkIdx = 0; luma8x8BlkIdx < 4; luma8x8BlkIdx++) {
        int bx = (luma8x8BlkIdx % 2) * 8;
        int by = (luma8x8BlkIdx / 2) * 8;
        int val = se.intra_pred_mode();

        bool    prev_intra8x8_pred_mode_flag = val == -1;
        uint8_t rem_intra8x8_pred_mode       = val;

        nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx - 1, by});
        nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx, by - 1});
        nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
        nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

        //get from array and decode
        if (pps.constrained_intra_pred_flag) {
            nbA.mb = nbA.mb && nbA.mb->is_intra_block ? nbA.mb : nullptr;
            nbB.mb = nbB.mb && nbB.mb->is_intra_block ? nbB.mb : nullptr;
        }

        bool dcPredModePredictedFlag = !nbA.mb || !nbB.mb;

        int scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };
        uint8_t intraMxMPredModeA = IntraPrediction::Intra_8x8_DC;
        uint8_t intraMxMPredModeB = IntraPrediction::Intra_8x8_DC;
        if (!dcPredModePredictedFlag) {
            if (nbA.mb->mb_type == I8MB)
                intraMxMPredModeA = nbA.mb->Intra8x8PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4] / 4];
            else if (nbA.mb->mb_type == I4MB)
                intraMxMPredModeA = nbA.mb->Intra4x4PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4]];
            if (nbB.mb->mb_type == I8MB)
                intraMxMPredModeB = nbB.mb->Intra8x8PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4] / 4];
            else if (nbB.mb->mb_type == I4MB)
                intraMxMPredModeB = nbB.mb->Intra4x4PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4]];
        }

        uint8_t predIntra8x8PredMode = min(intraMxMPredModeA, intraMxMPredModeB);
        if (prev_intra8x8_pred_mode_flag)
            mb.Intra8x8PredMode[luma8x8BlkIdx] = predIntra8x8PredMode;
        else if (rem_intra8x8_pred_mode < predIntra8x8PredMode)
            mb.Intra8x8PredMode[luma8x8BlkIdx] = rem_intra8x8_pred_mode;
        else
            mb.Intra8x8PredMode[luma8x8BlkIdx] = rem_intra8x8_pred_mode + 1;

        //mb.Intra4x4PredMode[luma8x8BlkIdx * 4    ] = mb.Intra8x8PredMode[luma8x8BlkIdx];
        //mb.Intra4x4PredMode[luma8x8BlkIdx * 4 + 1] = mb.Intra8x8PredMode[luma8x8BlkIdx];
        //mb.Intra4x4PredMode[luma8x8BlkIdx * 4 + 2] = mb.Intra8x8PredMode[luma8x8BlkIdx];
        //mb.Intra4x4PredMode[luma8x8BlkIdx * 4 + 3] = mb.Intra8x8PredMode[luma8x8BlkIdx];
    }
}


void Parser::Macroblock::parse_motion_info()
{
    if (slice.slice_type == B_slice && mb.mb_type == P8x8) {
        if (slice.direct_spatial_mv_pred_flag)
            this->get_direct_spatial(false);
        else
            this->get_direct_temporal(true);
    }

    this->parse_ref_pic_idx(LIST_0);
    if (slice.slice_type == B_slice)
        this->parse_ref_pic_idx(LIST_1);

    this->parse_motion_vectors(LIST_0);
    if (slice.slice_type == B_slice)
        this->parse_motion_vectors(LIST_1);

    int list_offset = slice.MbaffFrameFlag && mb.mb_field_decoding_flag ?
                      mb.mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = slice.listX[LIST_0 + list_offset];
    storable_picture **list1 = slice.listX[LIST_1 + list_offset];
    pic_motion_params **p_mv_info = slice.dec_picture->mv_info;
    // record reference picture Ids for deblocking decisions
    for (int j4 = 0; j4 < 4; j4++) {
        for (int i4 = 0; i4 < 4; ++i4) {
            auto mv_info = &p_mv_info[mb.mb.y * 4 + j4][mb.mb.x * 4 + i4];
            short ref_idx = mv_info->ref_idx[LIST_0];
            mv_info->ref_pic[LIST_0] = (ref_idx >= 0) ? list0[ref_idx] : NULL;
            if (slice.slice_type == B_slice) {
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

void Parser::Macroblock::parse_ref_pic_idx(int list)
{
    pic_motion_params **mv_info = slice.dec_picture->mv_info;

    int step_h0 = BLOCK_STEP[mb.mb_type][0];
    int step_v0 = BLOCK_STEP[mb.mb_type][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {      
            int k = 2 * (j0 >> 1) + (i0 >> 1);

            if ((mb.SubMbPredMode[k] == list || mb.SubMbPredMode[k] == BI_PRED) &&
                mb.SubMbType[k] != 0) {
                uint8_t refframe = se.ref_idx_l(list, i0, j0);
                for (int j = 0; j < step_v0; j++) {
                    for (int i = 0; i < step_h0; i++)
                        mv_info[mb.mb.y * 4 + j0 + j][mb.mb.x * 4 + i0 + i].ref_idx[list] = refframe;
                }
            }
        }
    }
}

void Parser::Macroblock::parse_motion_vectors(int list)
{
    pic_motion_params **mv_info = slice.dec_picture->mv_info;

    int step_h0 = BLOCK_STEP[mb.mb_type][0];
    int step_v0 = BLOCK_STEP[mb.mb_type][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int kk = 2 * (j0 >> 1) + (i0 >> 1);
            int step_h4 = BLOCK_STEP[mb.SubMbType[kk]][0];
            int step_v4 = BLOCK_STEP[mb.SubMbType[kk]][1];
            char cur_ref_idx = mv_info[mb.mb.y * 4 + j0][mb.mb.x * 4 + i0].ref_idx[list];

            if ((mb.SubMbPredMode[kk] == list || mb.SubMbPredMode[kk] == BI_PRED) &&
                (mb.SubMbType[kk] != 0)) { //has forward vector
                for (int j = j0; j < j0 + step_v0; j += step_v4) {
                    for (int i = i0; i < i0 + step_h0; i += step_h4)
                        this->parse_motion_vector(list, step_h4 * 4, step_v4 * 4, i, j, cur_ref_idx);
                }
            }
        }
    }
}

void Parser::Macroblock::parse_motion_vector(int list, int step_h4, int step_v4, int i, int j, char cur_ref_idx)
{
    //const uint8_t (*sub_mb_types)[5] = mb.p_Slice->slice_type == B_slice ?
    //                                   sub_mb_types_b_slice : sub_mb_types_p_slice;
    //int step_h4 = sub_mb_types[mb.sub_mb_type[kk]][3];
    //int step_v4 = sub_mb_types[mb.sub_mb_type[kk]][4];

    int16_t curr_mvd[2];
    for (int k = 0; k < 2; ++k)
        curr_mvd[k] = se.mvd_l(list, i, j, k);

    mv_t pred_mv = this->GetMVPredictor(cur_ref_idx, list, i, j, step_h4, step_v4);
    mv_t curr_mv;
    curr_mv.mv_x = curr_mvd[0] + pred_mv.mv_x;
    curr_mv.mv_y = curr_mvd[1] + pred_mv.mv_y;

    auto mv_info = slice.dec_picture->mv_info;
    int i4 = mb.mb.x * 4 + i;
    int j4 = mb.mb.y * 4 + j;
    for (int jj = 0; jj < step_v4 / 4; ++jj) {
        for (int ii = 0; ii < step_h4 / 4; ++ii) {
            mv_info[jj + j4][ii + i4].mv[list] = curr_mv;
            auto mvd = list == 0 ? mb.mvd_l0 : mb.mvd_l1;
            mvd[jj + j][ii + i][0] = curr_mvd[0];
            mvd[jj + j][ii + i][1] = curr_mvd[1];
        }
    }
}

void Parser::Macroblock::parse_cbp_qp()
{
#define IS_I16MB(MB)  ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)
#define IS_DIRECT(MB) ((MB)->mb_type == 0 && slice.slice_type == B_SLICE)
    // read CBP if not new intra mode
    if (!IS_I16MB(&mb)) {
        uint8_t coded_block_pattern = se.coded_block_pattern();
        mb.CodedBlockPatternLuma   = coded_block_pattern % 16;
        mb.CodedBlockPatternChroma = coded_block_pattern / 16;
        if (pps.entropy_coding_mode_flag) {
            if (!mb.CodedBlockPatternLuma && !mb.CodedBlockPatternChroma)
                slice.last_dquant = 0;
        }

        //============= Transform size flag for INTER MBs =============
        //-------------------------------------------------------------
        int need_transform_size_flag =
            ((mb.mb_type >= 1 && mb.mb_type <= 3) ||
             (IS_DIRECT(&mb) && sps.direct_8x8_inference_flag) ||
             mb.noSubMbPartSizeLessThan8x8Flag) &&
            (mb.mb_type != I8MB) && (mb.mb_type != I4MB) &&
            (mb.CodedBlockPatternLuma) && pps.transform_8x8_mode_flag;

        if (need_transform_size_flag)
            mb.transform_size_8x8_flag = se.transform_size_8x8_flag();
    }

    //=====   DQUANT   =====
    //----------------------
    // Delta quant only if nonzero coeffs
    if (IS_I16MB(&mb) || mb.CodedBlockPatternLuma != 0 || mb.CodedBlockPatternChroma != 0) {
        mb.mb_qp_delta = se.mb_qp_delta();

        if (pps.entropy_coding_mode_flag)        
            slice.last_dquant = mb.mb_qp_delta;

        assert(mb.mb_qp_delta >= -(26 + sps.QpBdOffsetY / 2) &&
               mb.mb_qp_delta <=  (25 + sps.QpBdOffsetY / 2));

        slice.SliceQpY =
            ((slice.SliceQpY + mb.mb_qp_delta + 52 + 2 * sps.QpBdOffsetY)
                % (52 + sps.QpBdOffsetY)) - sps.QpBdOffsetY;

        if (slice.dp_mode) {
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
#undef IS_DIRECT
#undef IS_I16MB

    slice.parser.update_qp(mb, slice.SliceQpY);
}

// Table 8-15 Specification of QPc as a function of qPi

static const uint8_t QP_SCALE_CR[52] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35,
    35, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 39
};

void Parser::update_qp(mb_t& mb, int qp)
{
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;
    pps_t& pps = *slice.active_pps;

    int QpOffset[2] = { pps.chroma_qp_index_offset, pps.second_chroma_qp_index_offset };

    mb.QpY          = qp;
    mb.qp_scaled[0] = qp + sps.QpBdOffsetY;

    for (int i = 0; i < 2; i++) {
        int8_t QpI = clip3(-(sps.QpBdOffsetC), 51, mb.QpY + QpOffset[i]);
        int8_t QpC = QpI < 30 ? QpI : QP_SCALE_CR[QpI];
        mb.QpC[i]           = QpC;
        mb.qp_scaled[i + 1] = QpC + sps.QpBdOffsetC;

        int8_t QsI = clip3(-(sps.QpBdOffsetC), 51, slice.QsY + QpOffset[i]);
        int8_t QsC = QsI < 30 ? QsI : QP_SCALE_CR[QsI];
        mb.QsC[i]           = QsC;
    }

    mb.TransformBypassModeFlag = (sps.qpprime_y_zero_transform_bypass_flag && mb.qp_scaled[0] == 0);
}


}
}
