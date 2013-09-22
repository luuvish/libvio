#include <math.h>

#include "global.h"
#include "slice.h"
#include "dpb.h"
#include "bitstream_cabac.h"
#include "data_partition.h"
#include "macroblock.h"
#include "fmo.h"
#include "image.h"
#include "neighbour.h"
#include "transform.h"
#include "mv_prediction.h"
#include "intra_prediction.h"
#include "inter_prediction.h"

#include "mb_read_syntax.h"


using namespace vio::h264;
using vio::h264::cabac_engine_t;


// Table 7-11 mb_t types for I slices
const uint8_t mb_types_i_slice[26][5] = {
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
const uint8_t mb_types_si_slice[1][5] = {
    { SI     , Intra_4x4  , NA, NA, NA }
    // [1..26] are mb_types_i_slice[0..25]
};

// Table 7-13 mb_t type value 0 to 4 for P and SP slices
const uint8_t mb_types_p_slice[5][6] = {
//  { P_Skip   , 1, Pred_L0,      NA, 16, 16 },
    { P_16x16  , 1, Pred_L0,      NA, 16, 16 },
    { P_16x8   , 2, Pred_L0, Pred_L0, 16,  8 },
    { P_8x16   , 2, Pred_L0, Pred_L0,  8, 16 },
    { P_8x8    , 4,      NA,      NA,  8,  8 },
    { P_8x8ref0, 4,      NA,      NA,  8,  8 }
    // [5..30] are mb_types_i_slice[0..25]
};

// Table 7-14 mb_t type value 0 to 22 for B slices
const uint8_t mb_types_b_slice[23][6] = {
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
const uint8_t sub_mb_types_p_slice[4][5] = {
    { P_8x8, 1, Pred_L0, 8, 8 },
    { P_8x4, 2, Pred_L0, 8, 4 },
    { P_4x8, 2, Pred_L0, 4, 8 },
    { P_4x4, 4, Pred_L0, 4, 4 }
};

// Table 7-18 Sub-macroblock types in B macroblocks
const uint8_t sub_mb_types_b_slice[14][5] = {
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

void macroblock_t::interpret_mb_mode()
{
    slice_t* slice = this->p_Slice;
    switch (slice->slice_type) {
    case P_SLICE: 
    case SP_SLICE:
        interpret_mb_mode_P(this);
        break;
    case B_SLICE:
        interpret_mb_mode_B(this);
        break;
    case I_SLICE: 
        interpret_mb_mode_I(this);
        break;
    case SI_SLICE: 
        interpret_mb_mode_SI(this);
        break;
    default:
        printf("Unsupported slice type\n");
        break;
    }
}


static void init_macroblock(mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    pic_motion_params** mv_info = slice->dec_picture->mv_info; 
    int slice_no = slice->current_slice_nr;

    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            auto mv = &mv_info[mb->mb.y * 4 + y][mb->mb.x * 4 + x];
            mv->slice_no = slice_no;
            for (int list = 0; list < 2; list++) {
                mv->ref_pic[list] = NULL;
                mv->ref_idx[list] = -1;
                mv->mv     [list] = zero_mv;
            }
        }
    }
}


static void skip_macroblock(mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    pps_t* pps = slice->active_pps;
    int list_offset = slice->MbaffFrameFlag && mb->mb_field_decoding_flag ?
                      mb->mbAddrX % 2 ? 4 : 2 : 0;

    int          a_ref_idx = 0;
    int          b_ref_idx = 0;
    MotionVector a_mv;
    MotionVector b_mv;

    PixelPos pos[4];
    get_neighbors(mb, pos, 0, 0, MB_BLOCK_SIZE);
    if (pos[0].available) {
        auto mv_info = &slice->dec_picture->mv_info[pos[0].pos_y][pos[0].pos_x];
        a_ref_idx = mv_info->ref_idx[LIST_0];
        a_mv      = mv_info->mv     [LIST_0];

        if (slice->MbaffFrameFlag) {
            auto mb_a = &mb->p_Vid->mb_data[pos[0].mb_addr];
            if (mb->mb_field_decoding_flag && !mb_a->mb_field_decoding_flag) {
                a_ref_idx *= 2;
                a_mv.mv_y /= 2;
            }
            if (!mb->mb_field_decoding_flag && mb_a->mb_field_decoding_flag) {
                a_ref_idx >>= 1;
                a_mv.mv_y *= 2;
            }
        }
    }
    if (pos[1].available) {
        auto mv_info = &slice->dec_picture->mv_info[pos[1].pos_y][pos[1].pos_x];
        b_ref_idx = mv_info->ref_idx[LIST_0];
        b_mv      = mv_info->mv     [LIST_0];

        if (slice->MbaffFrameFlag) {
            auto mb_b = &mb->p_Vid->mb_data[pos[1].mb_addr];
            if (mb->mb_field_decoding_flag && !mb_b->mb_field_decoding_flag) {
                b_ref_idx *= 2;
                b_mv.mv_y /= 2;
            }
            if (!mb->mb_field_decoding_flag && mb_b->mb_field_decoding_flag) {
                b_ref_idx >>= 1;
                b_mv.mv_y *= 2;
            }
        }
    }

    mb->CodedBlockPatternLuma   = 0;
    mb->CodedBlockPatternChroma = 0;
    if (!pps->entropy_coding_mode_flag)
        memset(mb->nz_coeff, 0, 3 * 16 * sizeof(uint8_t));

    pic_motion_params** mv_info = slice->dec_picture->mv_info;
    storable_picture* cur_pic = slice->listX[list_offset][0];

    bool zeroMotionA = !pos[0].available || (a_ref_idx == 0 && a_mv.mv_x == 0 && a_mv.mv_y == 0);
    bool zeroMotionB = !pos[1].available || (b_ref_idx == 0 && b_mv.mv_x == 0 && b_mv.mv_y == 0);
    MotionVector pred_mv;
    if (zeroMotionA || zeroMotionB)
        pred_mv = zero_mv;
    else
        GetMVPredictor(mb, pos, &pred_mv, 0, mv_info, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            auto mv = &mv_info[mb->mb.y * 4 + y][mb->mb.x * 4 + x];
            mv->ref_pic[LIST_0] = cur_pic;
            mv->ref_idx[LIST_0] = 0;
            mv->mv     [LIST_0] = pred_mv;
        }
    }
}




void macroblock_t::parse()
{
    slice_t *slice = this->p_Slice;
    pps_t *pps = slice->active_pps;
    int CurrMbAddr = this->mbAddrX;
    bool moreDataFlag = 1;

    data_partition_t *dp = &slice->partArr[0];

    if (slice->slice_type != I_slice && slice->slice_type != SI_slice) {
        if (pps->entropy_coding_mode_flag) {
            if (slice->prescan_skip_read) {
                slice->prescan_skip_read = false;
                this->mb_skip_flag = slice->prescan_skip_flag;
            } else {
                this->mb_skip_flag = parse_mb_skip_flag(this);
                if (this->mb_skip_flag)
                    slice->last_dquant = 0;
                this->ei_flag = 0;
            }
        } else {
            if (slice->mb_skip_run == -1)
                slice->mb_skip_run = parse_mb_skip_run(this);
            this->mb_skip_flag = (slice->mb_skip_run > 0);
            if (this->mb_skip_flag)
                this->ei_flag = 0;
        }

        this->mb_type                 = !this->mb_skip_flag;
        this->CodedBlockPatternLuma   = !this->mb_skip_flag;
        this->CodedBlockPatternChroma = 0;

        if (slice->MbaffFrameFlag && CurrMbAddr % 2 == 0) {
            if (pps->entropy_coding_mode_flag) {
                if (this->mb_skip_flag) {
                    //get next MB
                    ++slice->current_mb_nr;

                    mb_t *currMB;
                    currMB = &slice->mb_data[slice->current_mb_nr];
                    currMB->p_Vid    = slice->p_Vid;
                    currMB->p_Slice  = slice;
                    currMB->slice_nr = slice->current_slice_nr;
                    currMB->mb_field_decoding_flag = slice->mb_data[slice->current_mb_nr-1].mb_field_decoding_flag;
                    currMB->mbAddrX  = slice->current_mb_nr;

                    CheckAvailabilityOfNeighbors(currMB);
                    CheckAvailabilityOfNeighborsCABAC(currMB);

                    //check_next_mb
                    slice->last_dquant = 0;
                    slice->prescan_skip_read = true;
                    slice->prescan_skip_flag = parse_mb_skip_flag(currMB);
                    if (slice->prescan_skip_flag)
                        slice->last_dquant = 0;
                    this->ei_flag = 0;
                    if (!slice->prescan_skip_flag) {
                        slice->prescan_mb_field_decoding_read = true;
                        slice->prescan_mb_field_decoding_flag = parse_mb_field_decoding_flag(currMB);
                        slice->mb_data[slice->current_mb_nr-1].mb_field_decoding_flag = slice->prescan_mb_field_decoding_flag;
                    }

                    //reset
                    slice->current_mb_nr--;

                    CheckAvailabilityOfNeighborsCABAC(currMB);
                }
            } else {
                if (slice->mb_skip_run == 1)
                    this->mb_field_decoding_flag = dp->next_bits(1);
            }
        }

        if (pps->entropy_coding_mode_flag) {
            if (this->mb_skip_flag)
                slice->mb_skip_run = 0;
        } else
            slice->mb_skip_run--;

        moreDataFlag = !this->mb_skip_flag;
    }

    if (moreDataFlag) {
        bool prevMbSkipped = (CurrMbAddr % 2 == 1) ?
            slice->mb_data[CurrMbAddr - 1].mb_skip_flag : 0;
        if (slice->MbaffFrameFlag &&
            (CurrMbAddr % 2 == 0 || (CurrMbAddr % 2 == 1 && prevMbSkipped))) {
            if (slice->prescan_mb_field_decoding_read) {
                slice->prescan_mb_field_decoding_read = false;
                this->mb_field_decoding_flag = slice->prescan_mb_field_decoding_flag;
            } else
                this->mb_field_decoding_flag = parse_mb_field_decoding_flag(this);
        }
    }

    if (pps->entropy_coding_mode_flag)
        CheckAvailabilityOfNeighborsCABAC(this);


    if (moreDataFlag) {
        this->mb_type = parse_mb_type(this);
        this->ei_flag = 0;
    }

    if (slice->slice_type != I_slice && slice->slice_type != SI_slice) {
        if (slice->MbaffFrameFlag && this->mb_field_decoding_flag) {
            slice->num_ref_idx_l0_active_minus1 = ((slice->num_ref_idx_l0_active_minus1 + 1) << 1) - 1;
            slice->num_ref_idx_l1_active_minus1 = ((slice->num_ref_idx_l1_active_minus1 + 1) << 1) - 1;
        }
    }

    this->interpret_mb_mode();

    if (slice->slice_type != B_slice)
        this->noSubMbPartSizeLessThan8x8Flag = 1;
    slice->dec_picture->motion.mb_field_decoding_flag[this->mbAddrX] = this->mb_field_decoding_flag;

    if (this->mb_type == IPCM)
        this->parse_i_pcm();
    else if (this->mb_type == PSKIP || this->mb_type == BSKIP_DIRECT)
        this->parse_skip();
    else if (this->is_intra_block)
        this->parse_intra();
    else
        this->parse_inter();
}


void macroblock_t::parse_i_pcm()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    this->noSubMbPartSizeLessThan8x8Flag = 1;
    this->transform_size_8x8_flag = 0;

    init_macroblock(this);

    if (slice->dp_mode && slice->dpB_NotPresent) {
        for (int y = 0; y < 16; y++) {
            for (int x = 0; x < 16; x++)
                slice->cof[0][y][x] = 1 << (sps->BitDepthY - 1);
        }

        if (sps->chroma_format_idc != YUV400 && !sps->separate_colour_plane_flag) {
            for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
                for (int y = 0; y < sps->MbHeightC; y++) {
                    for (int x = 0; x < sps->MbWidthC; x++)
                        slice->cof[iCbCr + 1][y][x] = 1 << (sps->BitDepthC - 1);
                }
            }
        }
    } else {
        data_partition_t* dp = &slice->partArr[slice->dp_mode ? 1 : 0];

        if (dp->frame_bitoffset & 7)
            dp->f(8 - (dp->frame_bitoffset & 7));

        for (int y = 0; y < 16; y++) {
            for (int x = 0; x < 16; x++)
                slice->cof[0][y][x] = dp->f(sps->BitDepthY);
        }

        if (sps->chroma_format_idc != YUV400 && !sps->separate_colour_plane_flag) {
            for (int iCbCr = 0; iCbCr < 2; iCbCr++) {
                for (int y = 0; y < sps->MbHeightC; y++) {
                    for (int x = 0; x < sps->MbWidthC; x++)
                        slice->cof[iCbCr + 1][y][x] = dp->f(sps->BitDepthC);
                }
            }
        }

        if (pps->entropy_coding_mode_flag)
            dp->de_cabac.init(dp);
    }
}

void macroblock_t::parse_skip()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    pic_motion_params** mv_info = slice->dec_picture->mv_info;
    int slice_no = slice->current_slice_nr;

    this->transform_size_8x8_flag = 0;
    if (pps->constrained_intra_pred_flag)
        this->is_intra_block = 0;

    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            auto mv = &mv_info[this->mb.y * 4 + y][this->mb.x * 4 + x];
            mv->slice_no = slice_no;
            if (slice->slice_type != B_slice) {
                mv->ref_pic[LIST_1] = NULL;
                mv->ref_idx[LIST_1] = -1;
                mv->mv     [LIST_1] = zero_mv;
            }
        }
    }

    if (slice->slice_type == B_slice) {
        this->noSubMbPartSizeLessThan8x8Flag = sps->direct_8x8_inference_flag;

        if (slice->mb_skip_run >= 0) {
            if (pps->entropy_coding_mode_flag) {
                slice->is_reset_coeff = true;
                slice->mb_skip_run = -1;
            } else
                memset(this->nz_coeff, 0, 3 * 16 * sizeof(uint8_t));
        } else {
            this->parse_cbp_qp();
            this->residual();
        }
    } else
        skip_macroblock(this);
}

void macroblock_t::parse_intra()
{
    slice_t* slice = this->p_Slice;
    pps_t* pps = slice->active_pps;

    if (this->mb_type != I4MB)
        this->noSubMbPartSizeLessThan8x8Flag = 1;

    this->transform_size_8x8_flag = 0;
    if (this->mb_type == I4MB && pps->transform_8x8_mode_flag) {
        this->transform_size_8x8_flag = parse_transform_size_8x8_flag(this);
        if (this->transform_size_8x8_flag) {
            this->mb_type = I8MB;
            for (int i = 0; i < 4; i++) {
                this->SubMbType    [i] = I8MB;
                this->SubMbPredMode[i] = -1;
            }
        }
    }

    init_macroblock(this);
    this->parse_ipred_modes();
    this->parse_cbp_qp();

    this->residual();
}

void macroblock_t::parse_inter()
{
    static const char p_v2b8 [ 4] = { 4, 5, 6, 7 };
    static const char p_v2pd [ 4] = { 0, 0, 0, 0 };
    static const char b_v2b8 [13] = { 0, 4, 4, 4, 5, 6, 5, 6, 5, 6, 7, 7, 7 };
    static const char b_v2pd [13] = { 2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 1, 2 };

    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    this->noSubMbPartSizeLessThan8x8Flag = 1;
    this->transform_size_8x8_flag = 0;
    if (this->mb_type == P8x8) {
        for (int mbPartIdx = 0; mbPartIdx < 4; mbPartIdx++) {
            int val = parse_sub_mb_type(this);

            //const uint8_t (*sub_mb_types)[5] = slice->slice_type == B_slice ?
            //    sub_mb_types_b_slice : sub_mb_types_p_slice;

            if (slice->slice_type != B_slice) {
                assert(val < 4);
                this->SubMbType    [mbPartIdx] = p_v2b8[val];
                this->SubMbPredMode[mbPartIdx] = p_v2pd[val];
            } else {
                assert(val < 13);
                this->SubMbType    [mbPartIdx] = b_v2b8[val];
                this->SubMbPredMode[mbPartIdx] = b_v2pd[val];
            }

            //set noSubMbPartSizeLessThan8x8Flag for P8x8 mode
            this->noSubMbPartSizeLessThan8x8Flag &= 
                (this->SubMbType[mbPartIdx] == 0 && sps->direct_8x8_inference_flag) ||
                (this->SubMbType[mbPartIdx] == P8x8);
        }
    }

    if (pps->constrained_intra_pred_flag)
        this->is_intra_block = 0;

    init_macroblock(this);
    this->parse_motion_info();
    this->parse_cbp_qp();

    this->residual();
}


void macroblock_t::parse_ipred_modes()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    if (this->mb_type == I8MB)
        this->parse_ipred_8x8_modes();
    else if (this->mb_type == I4MB)
        this->parse_ipred_4x4_modes();

    if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
        this->intra_chroma_pred_mode = parse_intra_chroma_pred_mode(this);

        assert(this->intra_chroma_pred_mode >= intra_prediction_t::Intra_Chroma_DC ||
               this->intra_chroma_pred_mode <= intra_prediction_t::Intra_Chroma_Plane);
    }
}

void macroblock_t::parse_ipred_4x4_modes()
{
    slice_t* slice = this->p_Slice;
    pps_t* pps = slice->active_pps;

    for (int luma4x4BlkIdx = 0; luma4x4BlkIdx < 16; luma4x4BlkIdx++) {
        int bx = ((luma4x4BlkIdx / 4) % 2) * 8 + ((luma4x4BlkIdx % 4) % 2) * 4;
        int by = ((luma4x4BlkIdx / 4) / 2) * 8 + ((luma4x4BlkIdx % 4) / 2) * 4;
        int val = parse_intra_pred_mode(this);

        bool    prev_intra4x4_pred_mode_flag = val == -1;
        uint8_t rem_intra4x4_pred_mode       = val;

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
        uint8_t intraMxMPredModeA = intra_prediction_t::Intra_4x4_DC;
        uint8_t intraMxMPredModeB = intra_prediction_t::Intra_4x4_DC;
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

        uint8_t predIntra4x4PredMode = min(intraMxMPredModeA, intraMxMPredModeB);
        if (prev_intra4x4_pred_mode_flag)
            this->Intra4x4PredMode[luma4x4BlkIdx] = predIntra4x4PredMode;
        else if (rem_intra4x4_pred_mode < predIntra4x4PredMode)
            this->Intra4x4PredMode[luma4x4BlkIdx] = rem_intra4x4_pred_mode;
        else
            this->Intra4x4PredMode[luma4x4BlkIdx] = rem_intra4x4_pred_mode + 1;

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
        int val = parse_intra_pred_mode(this);

        bool    prev_intra8x8_pred_mode_flag = val == -1;
        uint8_t rem_intra8x8_pred_mode       = val;

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
        uint8_t intraMxMPredModeA = intra_prediction_t::Intra_8x8_DC;
        uint8_t intraMxMPredModeB = intra_prediction_t::Intra_8x8_DC;
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

        uint8_t predIntra8x8PredMode = min(intraMxMPredModeA, intraMxMPredModeB);
        if (prev_intra8x8_pred_mode_flag)
            this->Intra8x8PredMode[luma8x8BlkIdx] = predIntra8x8PredMode;
        else if (rem_intra8x8_pred_mode < predIntra8x8PredMode)
            this->Intra8x8PredMode[luma8x8BlkIdx] = rem_intra8x8_pred_mode;
        else
            this->Intra8x8PredMode[luma8x8BlkIdx] = rem_intra8x8_pred_mode + 1;

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
        slice->inter_prediction.update_direct_mv_info(this);   

    this->parse_ref_pic_idx(LIST_0);
    if (slice->slice_type == B_slice)
        this->parse_ref_pic_idx(LIST_1);

    this->parse_motion_vectors(LIST_0);
    if (slice->slice_type == B_slice)
        this->parse_motion_vectors(LIST_1);

    int list_offset = slice->MbaffFrameFlag && this->mb_field_decoding_flag ?
                      this->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = slice->listX[LIST_0 + list_offset];
    storable_picture **list1 = slice->listX[LIST_1 + list_offset];
    pic_motion_params **p_mv_info = slice->dec_picture->mv_info;
    // record reference picture Ids for deblocking decisions
    for (int j4 = 0; j4 < 4; j4++) {
        for (int i4 = 0; i4 < 4; ++i4) {
            auto mv_info = &p_mv_info[this->mb.y * 4 + j4][this->mb.x * 4 + i4];
            short ref_idx = mv_info->ref_idx[LIST_0];
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
    pic_motion_params **mv_info = slice->dec_picture->mv_info;

    int step_h0 = BLOCK_STEP[this->mb_type][0];
    int step_v0 = BLOCK_STEP[this->mb_type][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {      
            int k = 2 * (j0 >> 1) + (i0 >> 1);

            if ((this->SubMbPredMode[k] == list || this->SubMbPredMode[k] == BI_PRED) &&
                this->SubMbType[k] != 0) {
                uint8_t refframe = parse_ref_idx(this, list, i0, j0);
                for (int j = 0; j < step_v0; j++) {
                    for (int i = 0; i < step_h0; i++)
                        mv_info[this->mb.y * 4 + j0 + j][this->mb.x * 4 + i0 + i].ref_idx[list] = refframe;
                }
            }
        }
    }
}

void macroblock_t::parse_motion_vectors(int list)
{
    slice_t *slice = this->p_Slice;
    pic_motion_params **mv_info = slice->dec_picture->mv_info;

    int step_h0 = BLOCK_STEP[this->mb_type][0];
    int step_v0 = BLOCK_STEP[this->mb_type][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int kk = 2 * (j0 >> 1) + (i0 >> 1);
            int step_h4 = BLOCK_STEP[this->SubMbType[kk]][0];
            int step_v4 = BLOCK_STEP[this->SubMbType[kk]][1];
            char cur_ref_idx = mv_info[this->mb.y * 4 + j0][this->mb.x * 4 + i0].ref_idx[list];

            if ((this->SubMbPredMode[kk] == list || this->SubMbPredMode[kk] == BI_PRED) &&
                (this->SubMbType[kk] != 0)) { //has forward vector
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

    pic_motion_params **mv_info = slice->dec_picture->mv_info;
    //const uint8_t (*sub_mb_types)[5] = this->p_Slice->slice_type == B_slice ?
    //                                   sub_mb_types_b_slice : sub_mb_types_p_slice;
    //int step_h4 = sub_mb_types[this->sub_mb_type[kk]][3];
    //int step_v4 = sub_mb_types[this->sub_mb_type[kk]][4];

    PixelPos block[4]; // neighbor blocks
    MotionVector pred_mv;
    get_neighbors(this, block, BLOCK_SIZE * i, BLOCK_SIZE * j, step_h4);
    GetMVPredictor(this, block, &pred_mv, cur_ref_idx, mv_info, list, BLOCK_SIZE * i, BLOCK_SIZE * j, step_h4, step_v4);

    int16_t curr_mvd[2];
    for (int k = 0; k < 2; ++k)
        curr_mvd[k] = parse_mvd(this, list, i, j, k);

    MotionVector curr_mv;
    curr_mv.mv_x = curr_mvd[0] + pred_mv.mv_x;
    curr_mv.mv_y = curr_mvd[1] + pred_mv.mv_y;

    int i4 = this->mb.x * 4 + i;
    int j4 = this->mb.y * 4 + j;
    for (int jj = 0; jj < step_v4 / 4; ++jj) {
        for (int ii = 0; ii < step_h4 / 4; ++ii) {
            mv_info[jj + j4][ii + i4].mv[list] = curr_mv;
            auto mvd = list == 0 ? this->mvd_l0 : this->mvd_l1;
            mvd[jj + j][ii + i][0] = curr_mvd[0];
            mvd[jj + j][ii + i][1] = curr_mvd[1];
        }
    }
}

void macroblock_t::parse_cbp_qp()
{
    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;
    pps_t *pps = slice->active_pps;

#define IS_I16MB(MB)  ((MB)->mb_type == I16MB || (MB)->mb_type == IPCM)
#define IS_DIRECT(MB) ((MB)->mb_type == 0 && slice->slice_type == B_SLICE)
    // read CBP if not new intra mode
    if (!IS_I16MB(this)) {
        uint8_t coded_block_pattern = parse_coded_block_pattern(this);
        this->CodedBlockPatternLuma   = coded_block_pattern % 16;
        this->CodedBlockPatternChroma = coded_block_pattern / 16;
        if (pps->entropy_coding_mode_flag) {
            if (!this->CodedBlockPatternLuma && !this->CodedBlockPatternChroma)
                slice->last_dquant = 0;
        }

        //============= Transform size flag for INTER MBs =============
        //-------------------------------------------------------------
        int need_transform_size_flag =
            ((this->mb_type >= 1 && this->mb_type <= 3) ||
             (IS_DIRECT(this) && sps->direct_8x8_inference_flag) ||
             this->noSubMbPartSizeLessThan8x8Flag) &&
            (this->mb_type != I8MB) && (this->mb_type != I4MB) &&
            (this->CodedBlockPatternLuma) && pps->transform_8x8_mode_flag;

        if (need_transform_size_flag)
            this->transform_size_8x8_flag = parse_transform_size_8x8_flag(this);
    }

    //=====   DQUANT   =====
    //----------------------
    // Delta quant only if nonzero coeffs
    if (IS_I16MB(this) || this->CodedBlockPatternLuma != 0 || this->CodedBlockPatternChroma != 0) {
        this->mb_qp_delta = parse_mb_qp_delta(this);

        if (pps->entropy_coding_mode_flag)        
            slice->last_dquant = this->mb_qp_delta;

        assert(this->mb_qp_delta >= -(26 + sps->QpBdOffsetY / 2) &&
               this->mb_qp_delta <=  (25 + sps->QpBdOffsetY / 2));

        slice->SliceQpY =
            ((slice->SliceQpY + this->mb_qp_delta + 52 + 2 * sps->QpBdOffsetY)
                % (52 + sps->QpBdOffsetY)) - sps->QpBdOffsetY;

        if (slice->dp_mode) {
            if (!this->is_intra_block && slice->dpC_NotPresent)
                this->dpl_flag = 1;
            if (this->is_intra_block && slice->dpB_NotPresent) {
                this->ei_flag  = 1;
                this->dpl_flag = 1;
            }
            check_dp_neighbors(this);
            if (this->dpl_flag) {
                this->CodedBlockPatternLuma   = 0;
                this->CodedBlockPatternChroma = 0;
            }
        }
    }
#undef IS_DIRECT
#undef IS_I16MB

    this->update_qp(slice->SliceQpY);
}

// Table 8-15 Specification of QPc as a function of qPi

static const uint8_t QP_SCALE_CR[52] = {
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 34, 35,
    35, 36, 36, 37, 37, 37, 38, 38, 38, 39, 39, 39, 39
};

void macroblock_t::update_qp(int qp)
{
    VideoParameters* p_Vid = this->p_Vid;
    sps_t* sps = p_Vid->active_sps;
    pps_t* pps = p_Vid->active_pps;
    int QpOffset[2] = { pps->chroma_qp_index_offset, pps->second_chroma_qp_index_offset };

    this->QpY          = qp;
    this->qp_scaled[0] = qp + sps->QpBdOffsetY;

    for (int i = 0; i < 2; i++) {
        int8_t qpi = clip3(-(sps->QpBdOffsetC), 51, this->QpY + QpOffset[i]);
        int8_t qpc = qpi < 30 ? qpi : QP_SCALE_CR[qpi];

        this->QpC[i]           = qpc;
        this->qp_scaled[i + 1] = qpc + sps->QpBdOffsetC;
    }

    this->TransformBypassModeFlag = (sps->qpprime_y_zero_transform_bypass_flag && this->qp_scaled[0] == 0);
}
