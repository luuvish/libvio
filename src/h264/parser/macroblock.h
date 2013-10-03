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
 *  File      : macroblock.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _VIO_H264_MACROBLOCK_H_
#define _VIO_H264_MACROBLOCK_H_


namespace vio  {
namespace h264 {


enum {
    PSKIP        =  0,
    BSKIP_DIRECT =  0,

    P16x16       =  1,
    P16x8        =  2,
    P8x16        =  3,
    P8x8         =  4,
    P8x4         =  5,
    P4x8         =  6,
    P4x4         =  7,
};

enum {
    Intra_16x16 = 0,
    Intra_8x8,
    Intra_4x4,
    Intra_NxN,

    Direct  = 2,
    Pred_L0 = 0,
    Pred_L1,
    BiPred,

    NA = 0
};

enum {
    I_NxN   =  8,
    I_4x4   =  8,
    I_8x8   =  9,
    I_16x16 = 10,
    I_PCM   = 12,
    SI      = 11,

    P_Skip = 0,
    P_16x16,
    P_16x8,
    P_8x16,
    P_8x8     = 4,
    P_8x8ref0 = 12,
    P_8x4     = 5,
    P_4x8,
    P_4x4,

    B_Skip         = 0,
    B_Direct_16x16 = 0,
    B_16x16        = 1,
    B_16x8,
    B_8x16,
    B_Direct_8x8   = 0,
    B_8x8          = 4,
    B_8x4,
    B_4x8,
    B_4x4
};


struct macroblock_t {
    slice_t*    p_Slice;
    int         mbAddrX;

    struct {
        int32_t x;
        int32_t y;
    } mb;

    bool        is_intra_block;

    short       slice_nr;
    bool        ei_flag;
    bool        dpl_flag;

    int         allrefzero;


    bool        mb_skip_flag;
    bool        mb_field_decoding_flag;
    uint8_t     mb_type;
    bool        transform_size_8x8_flag;
    int8_t      mb_qp_delta;

    uint8_t     intra_chroma_pred_mode;
    uint8_t     ref_idx_l0[4];
    uint8_t     ref_idx_l1[4];
    int16_t     mvd_l0    [4][4][2];
    int16_t     mvd_l1    [4][4][2];

    uint8_t     SubMbType    [4];
    uint8_t     SubMbPredMode[4];
    bool        noSubMbPartSizeLessThan8x8Flag;
    uint8_t     Intra4x4PredMode[16];
    uint8_t     Intra8x8PredMode[ 4];
    uint8_t     Intra16x16PredMode;
    uint8_t     CodedBlockPatternLuma;
    uint8_t     CodedBlockPatternChroma;
    int8_t      QpY;
    int8_t      QpC[2];
    int8_t      QsC[2];
    uint8_t     qp_scaled[MAX_PLANE];
    bool        TransformBypassModeFlag;

    uint8_t     nz_coeff[3][4][4]; // cavlc
    uint64_t    cbp_bits[3];       // cabac
    uint64_t    cbp_blks[3];       // deblock

    bool        fieldMbInFrameFlag;
    bool        filterVerEdgeFlag[2][4];
    bool        filterHorEdgeFlag[2][5];
    uint8_t     strength_ver[4][16]; // bS
    uint8_t     strength_hor[5][16]; // bS

    void        create(slice_t& slice);
    void        init(slice_t& slice);
    bool        close(slice_t& slice);
};

using mb_t = macroblock_t;

    
}
}


#endif // _VIO_H264_MACROBLOCK_H_
