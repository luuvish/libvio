#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_

#include "global.h"
#include "dpb.h"
#include "transform.h"

struct slice_t;
struct video_par;
struct inp_par;

enum {
    Intra_4x4 = 0,
    Intra_8x8,
    Intra_16x16,

    Pred_L0 = 0,
    Pred_L1,
    BiPred,
    Direct,

    NA = 0
};

enum {
    I_NxN = 0,
    I_16x16_0_0_0,
    I_16x16_1_0_0,
    I_16x16_2_0_0,
    I_16x16_3_0_0,
    I_16x16_0_1_0,
    I_16x16_1_1_0,
    I_16x16_2_1_0,
    I_16x16_3_1_0,
    I_16x16_0_2_0,
    I_16x16_1_2_0,
    I_16x16_2_2_0,
    I_16x16_3_2_0,
    I_16x16_0_0_1,
    I_16x16_1_0_1,
    I_16x16_2_0_1,
    I_16x16_3_0_1,
    I_16x16_0_1_1,
    I_16x16_1_1_1,
    I_16x16_2_1_1,
    I_16x16_3_1_1,
    I_16x16_0_2_1,
    I_16x16_1_2_1,
    I_16x16_2_2_1,
    I_16x16_3_2_1,
    I_PCM,

    P_L0_16x16 = 0,
    P_L0_L0_16x8,
    P_L0_L0_8x16,
    P_8x8,
    P_8x8ref0,
    P_Skip,

    B_Direct_16x16 = 0,
    B_L0_16x16,
    B_L1_16x16,
    B_Bi_16x16,
    B_L0_L0_16x8,
    B_L0_L0_8x16,
    B_L1_L1_16x8,
    B_L1_L1_8x16,
    B_L0_L1_16x8,
    B_L0_L1_8x16,
    B_L1_L0_16x8,
    B_L1_L0_8x16,
    B_L0_Bi_16x8,
    B_L0_Bi_8x16,
    B_L1_Bi_16x8,
    B_L1_Bi_8x16,
    B_Bi_L0_16x8,
    B_Bi_L0_8x16,
    B_Bi_L1_16x8,
    B_Bi_L1_8x16,
    B_Bi_Bi_16x8,
    B_Bi_Bi_8x16,
    B_8x8,
    B_Skip
};

enum {
    P_L0_8x8 = 0,
    P_L0_8x4,
    P_L0_4x8,
    P_L0_4x4,

    B_Direct_8x8 = 0,
    B_L0_8x8,
    B_L1_8x8,
    B_Bi_8x8,
    B_L0_8x4,
    B_L0_4x8,
    B_L1_8x4,
    B_L1_4x8,
    B_Bi_8x4,
    B_Bi_4x8,
    B_L0_4x4,
    B_L1_4x4,
    B_Bi_4x4
};


enum {
    PSKIP        =  0,
    BSKIP_DIRECT =  0,
    P16x16       =  1,
    P16x8        =  2,
    P8x16        =  3,
    SMB8x8       =  4,
    SMB8x4       =  5,
    SMB4x8       =  6,
    SMB4x4       =  7,
    P8x8         =  8,
    I4MB         =  9,
    I16MB        = 10,
    IBLOCK       = 11,
    SI4MB        = 12,
    I8MB         = 13,
    IPCM         = 14,
    MAXMODE      = 15
};


struct macroblock_t {
    slice_t*    p_Slice;
    video_par*  p_Vid;
    inp_par*    p_Inp;
    int         mbAddrX;
    int         mbAddrA, mbAddrB, mbAddrC, mbAddrD;
    bool        mbAvailA, mbAvailB, mbAvailC, mbAvailD;
    BlockPos    mb;
    int         block_x;
    int         block_y;
    int         pix_x;
    int         pix_y;
    int         pix_c_x;
    int         pix_c_y;

    macroblock_t* mb_up;   //!< pointer to neighboring MB (CABAC)
    macroblock_t* mb_left; //!< pointer to neighboring MB (CABAC)

    macroblock_t* mbup;   // neighbors for loopfilter
    macroblock_t* mbleft; // neighbors for loopfilter

    bool        is_intra_block;
    uint8_t     DeblockCall;

    short       slice_nr;
    char        ei_flag;
    char        dpl_flag;


    bool        mb_skip_flag;
    bool        mb_field_decoding_flag;
    uint8_t     mb_type;
    bool        transform_size_8x8_flag;
    int8_t      mb_qp_delta;

    uint8_t     intra_chroma_pred_mode;
    uint8_t     ref_idx_l0 [4];
    uint8_t     ref_idx_l1 [4];
    int16_t     mvd_l0     [4][4][2];
    int16_t     mvd_l1     [4][4][2];
    uint8_t     sub_mb_type[4];

    bool        noSubMbPartSizeLessThan8x8Flag;
    uint8_t     NumMbPart;
    uint8_t     MbPartPredMode[2];
    uint8_t     MbPartWidth;
    uint8_t     MbPartHeight;
    uint8_t     Intra4x4PredMode[16];
    uint8_t     Intra8x8PredMode[ 4];
    uint8_t     Intra16x16PredMode;
    uint8_t     CodedBlockPatternLuma;
    uint8_t     CodedBlockPatternChroma;
    uint8_t     NumSubMbPart   [4];
    uint8_t     SubMbPredMode  [4];
    uint8_t     SubMbPartWidth [4];
    uint8_t     SubMbPartHeight[4];
    int8_t      QpY;
    int8_t      QpC[2];
    uint8_t     qp_scaled[MAX_PLANE];
    bool        TransformBypassModeFlag;

    uint8_t     nz_coeff[3][4][4]; // cavlc
    uint64_t    cbp_bits[3];       // cabac
    uint64_t    cbp_blks[3];       // deblock

    char        b8mode[4];
    char        b8pdir[4];


    short       DFDisableIdc;
    short       DFAlphaC0Offset;
    short       DFBetaOffset;

    bool        fieldMbInFrameFlag;
    bool        filterInternalEdgesFlag;
    bool        filterLeftMbEdgeFlag;
    bool        filterTopMbEdgeFlag;

    bool        mixedModeEdgeFlag;
    byte        strength_ver[4][4];  // bS
    byte        strength_hor[4][16]; // bS

    void        create(slice_t *slice);
    void        init(slice_t *slice);
    void        parse();
    void        decode();
    bool        close(slice_t *slice);

    void        parse_i_pcm();
    void        parse_skip();
    void        parse_intra();
    void        parse_inter();

    void        parse_ipred_modes();
    void        parse_ipred_4x4_modes();
    void        parse_ipred_8x8_modes();

    void        parse_motion_info();
    void        parse_ref_pic_idx(int list);
    void        parse_motion_vectors(int list);
    void        parse_motion_vector(int list, int step_h4, int step_v4, int i, int j, char cur_ref_idx);

    void        parse_cbp_qp();


    void        interpret_mb_mode();
    void        update_qp(int qp);

    void        residual       ();
    void        residual_luma  (ColorPlane pl);
    void        residual_chroma();
    void        residual_block_cavlc(uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                     ColorPlane pl, bool chroma, bool ac, int blkIdx);
    void        residual_block_cabac(uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                     ColorPlane pl, bool chroma, bool ac, int blkIdx);

    uint8_t     parse_coeff_token(int nC);
    uint8_t     parse_total_zeros(int yuv, int tzVlcIndex);
    uint8_t     parse_run_before(uint8_t zerosLeft);
};

using mb_t = macroblock_t;


#endif
