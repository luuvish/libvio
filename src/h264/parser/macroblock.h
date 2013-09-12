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

//! cbp structure
typedef struct cbp_s {
    int64_t blk;
    int64_t bits;
} CBPStructure;

//! mb_t
struct macroblock_t {
    slice_t*   p_Slice;                    //!< pointer to the current slice
    video_par* p_Vid;                      //!< pointer to VideoParameters
    inp_par*   p_Inp;
    int        mbAddrX;                    //!< current MB address
    int        mbAddrA, mbAddrB, mbAddrC, mbAddrD;
    Boolean    mbAvailA, mbAvailB, mbAvailC, mbAvailD;
    BlockPos   mb;
    int        block_x;
    int        block_y;
    int        pix_x;
    int        pix_y;
    int        pix_c_x;
    int        pix_c_y;

    int        subblock_x;
    int        subblock_y;

    int        qp;                    //!< QP luma
    int        qpc[2];                //!< QP chroma
    int        qp_scaled[MAX_PLANE];  //!< QP scaled for all comps.
    Boolean    is_intra_block;
    Boolean    DeblockCall;

    short      slice_nr;
    char       ei_flag;             //!< error indicator flag that enables concealment
    char       dpl_flag;            //!< error indicator flag that signals a missing data partition

    macroblock_t* mb_up;   //!< pointer to neighboring MB (CABAC)
    macroblock_t* mb_left; //!< pointer to neighboring MB (CABAC)

    macroblock_t* mbup;   // neighbors for loopfilter
    macroblock_t* mbleft; // neighbors for loopfilter

    // some storage of macroblock syntax elements for global access
    bool        mb_skip_flag;
    bool        mb_field_decoding_flag;
    uint8_t     mb_type;
    bool        transform_size_8x8_flag;
    uint8_t     coded_block_pattern;
    int8_t      mb_qp_delta;

    bool        prev_intra4x4_pred_mode_flag[16];
    uint8_t     rem_intra4x4_pred_mode      [16];
    bool        prev_intra8x8_pred_mode_flag[ 4];
    uint8_t     rem_intra8x8_pred_mode      [ 4];
    uint8_t     intra_chroma_pred_mode;

    uint8_t     ref_idx_l0 [4];
    uint8_t     ref_idx_l1 [4];
    uint32_t    mvd_l0     [4][4][2];
    uint32_t    mvd_l1     [4][4][2];
    uint8_t     sub_mb_type[4];

    bool        noSubMbPartSizeLessThan8x8Flag;
    uint8_t     NumMbPart;
    uint8_t     MbPartPredMode[2];
    uint8_t     MbPartWidth;
    uint8_t     MbPartHeight;
    uint8_t     Intra4x4PredMode[16];
    uint8_t     Intra8x8PredMode[4];
    uint8_t     Intra16x16PredMode;
    uint8_t     CodedBlockPatternLuma;
    uint8_t     CodedBlockPatternChroma;
    uint8_t     NumSubMbPart   [4];
    uint8_t     SubMbPredMode  [4];
    uint8_t     SubMbPartWidth [4];
    uint8_t     SubMbPartHeight[4];
    uint8_t     QpY;
    bool        TransformBypassModeFlag;

    uint8_t     nz_coeff[3][4][4];


    short       mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2]; //!< indices correspond to [forw,backw][block_y][block_x][x,y]
    int         cbp;
    CBPStructure s_cbp[3];

    char        b8mode[4];
    char        b8pdir[4];
    char        ipmode_DPCM;


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

    bool        NoMbPartLessThan8x8Flag;

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
