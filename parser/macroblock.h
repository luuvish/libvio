/*!
 ************************************************************************
 * \file macroblock.h
 *
 * \brief
 *    Arrays for macroblock encoding
 *
 * \author
 *    Inge Lille-Langoy               <inge.lille-langoy@telenor.com>
 *    Copyright (C) 1999 Telenor Satellite Services, Norway
 ************************************************************************
 */

#ifndef _MACROBLOCK_H_
#define _MACROBLOCK_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "mbuffer.h"
#include "transform.h"

struct slice_t;
struct video_par;
struct inp_par;

//! cbp structure
typedef struct cbp_s {
  int64         blk     ;
  int64         bits    ;
  int64         bits_8x8;
} CBPStructure;

//! Macroblock
typedef struct macroblock_dec {
    struct slice_t        *p_Slice;                    //!< pointer to the current slice
    struct video_par      *p_Vid;                      //!< pointer to VideoParameters
    struct inp_par        *p_Inp;
    int                    mbAddrX;                    //!< current MB address
    int                    mbAddrA, mbAddrB, mbAddrC, mbAddrD;
    Boolean                mbAvailA, mbAvailB, mbAvailC, mbAvailD;
    BlockPos               mb;
    int                    block_x;
    int                    block_y;
    int                    block_y_aff;
    int                    pix_x;
    int                    pix_y;
    int                    pix_c_x;
    int                    pix_c_y;

    int                    subblock_x;
    int                    subblock_y;

    int                    qp;                    //!< QP luma
    int                    qpc[2];                //!< QP chroma
    int                    qp_scaled[MAX_PLANE];  //!< QP scaled for all comps.
    Boolean                is_lossless;
    Boolean                is_intra_block;
    Boolean                is_v_block;
    Boolean                DeblockCall;

    short                  slice_nr;
    char                   ei_flag;             //!< error indicator flag that enables concealment
    char                   dpl_flag;            //!< error indicator flag that signals a missing data partition
    short                  delta_quant;         //!< for rate control
    short                  list_offset;

    struct macroblock_dec *mb_up;   //!< pointer to neighboring MB (CABAC)
    struct macroblock_dec *mb_left; //!< pointer to neighboring MB (CABAC)

    struct macroblock_dec *mbup;   // neighbors for loopfilter
    struct macroblock_dec *mbleft; // neighbors for loopfilter

    // some storage of macroblock syntax elements for global access
    short                  mb_type;
    short                  mvd[2][BLOCK_MULTIPLE][BLOCK_MULTIPLE][2]; //!< indices correspond to [forw,backw][block_y][block_x][x,y]
    int                    cbp;
    CBPStructure           s_cbp[3];

    int                    i16mode;
    char                   b8mode[4];
    char                   b8pdir[4];
    char                   ipmode_DPCM;
    char                   c_ipred_mode;       //!< chroma intra prediction mode
    char                   skip_flag;
    short                  DFDisableIdc;
    short                  DFAlphaC0Offset;
    short                  DFBetaOffset;

    Boolean                mb_field;
    //Flag for MBAFF deblocking;
    byte                   mixedModeEdgeFlag;

    // deblocking strength indices
    byte                   strength_ver[4][4];
    byte                   strength_hor[4][16];

    Boolean                luma_transform_size_8x8_flag;
    Boolean                NoMbPartLessThan8x8Flag;

    void (*itrans_4x4)(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
    void (*itrans_8x8)(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);

    void (*GetMVPredictor)(struct macroblock_dec *currMB, PixelPos *block, 
        MotionVector *pmv, short ref_frame, struct pic_motion_params **mv_info,
        int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y);

    int  (*read_and_store_CBP_block_bit)(struct macroblock_dec *currMB,
        DecodingEnvironmentPtr dep_dp, int type);
    char (*readRefPictureIdx)           (struct macroblock_dec *currMB,
        struct syntaxelement_dec *currSE, struct datapartition_dec *dP,
        char b8mode, int list);

    void (*read_comp_coeff_4x4_CABAC)(struct macroblock_dec *currMB,
        struct syntaxelement_dec *currSE, ColorPlane pl, int (*InvLevelScale4x4)[4], int qp_per, int cbp);
    void (*read_comp_coeff_8x8_CABAC)(struct macroblock_dec *currMB,
        struct syntaxelement_dec *currSE, ColorPlane pl);

    void (*read_comp_coeff_4x4_CAVLC)(struct macroblock_dec *currMB,
        ColorPlane pl, int (*InvLevelScale4x4)[4], int qp_per, int cbp, byte **nzcoeff);
    void (*read_comp_coeff_8x8_CAVLC)(struct macroblock_dec *currMB,
        ColorPlane pl, int (*InvLevelScale8x8)[8], int qp_per, int cbp, byte **nzcoeff);
} Macroblock;

void setup_slice_methods(struct slice_t *currSlice);

void start_macroblock     (struct slice_t *currSlice, Macroblock **currMB);
Boolean  exit_macroblock  (struct slice_t *currSlice, int eos_bit);
void update_qp            (Macroblock *currMB, int qp);


#ifdef __cplusplus
}
#endif

#endif
