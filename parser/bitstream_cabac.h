
/*!
 ***************************************************************************
 * \file
 *    cabac.h
 *
 * \brief
 *    Header file for entropy coding routines
 *
 * \author
 *    Detlev Marpe                                                         \n
 *    Copyright (C) 2000 HEINRICH HERTZ INSTITUTE All Rights Reserved.
 *
 * \date
 *    21. Oct 2000 (Changes by Tobias Oelbaum 28.08.2001)
 ***************************************************************************
 */

#ifndef _CABAC_H_
#define _CABAC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "bitstream.h"

/***********************************************************************
 * D a t a    t y p e s   f o r  C A B A C
 ***********************************************************************
 */

//! struct for context management
typedef struct {
    uint16        state; // index into state-table CP
    unsigned char MPS;   // Least Probable Symbol 0/1 CP
    unsigned char dummy; // for alignment
} BiContextType;

typedef BiContextType *BiContextTypePtr;


/**********************************************************************
 * C O N T E X T S   F O R   T M L   S Y N T A X   E L E M E N T S
 **********************************************************************
 */

#define NUM_MB_TYPE_CTX  11
#define NUM_B8_TYPE_CTX   9
#define NUM_MV_RES_CTX   10
#define NUM_REF_NO_CTX    6
#define NUM_DELTA_QP_CTX  4
#define NUM_MB_AFF_CTX    4

typedef struct motion_info_context_t {
    BiContextType mb_type_contexts [3][NUM_MB_TYPE_CTX];
    BiContextType b8_type_contexts [2][NUM_B8_TYPE_CTX];
    BiContextType mv_res_contexts  [2][NUM_MV_RES_CTX];
    BiContextType ref_no_contexts  [2][NUM_REF_NO_CTX];
    BiContextType delta_qp_contexts   [NUM_DELTA_QP_CTX];
    BiContextType mb_aff_contexts     [NUM_MB_AFF_CTX];
} MotionInfoContexts;

#define NUM_BLOCK_TYPES        22

#define NUM_TRANSFORM_SIZE_CTX  3
#define NUM_IPR_CTX             2
#define NUM_CIPR_CTX            4
#define NUM_CBP_CTX             4
#define NUM_BCBP_CTX            4
#define NUM_MAP_CTX            15
#define NUM_LAST_CTX           15
#define NUM_ONE_CTX             5
#define NUM_ABS_CTX             5

typedef struct texture_info_context_t {
    BiContextType transform_size_contexts   [NUM_TRANSFORM_SIZE_CTX];
    BiContextType ipr_contexts              [NUM_IPR_CTX];
    BiContextType cipr_contexts             [NUM_CIPR_CTX];
    BiContextType cbp_contexts           [3][NUM_CBP_CTX];
    BiContextType bcbp_contexts             [NUM_BLOCK_TYPES][NUM_BCBP_CTX];
    BiContextType map_contexts           [2][NUM_BLOCK_TYPES][NUM_MAP_CTX];
    BiContextType last_contexts          [2][NUM_BLOCK_TYPES][NUM_LAST_CTX];
    BiContextType one_contexts              [NUM_BLOCK_TYPES][NUM_ONE_CTX];
    BiContextType abs_contexts              [NUM_BLOCK_TYPES][NUM_ABS_CTX];
} TextureInfoContexts;

//*********************** end of data type definition for CABAC *******************

// structures that will be declared somewhere else
struct storable_picture;
struct datapartition_dec;
struct syntaxelement_dec;
struct slice_t;
struct macroblock_t;

MotionInfoContexts*  create_contexts_MotionInfo(void);
TextureInfoContexts* create_contexts_TextureInfo(void);
void delete_contexts_MotionInfo(MotionInfoContexts *enco_ctx);
void delete_contexts_TextureInfo(TextureInfoContexts *enco_ctx);

void read_CBP_CABAC                  (struct macroblock_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp);
void readRunLevel_CABAC              (struct macroblock_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp);
void read_dQuant_CABAC               (struct macroblock_t *currMB, SyntaxElement *se, DecodingEnvironment *dep_dp);


int  cabac_startcode_follows(struct slice_t *currSlice, int eos_bit);
int  uvlc_startcode_follows(struct slice_t *currSlice, int dummy);

int  readSyntaxElement_CABAC         (struct macroblock_t *currMB, SyntaxElement *se, DataPartition *this_dataPart);

int  read_and_store_CBP_block_bit(struct macroblock_t *currMB, DecodingEnvironment *dep_dp, int type);


void  init_contexts(struct slice_t *currslice);

#ifdef __cplusplus
}
#endif

#endif /* _CABAC_H_ */
