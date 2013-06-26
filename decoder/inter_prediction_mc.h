
/*!
 *************************************************************************************
 * \file mc_prediction.h
 *
 * \brief
 *    definitions for motion compensated prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */

#ifndef _MC_PREDICTION_H_
#define _MC_PREDICTION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "dpb.h"

struct slice_t;
struct macroblock_dec;

void get_block_luma(StorablePicture *curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y, imgpel **block,
                           int shift_x,int maxold_x,int maxold_y, ColorPlane pl, struct macroblock_dec *currMB);

void prepare_direct_params(struct macroblock_dec *currMB, StorablePicture *dec_picture, MotionVector *pmvl0, MotionVector *pmvl1,char *l0_rFrame, char *l1_rFrame);
void perform_mc           (struct macroblock_dec *currMB, ColorPlane pl, StorablePicture *dec_picture, int pred_dir, int i, int j, int block_size_x, int block_size_y);

int get_direct8x8temporal(struct macroblock_dec *currMB, StorablePicture *dec_picture, int block8x8);
int get_direct4x4temporal(struct macroblock_dec *currMB, StorablePicture *dec_picture, int block8x8);
int get_direct8x8spatial_eq(struct macroblock_dec *currMB, StorablePicture *dec_picture, int block8x8, MotionVector *pmvl0, MotionVector *pmvl1, char l0_rFrame, char l1_rFrame);
int get_direct8x8spatial_ne(struct macroblock_dec *currMB, StorablePicture *dec_picture, int block8x8, MotionVector *pmvl0, MotionVector *pmvl1, char l0_rFrame, char l1_rFrame);
int get_direct4x4spatial(struct macroblock_dec *currMB, StorablePicture *dec_picture, int block8x8, MotionVector *pmvl0, MotionVector *pmvl1, char l0_rFrame, char l1_rFrame);
int get_inter8x8(struct macroblock_dec *currMB, StorablePicture *dec_picture, int block8x8);


#ifdef __cplusplus
}
#endif

#endif
