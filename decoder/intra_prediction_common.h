/*!
 *************************************************************************************
 * \file intra16x16_pred.h
 *
 * \brief
 *    definitions for intra 16x16 prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */

#ifndef _INTRA16x16_PRED_H_
#define _INTRA16x16_PRED_H_

#ifdef __cplusplus
extern "C" {
#endif

int intra_pred_16x16_mbaff (Macroblock *currMB, ColorPlane pl, int predmode);
int intra_pred_16x16_normal(Macroblock *currMB, ColorPlane pl, int predmode);
int intra_pred_8x8_normal  (Macroblock *currMB, ColorPlane pl, int ioff, int joff);
int intra_pred_8x8_mbaff   (Macroblock *currMB, ColorPlane pl, int ioff, int joff);
int intra_pred_4x4(Macroblock *currMB, ColorPlane pl, int ioff, int joff, int img_block_x, int img_block_y);

void intra_pred_chroma      (Macroblock *currMB);
void intra_pred_chroma_mbaff(Macroblock *currMB);

#ifdef __cplusplus
}
#endif

#endif
