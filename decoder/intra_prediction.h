
/*!
 *************************************************************************************
 * \file mb_prediction.h
 *
 * \brief
 *    Functions for macroblock prediction
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>  
 *************************************************************************************
 */

#ifndef _INTRA_PREDICTION_H_
#define _INTRA_PREDICTION_H_

#ifdef __cplusplus
extern "C" {
#endif

void set_intra_prediction_modes(Slice *currSlice);

int mb_pred_intra4x4  (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int mb_pred_intra16x16(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture);
int mb_pred_intra8x8  (Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture);
int mb_pred_ipcm      (Macroblock *currMB);

#ifdef __cplusplus
}
#endif

#endif /* _INTRA_PREDICTION_H_ */
