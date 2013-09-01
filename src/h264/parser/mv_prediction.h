/*!
 *************************************************************************************
 * \file mv_prediction.h
 *
 * \brief
 *    Declarations for Motion Vector Prediction
 *
 * \author
 *     Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Alexis Michael Tourapis         <alexismt@ieee.org>  
 *************************************************************************************
 */

#ifndef _MV_PREDICTION_H_
#define _MV_PREDICTION_H_

void GetMVPredictor(mb_t *currMB, PixelPos *block, MotionVector *pmv,
                    short ref_frame, PicMotionParams **mv_info,
                    int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y);

#endif
