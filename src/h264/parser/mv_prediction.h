#ifndef _MV_PREDICTION_H_
#define _MV_PREDICTION_H_

void GetMVPredictor(mb_t *currMB, PixelPos *block, MotionVector *pmv,
                    short ref_frame, pic_motion_params **mv_info,
                    int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y);

#endif
