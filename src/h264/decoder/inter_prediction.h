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
 *  File      : inter_prediction.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _INTER_PREDICTION_H_
#define _INTER_PREDICTION_H_


namespace vio  {
namespace h264 {


struct inter_prediction_t {
    void motion_compensation(mb_t* mb);

    void get_block_luma(storable_picture *curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y, imgpel **block,
                        int shift_x,int maxold_x,int maxold_y, ColorPlane pl, mb_t* mb);

    void perform_mc           (mb_t* mb, ColorPlane pl, int pred_dir, int i, int j, int block_size_x, int block_size_y);

    int  get_colocated_info   (mb_t* mb, storable_picture* list1, int i, int j);
    void set_direct_references(mb_t* mb, const PixelPos* pix, char* l0_rFrame, char* l1_rFrame, pic_motion_params** mv_info);
    void prepare_direct_params(mb_t* mb, MotionVector* pmvl0, MotionVector* pmvl1,char* l0_rFrame, char* l1_rFrame);
    void get_direct_temporal(mb_t* mb, bool dir=true);
    void get_direct_spatial (mb_t* mb, bool dir=true);
    int  get_inter8x8       (mb_t* mb, int block8x8);
};


}
}


#endif /* _INTER_PREDICTION_H_ */
