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

#ifndef _VIO_H264_INTER_PREDICTION_H_
#define _VIO_H264_INTER_PREDICTION_H_


namespace vio  {
namespace h264 {


class InterPrediction {
public:
    InterPrediction();
    ~InterPrediction();

    void        init(slice_t& slice);

    void        get_block_luma(storable_picture* curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y,
                    px_t block[16][16], int shift_x, int maxold_x, int maxold_y, ColorPlane pl, mb_t* mb);

    void        perform_mc(mb_t* mb, ColorPlane pl, int pred_dir, int i, int j, int block_size_x, int block_size_y);
    void        motion_compensation(mb_t* mb);

protected:
    void        get_block_chroma(storable_picture* curr_ref, int x_pos, int y_pos,
                    int maxold_x, int maxold_y, int block_size_x, int vert_block_size,
                    px_t block1[16][16], px_t block2[16][16], mb_t* mb);

    void        mc_prediction(px_t* mb_pred,
                    px_t block[16][16], int block_size_y, int block_size_x,
                    mb_t* mb, ColorPlane pl, short l0_refframe, int pred_dir);
    void        bi_prediction(px_t* mb_pred, 
                    px_t block_l0[16][16], px_t block_l1[16][16], int block_size_y, int block_size_x,
                    mb_t* mb, ColorPlane pl, short l0_refframe, short l1_refframe);

    void        check_motion_vector_range(mb_t& mb, const mv_t *mv, slice_t *pSlice);
    int         CheckVertMV(mb_t *currMB, int vec_y, int block_size_y);
};


}
}


#endif // _VIO_H264_INTER_PREDICTION_H_
