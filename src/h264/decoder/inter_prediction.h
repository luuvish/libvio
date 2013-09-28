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

    void        fill_wp_params(slice_t *currSlice);

    void        motion_compensation(mb_t* mb);

    void        get_block_luma(storable_picture *curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y, imgpel **block,
                        int shift_x,int maxold_x,int maxold_y, ColorPlane pl, mb_t* mb);

    void        perform_mc(mb_t* mb, ColorPlane pl, int pred_dir, int i, int j, int block_size_x, int block_size_y);

    void        set_chroma_vector(mb_t& mb);

protected:
	void        check_motion_vector_range(mb_t& mb, const mv_t *mv, slice_t *pSlice);
	int         CheckVertMV(mb_t *currMB, int vec_y, int block_size_y);

private:
    char        chroma_vector_adjustment[6][32];
    int         max_mb_vmv_r;

    int   **    tmp_res;
    imgpel**    tmp_block_l0;
    imgpel**    tmp_block_l1;
    imgpel**    tmp_block_l2;
    imgpel**    tmp_block_l3;
};


}
}


#endif // _VIO_H264_INTER_PREDICTION_H_
