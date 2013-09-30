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
 *  File      : decoder.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _VIO_H264_DECODER_H_
#define _VIO_H264_DECODER_H_


namespace vio  {
namespace h264 {


class Transform;
class IntraPrediction;
class InterPrediction;
class Deblock;


class Decoder {
public:
    Decoder();
    ~Decoder();

    void        init(slice_t& slice);

    void        assign_quant_params(slice_t& slice);

    void        decode(mb_t& mb);

    void        coeff_luma_dc  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_luma_ac  (mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);
    void        coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr);

    void        transform_luma_dc  (mb_t* mb, ColorPlane pl);
    void        transform_chroma_dc(mb_t* mb, ColorPlane pl);

    void        deblock_filter(slice_t& slice);

    // called in erc_do_p.cpp
    void        get_block_luma(storable_picture *curr_ref, int x_pos, int y_pos,
                               int block_size_x, int block_size_y, imgpel block[16][16],
                               int shift_x, int maxold_x, int maxold_y, ColorPlane pl, mb_t* mb);

protected:
    void        decode_one_component(mb_t& mb, ColorPlane curr_plane);
    void        mb_pred_ipcm        (mb_t& mb, ColorPlane curr_plane);
    void        mb_pred_intra       (mb_t& mb, ColorPlane curr_plane);
    void        mb_pred_inter       (mb_t& mb, ColorPlane curr_plane);

public:
    IntraPrediction* intra_prediction;
    InterPrediction* inter_prediction;
    Transform*       transform;
    Deblock*         deblock;
};


}
}


#endif // _VIO_H264_DECODER_H_
