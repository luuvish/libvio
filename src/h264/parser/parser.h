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
 *  File      : parser.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _VIO_H264_PARSER_H_
#define _VIO_H264_PARSER_H_


namespace vio {
namespace h264 {


class Parser {
public:
    void        parse(sps_t& sps);
    void        parse(pps_t& pps);
    void        parse(slice_t& slice);
    void        parse(mb_t& mb);

protected:
    void        parse_i_pcm(mb_t& mb);
    void        parse_skip(mb_t& mb);
    void        parse_intra(mb_t& mb);
    void        parse_inter(mb_t& mb);

    void        parse_ipred_modes(mb_t& mb);
    void        parse_ipred_4x4_modes(mb_t& mb);
    void        parse_ipred_8x8_modes(mb_t& mb);

    void        parse_motion_info(mb_t& mb);
    void        parse_ref_pic_idx(mb_t& mb, int list);
    void        parse_motion_vectors(mb_t& mb, int list);
    void        parse_motion_vector(mb_t& mb, int list, int step_h4, int step_v4, int i, int j, char cur_ref_idx);

    void        parse_cbp_qp(mb_t& mb);

    uint32_t    parse_mb_skip_run            (mb_t& mb);
    bool        parse_mb_skip_flag           (mb_t& mb);
    bool        parse_mb_field_decoding_flag (mb_t& mb);
    uint32_t    parse_mb_type                (mb_t& mb);
    uint8_t     parse_sub_mb_type            (mb_t& mb);

    bool        parse_transform_size_8x8_flag(mb_t& mb);
    int8_t      parse_intra_pred_mode        (mb_t& mb);
    uint8_t     parse_intra_chroma_pred_mode (mb_t& mb);
    uint8_t     parse_ref_idx                (mb_t& mb, uint8_t list, uint8_t x0, uint8_t y0);
    int16_t     parse_mvd                    (mb_t& mb, uint8_t list, uint8_t x0, uint8_t y0, uint8_t xy);
    uint8_t     parse_coded_block_pattern    (mb_t& mb);
    int8_t      parse_mb_qp_delta            (mb_t& mb);

    void        residual       (mb_t& mb);
    void        residual_luma  (mb_t& mb, ColorPlane pl);
    void        residual_chroma(mb_t& mb);
    void        residual_block_cavlc(mb_t& mb, uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                     ColorPlane pl, bool chroma, bool ac, int blkIdx);
    void        residual_block_cabac(mb_t& mb, uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                     ColorPlane pl, bool chroma, bool ac, int blkIdx);

    uint8_t     parse_coeff_token(mb_t& mb, int nC);
    uint8_t     parse_total_zeros(mb_t& mb, int yuv, int tzVlcIndex);
    uint8_t     parse_run_before(mb_t& mb, uint8_t zerosLeft);

private:
};


}
}


#endif // _VIO_H264_PARSER_H_
