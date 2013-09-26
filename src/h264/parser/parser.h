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


#include "bitstream_cabac.h"
#include "data_partition.h"


namespace vio  {
namespace h264 {


class Parser {
public:
    void        parse(sps_t& sps);
    void        parse(pps_t& pps);
    void        parse(slice_t& slice);
    void        parse(mb_t& mb);

protected:
    class SyntaxElement {
    public:
        SyntaxElement(mb_t& mb);
        ~SyntaxElement();

        uint32_t    mb_skip_run();
        bool        mb_skip_flag();
        bool        mb_field_decoding_flag();
        uint8_t     mb_type();
        uint8_t     mb_type_i_slice();
        uint8_t     mb_type_p_slice();
        uint8_t     mb_type_b_slice();
        uint8_t     sub_mb_type();
        uint8_t     sub_mb_type_p_slice();
        uint8_t     sub_mb_type_b_slice();

        bool        transform_size_8x8_flag();
        int8_t      intra_pred_mode();
        uint8_t     intra_chroma_pred_mode();
        uint8_t     ref_idx(uint8_t list, uint8_t x0, uint8_t y0);
        int16_t     mvd(uint8_t list, uint8_t x0, uint8_t y0, uint8_t xy);
        uint8_t     coded_block_pattern();
        int8_t      mb_qp_delta();

    private:
        sps_t&      sps;
        pps_t&      pps;
        slice_t&    slice;
        mb_t&       mb;

        data_partition_t& cavlc;
        cabac_engine_t&   cabac;
        cabac_contexts_t& contexts;
    };

    class Residual {
    public:
        Residual(mb_t& mb);
        ~Residual();

        void        residual       ();
        void        residual_luma  (ColorPlane pl);
        void        residual_chroma();
        void        residual_block_cavlc(uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                         ColorPlane pl, bool chroma, bool ac, int blkIdx);
        void        residual_block_cabac(uint8_t ctxBlockCat, uint8_t startIdx, uint8_t endIdx, uint8_t maxNumCoeff,
                                         ColorPlane pl, bool chroma, bool ac, int blkIdx);

        uint8_t     parse_coeff_token(int nC);
        uint8_t     parse_total_zeros(int yuv, int tzVlcIndex);
        uint8_t     parse_run_before(uint8_t zerosLeft);

    private:
        sps_t&      sps;
        pps_t&      pps;
        slice_t&    slice;
        mb_t&       mb;
    };

    class Macroblock {
    public:
        Macroblock(mb_t& mb);
        ~Macroblock();

        void        parse();
        void        parse_i_pcm();
        void        parse_skip();
        void        parse_intra();
        void        parse_inter();

        void        parse_ipred_modes();
        void        parse_ipred_4x4_modes();
        void        parse_ipred_8x8_modes();

        void        parse_motion_info();
        void        parse_ref_pic_idx(int list);
        void        parse_motion_vectors(int list);
        void        parse_motion_vector(int list, int step_h4, int step_v4, int i, int j, char cur_ref_idx);

        void        parse_cbp_qp();

    private:
        sps_t&      sps;
        pps_t&      pps;
        slice_t&    slice;
        mb_t&       mb;

        SyntaxElement se;
        Residual      re;
    };

private:
};


}
}


#endif // _VIO_H264_PARSER_H_
