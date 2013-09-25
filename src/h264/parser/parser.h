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

        SyntaxElement& operator = (const SyntaxElement& se);

        uint32_t    parse_mb_skip_run            ();
        bool        parse_mb_skip_flag           ();
        bool        parse_mb_field_decoding_flag ();
        uint32_t    parse_mb_type                ();
        uint8_t     parse_sub_mb_type            ();

        bool        parse_transform_size_8x8_flag();
        int8_t      parse_intra_pred_mode        ();
        uint8_t     parse_intra_chroma_pred_mode ();
        uint8_t     parse_ref_idx                (uint8_t list, uint8_t x0, uint8_t y0);
        int16_t     parse_mvd                    (uint8_t list, uint8_t x0, uint8_t y0, uint8_t xy);
        uint8_t     parse_coded_block_pattern    ();
        int8_t      parse_mb_qp_delta            ();

    private:
        sps_t&      sps;
        pps_t&      pps;
        slice_t&    slice;
        mb_t&       mb;
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
