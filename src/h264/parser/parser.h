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
#include "frame_buffer.h"


namespace vio  {
namespace h264 {


//! Data Partitioning Modes
enum {
    PAR_DP_1,   //!< no data partitioning is supported
    PAR_DP_3    //!< data partitioning with 3 partitions
};

class Parser {
public:
    void        init(slice_t& slice);

    void        parse(sps_t& sps);
    void        parse(pps_t& pps);
    void        parse(slice_t& slice);
    void        parse(mb_t& mb);

public:
    uint32_t    current_mb_nr;

    int         dp_mode;

    data_partition_t partArr[3];
    cabac_contexts_t mot_ctx;

    int         mb_skip_run;

    bool        prescan_skip_read;
    bool        prescan_skip_flag;
    bool        prescan_mb_field_decoding_read;
    bool        prescan_mb_field_decoding_flag;

    int         last_dquant;
    int8_t      QpY;

    bool        is_reset_coeff;
    bool        is_reset_coeff_cr;

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
        uint8_t     ref_idx_l(uint8_t list, uint8_t x0, uint8_t y0);
        int16_t     mvd_l(uint8_t list, uint8_t x0, uint8_t y0, uint8_t xy);
        uint8_t     coded_block_pattern();
        int8_t      mb_qp_delta();

        uint8_t     coeff_token(int nC);
        uint8_t     total_zeros(int yuv, int tzVlcIndex);
        uint8_t     run_before(uint8_t zerosLeft);
        bool        coded_block_flag();
        bool        significant_coeff_flag();
        bool        last_significant_coeff_flag();
        uint16_t    coeff_abs_level_minus1();

    private:
        const sps_t& sps;
        const pps_t& pps;
        slice_t&    slice;
        mb_t&       mb;

        CtxIdxInc ctxidx;

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

    private:
        const sps_t& sps;
        const pps_t& pps;
        slice_t&    slice;
        mb_t&       mb;

        SyntaxElement se;
    };

    class Macroblock {
    public:
        Macroblock(mb_t& mb);
        ~Macroblock();

        void        parse();

        void        mb_type         (uint8_t mb_type);
        void        mb_type_i_slice (uint8_t mb_type);
        void        mb_type_si_slice(uint8_t mb_type);
        void        mb_type_p_slice (uint8_t mb_type);
        void        mb_type_b_slice (uint8_t mb_type);
        void        sub_mb_type();

        void        parse_i_pcm();

        void        mb_pred_intra();
        void        mb_pred_inter();
        void        ref_idx_l(int list);
        void        mvd_l    (int list);

        void        coded_block_pattern();
        void        mb_qp_delta();
        void        parse_cbp_qp();

        mv_t        GetMVPredictor(char ref_frame, int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y);
        void        skip_macroblock();

        int         get_colocated_info (int i, int j);
        void        get_direct_temporal();
        void        get_direct_spatial ();

        void        update_qp(int qp);

        void        erc_dpl();

    private:
        const sps_t& sps;
        const pps_t& pps;
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
