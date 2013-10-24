/*
 * =============================================================================
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
 * =============================================================================
 *
 *  File      : interpret.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * =============================================================================
 */

#ifndef __VIO_H264_INTERPRET_H__
#define __VIO_H264_INTERPRET_H__

#include <cstdint>
#include "sets.h"
#include "picture.h"
#include "neighbour.h"
#include "bitstream_cabac.h"


namespace vio  {
namespace h264 {


class Interpreter : public nal_unit_t {
public:
    int         frame_bitoffset;

    Interpreter(uint32_t size=MAX_NAL_UNIT_SIZE);
    Interpreter(const nal_unit_t& nal);

    Interpreter& operator=(const nal_unit_t& nal);

public:
    bool        byte_aligned            (void);
    bool        more_data_in_byte_stream(void);
    bool        more_rbsp_data          (void);
    bool        more_rbsp_trailing_data (void);

    uint32_t    next_bits(uint8_t n);
    uint32_t    read_bits(uint8_t n);

    uint32_t    u (uint8_t n,   const char* name="");
    int32_t     i (uint8_t n,   const char* name="");
    uint32_t    f (uint8_t n,   const char* name="");
    uint32_t    b (uint8_t n=8, const char* name="");

    uint32_t    ue(const char* name="");
    int32_t     se(const char* name="");
    uint32_t    ae(const char* name="");
    uint32_t    ce(const char* name="");
    uint32_t    me(const char* name="");
    uint32_t    te(const char* name="");

public:
    void        seq_parameter_set_rbsp(sps_t& sps);
    void        seq_parameter_set_data(sps_t& sps);
    void        scaling_list(int* scalingList, int sizeOfScalingList,
                             bool& useDefaultScalingMatrixFlag);
    void        seq_parameter_set_extension_rbsp(sps_ext_t& sps_ext);
    void        subset_seq_parameter_set_rbsp   (sub_sps_t& sub_sps);
    void        pic_parameter_set_rbsp          (VideoParameters *p_Vid, pps_t& pps);

    void        sei_rbsp                  (void);
    void        sei_message               (void);
    void        access_unit_delimiter_rbsp(void);
    void        end_of_seq_rbsp           (void);
    void        end_of_stream_rbsp        (void);
    void        filler_data_rbsp          (void);

    void        slice_layer_without_partitioning_rbsp(void);
    void        slice_data_partition_a_layer_rbsp    (void);
    void        slice_data_partition_b_layer_rbsp    (void);
    void        slice_data_partition_c_layer_rbsp    (void);
    void        rbsp_slice_trailing_bits             (void);
    void        rbsp_trailing_bits                   (void);
    void        prefix_nal_unit_rbsp                 (void);
    void        slice_layer_extension_rbsp           (void);

    void        slice_header                 (slice_t& slice);
    void        ref_pic_list_modification    (slice_t& slice);
    void        ref_pic_list_mvc_modification(slice_t& slice);
    void        pred_weight_table            (slice_t& slice);
    void        dec_ref_pic_marking          (slice_t& slice);
    void        slice_data();

    void        vui_parameters(vui_t& vui);
    void        hrd_parameters(hrd_t& hrd);

    void        seq_parameter_set_svc_extension (sps_svc_t& sps_svc);
    void        svc_vui_parameters_extension    (svc_vui_t& svc_vui);
    void        seq_parameter_set_mvc_extension (sps_mvc_t& sps_mvc);
    void        mvc_vui_parameters_extension    (mvc_vui_t& mvc_vui);
    void        seq_parameter_set_mvcd_extension(sps_mvcd_t& sps_mvcd);
};


struct cabac_context_t;

struct cabac_engine_t {
    Interpreter* dp;
    uint16_t    codIRange;
    uint16_t    codIOffset;

    void        init(Interpreter* dp);

    bool        decode_decision (cabac_context_t* ctx);
    bool        decode_bypass   ();
    bool        decode_terminate();
    void        renormD();

    uint32_t    u  (cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx);
    uint32_t    tu (cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx, uint32_t cMax);
    int32_t     ueg(cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx, uint32_t cMax, uint8_t k);
    uint32_t    fl (cabac_context_t* ctx, uint8_t* ctxIdxIncs, uint8_t maxBinIdxCtx, uint32_t cMax);
};


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

    Interpreter partArr[3];
    cabac_engine_t   cabac[3];
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

        Interpreter& cavlc;
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


#endif // __VIO_H264_INTERPRET_H__
