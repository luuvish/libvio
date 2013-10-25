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
    int         frame_bitoffset;
private:
    nal_unit_t  nal;

public:
    VideoParameters *p_Vid;
    slice_t*         slice;
};

class InterpreterRbsp : public Interpreter {
public:
    InterpreterRbsp(uint32_t size=MAX_NAL_UNIT_SIZE);
    InterpreterRbsp(const nal_unit_t& nal);

    InterpreterRbsp& operator=(const nal_unit_t& nal);

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

class InterpreterSEI : public InterpreterRbsp {
public:
    InterpreterSEI(const InterpreterRbsp& rbsp);
    ~InterpreterSEI();

public:
    void        sei_payload(uint32_t payloadType, uint32_t payloadSize);

    void        buffering_period                        (uint32_t payloadSize);
    void        pic_timing                              (uint32_t payloadSize);
    void        pan_scan_rect                           (uint32_t payloadSize);
    void        filler_payload                          (uint32_t payloadSize);
    void        user_data_registered_itu_t_t35          (uint32_t payloadSize);
    void        user_data_unregistered                  (uint32_t payloadSize);
    void        recovery_point                          (uint32_t payloadSize);
    void        dec_ref_pic_marking_repetition          (uint32_t payloadSize);
    void        spare_pic                               (uint32_t payloadSize);
    void        scene_info                              (uint32_t payloadSize);
    void        sub_seq_info                            (uint32_t payloadSize);
    void        sub_seq_layer_characteristics           (uint32_t payloadSize);
    void        sub_seq_characteristics                 (uint32_t payloadSize);
    void        full_frame_freeze                       (uint32_t payloadSize);
    void        full_frame_freeze_release               (uint32_t payloadSize);
    void        full_frame_snapshot                     (uint32_t payloadSize);
    void        progressive_refinement_segment_start    (uint32_t payloadSize);
    void        progressive_refinement_segment_end      (uint32_t payloadSize);
    void        motion_constrained_slice_group_set      (uint32_t payloadSize);
    void        film_grain_characteristics              (uint32_t payloadSize);
    void        deblocking_filter_display_preference    (uint32_t payloadSize);
    void        stereo_video_info                       (uint32_t payloadSize);
    void        post_filter_hint                        (uint32_t payloadSize);
    void        tone_mapping_info                       (uint32_t payloadSize);
    void        scalability_info                        (uint32_t payloadSize);
    void        sub_pic_scalable_layer                  (uint32_t payloadSize);
    void        non_required_layer_rep                  (uint32_t payloadSize);
    void        priority_layer_info                     (uint32_t payloadSize);
    void        layers_not_present                      (uint32_t payloadSize);
    void        layer_dependency_change                 (uint32_t payloadSize);
    void        scalable_nesting                        (uint32_t payloadSize);
    void        base_layer_temporal_hrd                 (uint32_t payloadSize);
    void        quality_layer_integrity_check           (uint32_t payloadSize);
    void        redundant_pic_property                  (uint32_t payloadSize);
    void        tl0_dep_rep_index                       (uint32_t payloadSize);
    void        tl_switching_point                      (uint32_t payloadSize);
    void        parallel_decoding_info                  (uint32_t payloadSize);
    void        mvc_scalable_nesting                    (uint32_t payloadSize);
    void        view_scalability_info                   (uint32_t payloadSize);
    void        multiview_scene_info                    (uint32_t payloadSize);
    void        multiview_acquisition_info              (uint32_t payloadSize);
    void        non_required_view_component             (uint32_t payloadSize);
    void        view_dependency_change                  (uint32_t payloadSize);
    void        operation_points_not_present            (uint32_t payloadSize);
    void        base_view_temporal_hrd                  (uint32_t payloadSize);
    void        frame_packing_arrangement               (uint32_t payloadSize);
    void        multiview_view_position                 (uint32_t payloadSize);
    void        display_orientation                     (uint32_t payloadSize);
    void        mvcd_scalable_nesting                   (uint32_t payloadSize);
    void        mvcd_view_scalability_info              (uint32_t payloadSize);
    void        depth_representation_info               (uint32_t payloadSize);
    void        three_dimensional_reference_display_info(uint32_t payloadSize);
    void        depth_timing                            (uint32_t payloadSize);
    void        depth_sampling_info                     (uint32_t payloadSize);
    void        reserved_sei_message                    (uint32_t payloadSize);
};


struct cabac_context_t;

struct cabac_engine_t {
    InterpreterRbsp* dp;
    uint16_t    codIRange;
    uint16_t    codIOffset;

    void        init(InterpreterRbsp* dp);

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

    InterpreterRbsp partArr[3];
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
