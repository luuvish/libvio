#ifndef _PARSET_H_
#define _PARSET_H_


namespace vio { namespace h264 {
struct nal_unit_t;
struct data_partition_t;
}}
struct VideoParameters;
struct slice_t;


using vio::h264::nal_unit_t;
using vio::h264::data_partition_t;


#define MAX_NUM_SPS  32
#define MAX_NUM_PPS 256

enum {
    CHROMA_FORMAT_400 = 0,
    CHROMA_FORMAT_420 = 1,
    CHROMA_FORMAT_422 = 2,
    CHROMA_FORMAT_444 = 3
};

#define MAX_NUM_HRD_CPB_CNT 32

#define MAX_NUM_REF_FRAMES_IN_POC_CYCLE 256
#define MAX_NUM_SLICE_GROUPS 8


typedef struct hrd_parameters_t {
    uint8_t     cpb_cnt_minus1;                                       // ue(v)
    uint8_t     bit_rate_scale;                                       // u(4)
    uint8_t     cpb_size_scale;                                       // u(4)
    uint32_t    bit_rate_value_minus1[MAX_NUM_HRD_CPB_CNT];           // ue(v)
    uint32_t    cpb_size_value_minus1[MAX_NUM_HRD_CPB_CNT];           // ue(v)
    bool        cbr_flag             [MAX_NUM_HRD_CPB_CNT];           // u(1)
    uint8_t     initial_cpb_removal_delay_length_minus1;              // u(5)
    uint8_t     cpb_removal_delay_length_minus1;                      // u(5)
    uint8_t     dpb_output_delay_length_minus1;                       // u(5)
    uint8_t     time_offset_length;                                   // u(5)
} hrd_t;

typedef struct vui_parameters_t {
    bool        aspect_ratio_info_present_flag;                       // u(1)
    uint8_t     aspect_ratio_idc;                                     // u(8)
    uint16_t    sar_width;                                            // u(16)
    uint16_t    sar_height;                                           // u(16)
    bool        overscan_info_present_flag;                           // u(1)
    bool        overscan_appropriate_flag;                            // u(1)
    bool        video_signal_type_present_flag;                       // u(1)
    uint8_t     video_format;                                         // u(3)
    bool        video_full_range_flag;                                // u(1)
    bool        colour_description_present_flag;                      // u(1)
    uint8_t     colour_primaries;                                     // u(8)
    uint8_t     transfer_characteristics;                             // u(8)
    uint8_t     matrix_coefficients;                                  // u(8)
    bool        chroma_loc_info_present_flag;                         // u(1)
    uint8_t     chroma_sample_loc_type_top_field;                     // ue(v)
    uint8_t     chroma_sample_loc_type_bottom_field;                  // ue(v)
    bool        timing_info_present_flag;                             // u(1)
    uint32_t    num_units_in_tick;                                    // u(32)
    uint32_t    time_scale;                                           // u(32)
    bool        fixed_frame_rate_flag;                                // u(1)
    bool        nal_hrd_parameters_present_flag;                      // u(1)
    hrd_t       nal_hrd_parameters;
    bool        vcl_hrd_parameters_present_flag;                      // u(1)
    hrd_t       vcl_hrd_parameters;
    bool        low_delay_hrd_flag;                                   // u(1)
    bool        pic_struct_present_flag;                              // u(1)
    bool        bitstream_restriction_flag;                           // u(1)
    bool        motion_vectors_over_pic_boundaries_flag;              // u(1)
    uint8_t     max_bytes_per_pic_denom;                              // ue(v)
    uint8_t     max_bits_per_mb_denom;                                // ue(v)
    uint8_t     log2_max_mv_length_horizontal;                        // ue(v)
    uint8_t     log2_max_mv_length_vertical;                          // ue(v)
    uint8_t     max_num_reorder_frames;                               // ue(v)
    uint8_t     max_dec_frame_buffering;                              // ue(v)
} vui_t;



//AVC Profile IDC definitions
enum {
    FREXT_CAVLC444 = 44,       //!< YUV 4:4:4/14 "CAVLC 4:4:4"
    BASELINE       = 66,       //!< YUV 4:2:0/8  "Baseline"
    MAIN           = 77,       //!< YUV 4:2:0/8  "Main"
    EXTENDED       = 88,       //!< YUV 4:2:0/8  "Extended"
    FREXT_HP       = 100,      //!< YUV 4:2:0/8  "High"
    FREXT_Hi10P    = 110,      //!< YUV 4:2:0/10 "High 10"
    FREXT_Hi422    = 122,      //!< YUV 4:2:2/10 "High 4:2:2"
    FREXT_Hi444    = 244,      //!< YUV 4:4:4/14 "High 4:4:4"
    MVC_HIGH       = 118,      //!< YUV 4:2:0/8  "Multiview High"
    STEREO_HIGH    = 128       //!< YUV 4:2:0/8  "Stereo High"
};

// A.2 Profiles
enum {
    CAVLC444_Profile =  44,
    Baseline_Profile =  66,
    Main_Profile     =  77,
    Extended_Profile =  88,
    High_Profile     = 100,
    High10_Profile   = 110,
    High422_Profile  = 122,
    High444_Profile  = 244
};

typedef struct seq_parameter_set_t {
    bool        Valid;                  // indicates the parameter set is valid

    uint8_t     profile_idc;                                          // u(8)
    bool        constraint_set0_flag;                                 // u(1)
    bool        constraint_set1_flag;                                 // u(1)
    bool        constraint_set2_flag;                                 // u(1)
    bool        constraint_set3_flag;                                 // u(1)
    bool        constraint_set4_flag;                                 // u(1)
    bool        constraint_set5_flag;                                 // u(1)
    uint8_t     level_idc;                                            // u(8)
    uint8_t     seq_parameter_set_id;                                 // ue(v)
    uint8_t     chroma_format_idc;                                    // ue(v)
    bool        separate_colour_plane_flag;                           // u(1)
    uint8_t     bit_depth_luma_minus8;                                // ue(v)
    uint8_t     bit_depth_chroma_minus8;                              // ue(v)
    bool        qpprime_y_zero_transform_bypass_flag;                 // u(1)
    bool        seq_scaling_matrix_present_flag;                      // u(1)
    bool        seq_scaling_list_present_flag[12];                    // u(1)
    uint8_t     log2_max_frame_num_minus4;                            // ue(v)
    uint8_t     pic_order_cnt_type;                                   // ue(v)
    uint8_t     log2_max_pic_order_cnt_lsb_minus4;                    // ue(v)
    bool        delta_pic_order_always_zero_flag;                     // u(1)
    int32_t     offset_for_non_ref_pic;                               // se(v)
    int32_t     offset_for_top_to_bottom_field;                       // se(v)
    uint8_t     num_ref_frames_in_pic_order_cnt_cycle;                // ue(v)
    int32_t     offset_for_ref_frame[MAX_NUM_REF_FRAMES_IN_POC_CYCLE];// se(v)
    uint8_t     max_num_ref_frames;                                   // ue(v)
    bool        gaps_in_frame_num_value_allowed_flag;                 // u(1)
    uint32_t    pic_width_in_mbs_minus1;                              // ue(v)
    uint32_t    pic_height_in_map_units_minus1;                       // ue(v)
    bool        frame_mbs_only_flag;                                  // u(1)
    bool        mb_adaptive_frame_field_flag;                         // u(1)
    bool        direct_8x8_inference_flag;                            // u(1)
    bool        frame_cropping_flag;                                  // u(1)
    uint32_t    frame_crop_left_offset;                               // ue(v)
    uint32_t    frame_crop_right_offset;                              // ue(v)
    uint32_t    frame_crop_top_offset;                                // ue(v)
    uint32_t    frame_crop_bottom_offset;                             // ue(v)
    bool        vui_parameters_present_flag;                          // u(1)
    vui_t       vui_parameters;

    uint8_t     ChromaArrayType;
    uint8_t     SubWidthC;
    uint8_t     SubHeightC;
    uint8_t     MbWidthC;
    uint8_t     MbHeightC;
    uint8_t     BitDepthY;
    uint8_t     BitDepthC;
    uint8_t     QpBdOffsetY;
    uint8_t     QpBdOffsetC;
    uint16_t    RawMbBits;
    int         ScalingList4x4[6][16];
    int         ScalingList8x8[6][64];
    bool        UseDefaultScalingMatrix4x4Flag[6];
    bool        UseDefaultScalingMatrix8x8Flag[6];
    uint32_t    MaxFrameNum;
    uint32_t    MaxPicOrderCntLsb;
    int32_t     ExpectedDeltaPerPicOrderCntCycle;
    uint32_t    PicWidthInMbs;
    uint32_t    PicWidthInSamplesL;
    uint32_t    PicWidthInSamplesC;
    uint32_t    PicHeightInMapUnits;
    uint32_t    PicSizeInMapUnits;
    uint32_t    FrameHeightInMbs;
    uint8_t     CropUnitX;
    uint8_t     CropUnitY;
    uint8_t     MaxDpbFrames;
} sps_t;

typedef struct seq_parameter_set_extension_t {
    uint8_t     seq_parameter_set_id;                                 // ue(v)
    uint8_t     aux_format_idc;                                       // ue(v)
    uint8_t     bit_depth_aux_minus8;                                 // ue(v)
    bool        alpha_incr_flag;                                      // u(1)
    uint16_t    alpha_opaque_value;                                   // u(v)
    uint16_t    alpha_transparent_value;                              // u(v)
    bool        additional_extension_flag;                            // u(1)
} sps_ext_t;

bool operator==(const sps_t& l, const sps_t& r);


#if (MVC_EXTENSION_ENABLE)
typedef struct seq_parameter_set_svc_extension_t {
    bool        inter_layer_deblocking_filter_control_present_flag;   // u(1)
    uint8_t     extended_spatial_scalability_idc;                     // u(2)
    bool        chroma_phase_x_plus1_flag;                            // u(1)
    uint8_t     chroma_phase_y_plus1;                                 // u(2)
    bool        seq_ref_layer_chroma_phase_x_plus1_flag;              // u(1)
    uint8_t     seq_ref_layer_chroma_phase_y_plus1;                   // u(2)
    int16_t     seq_scaled_ref_layer_left_offset;                     // se(v)
    int16_t     seq_scaled_ref_layer_top_offset;                      // se(v)
    int16_t     seq_scaled_ref_layer_right_offset;                    // se(v)
    int16_t     seq_scaled_ref_layer_bottom_offset;                   // se(v)
    bool        seq_tcoeff_level_prediction_flag;                     // u(1)
    bool        adaptive_tcoeff_level_prediction_flag;                // u(1)
    bool        slice_header_restriction_flag;                        // u(1)
} sps_svc_t;

typedef struct svc_vui_parameters_extension_t {
    uint8_t     vui_ext_num_entries_minus1;                           // ue(v)
    struct {
        uint8_t     vui_ext_dependency_id;                            // u(3)
        uint8_t     vui_ext_quality_id;                               // u(4)
        uint8_t     vui_ext_temporal_id;                              // u(3)
        bool        vui_ext_timing_info_present_flag;                 // u(1)
        uint32_t    vui_ext_num_units_in_tick;                        // u(32)
        uint32_t    vui_ext_time_scale;                               // u(32)
        bool        vui_ext_fixed_frame_rate_flag;                    // u(1)
        bool        vui_ext_nal_hrd_parameters_present_flag;          // u(1)
        hrd_t       vui_ext_nal_parameters;
        bool        vui_ext_vcl_hrd_parameters_present_flag;          // u(1)
        hrd_t       vui_ext_vcl_parameters;
        bool        vui_ext_low_delay_hrd_flag;                       // u(1)
        bool        vui_ext_pic_struct_present_flag;                  // u(1)
    } entries[1024];
} svc_vui_t;

typedef struct seq_parameter_set_mvc_extension_t {
    uint32_t    num_views_minus1;                                     // ue(v)
    struct {
        uint16_t    view_id;                                          // ue(v)
        uint8_t     num_anchor_refs_l0;                               // ue(v)
        uint16_t    anchor_ref_l0[15];                                // ue(v)
        uint8_t     num_anchor_refs_l1;                               // ue(v)
        uint16_t    anchor_ref_l1[15];                                // ue(v)
        uint8_t     num_non_anchor_refs_l0;                           // ue(v)
        uint16_t    non_anchor_ref_l0[15];                            // ue(v)
        uint8_t     num_non_anchor_refs_l1;                           // ue(v)
        uint16_t    non_anchor_ref_l1[15];                            // ue(v)
    } views[1024];
    uint32_t    num_level_values_signalled_minus1;                    // ue(v)
    struct {
        uint8_t     level_idc;                                        // u(8)
        uint16_t    num_applicable_ops_minus1;                        // ue(v)
        struct {
            uint8_t     applicable_op_temporal_id;                    // u(3)
            uint16_t    applicable_op_num_target_views_minus1;        // ue(v)
            uint16_t    applicable_op_target_view_id[1024];           // ue(v)
            uint16_t    applicable_op_num_views_minus1;               // ue(v)
        } applicable_ops[1024];
    } level_values_signalled[64];
} sps_mvc_t;

typedef struct mvc_vui_parameters_t {
    int         num_ops_minus1;
    char *      temporal_id;
    int  *      num_target_output_views_minus1;
    int **      view_id;
    char *      timing_info_present_flag;
    int  *      num_units_in_tick;
    int  *      time_scale;
    char *      fixed_frame_rate_flag;
    char *      nal_hrd_parameters_present_flag;
    char *      vcl_hrd_parameters_present_flag;
    char *      low_delay_hrd_flag;
    char *      pic_struct_present_flag;

    //hrd parameters;
    char        cpb_cnt_minus1;
    char        bit_rate_scale;
    char        cpb_size_scale;
    int         bit_rate_value_minus1[32];
    int         cpb_size_value_minus1[32];
    char        cbr_flag[32];
    char        initial_cpb_removal_delay_length_minus1;
    char        cpb_removal_delay_length_minus1;
    char        dpb_output_delay_length_minus1;
    char        time_offset_length;
} mvc_vui_t;

typedef struct seq_parameter_set_mvcd_t {

} sps_mvcd_t;

typedef struct subset_seq_parameter_set_t {
    bool        Valid;

    sps_t       sps;

    //sps_svc_t   sps_svc;
    bool        svc_vui_parameters_present_flag;                      // u(1)
    //svc_vui_t   svc_vui_parameters;

    bool        bit_equal_to_one;                                     // f(1)

    //sps_mvc_t   sps_mvc;
    bool        mvc_vui_parameters_present_flag;                      // u(1)
    //mvc_vui_t   mvc_vui_parameters;

    sps_mvcd_t  sps_mvcd;


// seq_parameter_set_mvc_extension_t
    int         num_views_minus1;
    int *       view_id;
    int *       num_anchor_refs_l0;
    int**       anchor_ref_l0;
    int *       num_anchor_refs_l1;
    int**       anchor_ref_l1;

    int *       num_non_anchor_refs_l0;
    int**       non_anchor_ref_l0;
    int *       num_non_anchor_refs_l1;
    int**       non_anchor_ref_l1;
   
    int         num_level_values_signalled_minus1;
    int  *      level_idc;
    int  *      num_applicable_ops_minus1;
    int **      applicable_op_temporal_id;
    int **      applicable_op_num_target_views_minus1;
    int***      applicable_op_target_view_id;
    int **      applicable_op_num_views_minus1;

    mvc_vui_t   MVCVUIParams;
} sub_sps_t;
#endif


typedef struct pic_parameter_set_t {
    bool        Valid;

    uint8_t     pic_parameter_set_id;                                 // ue(v)
    uint8_t     seq_parameter_set_id;                                 // ue(v)
    bool        entropy_coding_mode_flag;                             // u(1)
    bool        bottom_field_pic_order_in_frame_present_flag;         // u(1)
    uint8_t     num_slice_groups_minus1;                              // ue(v)
    uint8_t     slice_group_map_type;                                 // ue(v)
    uint32_t    run_length_minus1[MAX_NUM_SLICE_GROUPS];              // ue(v)
    uint32_t    top_left         [MAX_NUM_SLICE_GROUPS];              // ue(v)
    uint32_t    bottom_right     [MAX_NUM_SLICE_GROUPS];              // ue(v)
    bool        slice_group_change_direction_flag;                    // u(1)
    uint32_t    slice_group_change_rate_minus1;                       // ue(v)
    uint32_t    pic_size_in_map_units_minus1;                         // ue(v)
    uint8_t    *slice_group_id;                                       // u(v)
    uint8_t     num_ref_idx_l0_default_active_minus1;                 // ue(v)
    uint8_t     num_ref_idx_l1_default_active_minus1;                 // ue(v)
    bool        weighted_pred_flag;                                   // u(1)
    uint8_t     weighted_bipred_idc;                                  // u(2)
    int8_t      pic_init_qp_minus26;                                  // se(v)
    int8_t      pic_init_qs_minus26;                                  // se(v)
    int8_t      chroma_qp_index_offset;                               // se(v)
    bool        deblocking_filter_control_present_flag;               // u(1)
    bool        constrained_intra_pred_flag;                          // u(1)
    bool        redundant_pic_cnt_present_flag;                       // u(1)
    bool        transform_8x8_mode_flag;                              // u(1)
    bool        pic_scaling_matrix_present_flag;                      // u(1)
    bool        pic_scaling_list_present_flag[12];                    // u(1)
    int8_t      second_chroma_qp_index_offset;                        // se(v)

    uint32_t    SliceGroupChangeRate;
    int         ScalingList4x4[6][16];
    int         ScalingList8x8[6][64];
    bool        UseDefaultScalingMatrix4x4Flag[6];
    bool        UseDefaultScalingMatrix8x8Flag[6];

    pic_parameter_set_t();
    ~pic_parameter_set_t();
} pps_t;

bool operator==(const pps_t& l, const pps_t& r);




// 7.3.2.2 Picture parameter set RBSP syntax
void pic_parameter_set_rbsp(VideoParameters* p_Vid, data_partition_t *p, pps_t *pps);

// 7.3.2.5 End of sequence RBSP syntax
void end_of_seq_rbsp();
// 7.3.2.6 End of stream RBSP syntax
void end_of_stream_rbsp();


#if (MVC_EXTENSION_ENABLE)
typedef struct nalunitheadermvcext_tag {
    unsigned int non_idr_flag;
    unsigned int priority_id;
    unsigned int view_id;
    unsigned int temporal_id;
    unsigned int anchor_pic_flag;
    unsigned int inter_view_flag;
    unsigned int reserved_one_bit;
    unsigned int iPrefixNALU;
} NALUnitHeaderMVCExt_t;
#endif


#if (MVC_EXTENSION_ENABLE)
struct nalunitheadermvcext_tag;
void nal_unit_header_mvc_extension(struct nalunitheadermvcext_tag *NaluHeaderMVCExt, struct data_partition_t *bitstream);
void nal_unit_header_svc_extension();
void prefix_nal_unit_svc();
#endif


void MakePPSavailable (VideoParameters* p_Vid, int id, pps_t *pps);

void ProcessPPS (VideoParameters* p_Vid, nal_unit_t *nalu);

void CleanUpPPS(VideoParameters* p_Vid);

void activate_sps (VideoParameters* p_Vid, sps_t *sps);
void activate_pps (VideoParameters* p_Vid, pps_t *pps);

void UseParameterSet (struct slice_t *currSlice);

#if (MVC_EXTENSION_ENABLE)
void ProcessSubsetSPS (VideoParameters* p_Vid, nal_unit_t *nalu);
void init_subset_sps_list(sub_sps_t *subset_sps_list, int iSize);
void reset_subset_sps(sub_sps_t *subset_sps);
#endif

#endif /* _PARSET_H_ */
