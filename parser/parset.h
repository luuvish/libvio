
/*!
 **************************************************************************************
 * \file
 *    parset.h
 * \brief
 *    Picture and Sequence Parameter Sets, decoder operations
 * 
 * \date 25 November 2002
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *      - Stephan Wenger        <stewe@cs.tu-berlin.de>
 ***************************************************************************************
 */


#ifndef _PARSET_H_
#define _PARSET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "bitstream_nal.h"

#include "defines.h"
struct video_par;
struct slice_t;


#define MAXIMUMPARSETRBSPSIZE   1500
#define MAXIMUMPARSETNALUSIZE   1500

#define MAXSPS  32
#define MAXPPS  256

#define MAX_NUM_HRD_CPB_CNT 32
#define MAX_NUM_REF_FRAMES 256
#define MAX_NUM_SLICE_GROUPS 8

#define Extended_SAR 255


typedef struct hrd_parameters_t {
    uint8_t         cpb_cnt_minus1;                                   // ue(v)
    uint8_t         bit_rate_scale;                                   // u(4)
    uint8_t         cpb_size_scale;                                   // u(4)
    uint32_t        bit_rate_value_minus1[MAX_NUM_HRD_CPB_CNT];       // ue(v)
    uint32_t        cpb_size_value_minus1[MAX_NUM_HRD_CPB_CNT];       // ue(v)
    bool            cbr_flag             [MAX_NUM_HRD_CPB_CNT];       // u(1)
    uint8_t         initial_cpb_removal_delay_length_minus1;          // u(5)
    uint8_t         cpb_removal_delay_length_minus1;                  // u(5)
    uint8_t         dpb_output_delay_length_minus1;                   // u(5)
    uint8_t         time_offset_length;                               // u(5)
} hrd_t;

typedef struct vui_parameters_t {
    bool            aspect_ratio_info_present_flag;                   // u(1)
    uint8_t         aspect_ratio_idc;                                 // u(8)
    uint16_t        sar_width;                                        // u(16)
    uint16_t        sar_height;                                       // u(16)
    bool            overscan_info_present_flag;                       // u(1)
    bool            overscan_appropriate_flag;                        // u(1)
    bool            video_signal_type_present_flag;                   // u(1)
    uint8_t         video_format;                                     // u(3)
    bool            video_full_range_flag;                            // u(1)
    bool            colour_description_present_flag;                  // u(1)
    uint8_t         colour_primaries;                                 // u(8)
    uint8_t         transfer_characteristics;                         // u(8)
    uint8_t         matrix_coefficients;                              // u(8)
    bool            chroma_loc_info_present_flag;                     // u(1)
    uint8_t         chroma_sample_loc_type_top_field;                 // ue(v)
    uint8_t         chroma_sample_loc_type_bottom_field;              // ue(v)
    bool            timing_info_present_flag;                         // u(1)
    uint32_t        num_units_in_tick;                                // u(32)
    uint32_t        time_scale;                                       // u(32)
    bool            fixed_frame_rate_flag;                            // u(1)
    bool            nal_hrd_parameters_present_flag;                  // u(1)
    hrd_t           nal_hrd_parameters;
    bool            vcl_hrd_parameters_present_flag;                  // u(1)
    hrd_t           vcl_hrd_parameters;
    bool            low_delay_hrd_flag;                               // u(1)
    bool            pic_struct_present_flag;                          // u(1)
    bool            bitstream_restriction_flag;                       // u(1)
    bool            motion_vectors_over_pic_boundaries_flag;          // u(1)
    uint8_t         max_bytes_per_pic_denom;                          // ue(v)
    uint8_t         max_bits_per_mb_denom;                            // ue(v)
    uint8_t         log2_max_mv_length_horizontal;                    // ue(v)
    uint8_t         log2_max_mv_length_vertical;                      // ue(v)
    uint8_t         max_num_reorder_frames;                           // ue(v)
    uint8_t         max_dec_frame_buffering;                          // ue(v)
} vui_t;


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
    bool            Valid;                  // indicates the parameter set is valid

    uint8_t         profile_idc;                                      // u(8)
    bool            constraint_set0_flag;                             // u(1)
    bool            constraint_set1_flag;                             // u(1)
    bool            constraint_set2_flag;                             // u(1)
    bool            constraint_set3_flag;                             // u(1)
    bool            constraint_set4_flag;                             // u(1)
    bool            constraint_set5_flag;                             // u(1)
    uint8_t         level_idc;                                        // u(8)
    uint8_t         seq_parameter_set_id;                             // ue(v)
    uint8_t         chroma_format_idc;                                // ue(v)
    bool            separate_colour_plane_flag;                       // u(1)
    uint8_t         bit_depth_luma_minus8;                            // ue(v)
    uint8_t         bit_depth_chroma_minus8;                          // ue(v)
    bool            qpprime_y_zero_transform_bypass_flag;             // u(1)
    bool            seq_scaling_matrix_present_flag;                  // u(1)
    bool            seq_scaling_list_present_flag[12];                // u(1)

    int             ScalingList4x4[6][16];
    int             ScalingList8x8[6][64];
    bool            UseDefaultScalingMatrix4x4Flag[6];
    bool            UseDefaultScalingMatrix8x8Flag[6];

    uint8_t         log2_max_frame_num_minus4;                        // ue(v)
    uint8_t         pic_order_cnt_type;                               // ue(v)
    uint8_t         log2_max_pic_order_cnt_lsb_minus4;                // ue(v)
    bool            delta_pic_order_always_zero_flag;                 // u(1)
    int32_t         offset_for_non_ref_pic;                           // se(v)
    int32_t         offset_for_top_to_bottom_field;                   // se(v)
    uint8_t         num_ref_frames_in_pic_order_cnt_cycle;            // ue(v)
    int32_t         offset_for_ref_frame[MAX_NUM_REF_FRAMES];         // se(v)
    uint8_t         max_num_ref_frames;                               // ue(v)
    bool            gaps_in_frame_num_value_allowed_flag;             // u(1)
    uint32_t        pic_width_in_mbs_minus1;                          // ue(v)
    uint32_t        pic_height_in_map_units_minus1;                   // ue(v)
    bool            frame_mbs_only_flag;                              // u(1)
    bool            mb_adaptive_frame_field_flag;                     // u(1)
    bool            direct_8x8_inference_flag;                        // u(1)
    bool            frame_cropping_flag;                              // u(1)
    uint32_t        frame_crop_left_offset;                           // ue(v)
    uint32_t        frame_crop_right_offset;                          // ue(v)
    uint32_t        frame_crop_top_offset;                            // ue(v)
    uint32_t        frame_crop_bottom_offset;                         // ue(v)
    bool            vui_parameters_present_flag;                      // u(1)
    vui_t           vui_parameters;

    uint8_t         ChromaArrayType;
    uint8_t         SubWidthC;
    uint8_t         SubHeightC;
    uint8_t         MbWidthC;
    uint8_t         MbHeightC;
    uint8_t         BitDepthY;
    uint8_t         BitDepthC;
    uint8_t         QpBdOffsetY;
    uint8_t         QpBdOffsetC;
    uint16_t        RawMbBits;

    uint32_t        MaxFrameNum;
    uint32_t        MaxPicOrderCntLsb;
    uint32_t        PicWidthInMbs;
    uint32_t        PicWidthInSamplesL;
    uint32_t        PicWidthInSamplesC;
    uint32_t        PicHeightInMapUnits;
    uint32_t        PicSizeInMapUnits;
    uint32_t        FrameHeightInMbs;
    uint8_t         CropUnitX;
    uint8_t         CropUnitY;

#if (MVC_EXTENSION_ENABLE)
    int max_dec_frame_buffering;
#endif
} sps_t;



#if (MVC_EXTENSION_ENABLE)
typedef struct mvcvui_tag {
    int   num_ops_minus1;
    char *temporal_id;
    int  *num_target_output_views_minus1;
    int **view_id;
    char *timing_info_present_flag;
    int  *num_units_in_tick;
    int  *time_scale;
    char *fixed_frame_rate_flag;
    char *nal_hrd_parameters_present_flag;
    char *vcl_hrd_parameters_present_flag;
    char *low_delay_hrd_flag;
    char *pic_struct_present_flag;

    //hrd parameters;
    char cpb_cnt_minus1;
    char bit_rate_scale;
    char cpb_size_scale;
    int  bit_rate_value_minus1[32];
    int  cpb_size_value_minus1[32];
    char cbr_flag[32];
    char initial_cpb_removal_delay_length_minus1;
    char cpb_removal_delay_length_minus1;
    char dpb_output_delay_length_minus1;
    char time_offset_length;
} MVCVUI_t;

typedef struct subset_seq_parameter_set_rbsp_t {
    Boolean      Valid;                  // indicates the parameter set is valid

    sps_t sps;

    unsigned int bit_equal_to_one;
    int          num_views_minus1;
    int         *view_id;
    int         *num_anchor_refs_l0;
    int        **anchor_ref_l0;
    int         *num_anchor_refs_l1;
    int        **anchor_ref_l1;

    int         *num_non_anchor_refs_l0;
    int        **non_anchor_ref_l0;
    int         *num_non_anchor_refs_l1;
    int        **non_anchor_ref_l1;
   
    int          num_level_values_signalled_minus1;
    int         *level_idc;
    int         *num_applicable_ops_minus1;
    int        **applicable_op_temporal_id;
    int        **applicable_op_num_target_views_minus1;
    int       ***applicable_op_target_view_id;
    int        **applicable_op_num_views_minus1;

    unsigned int mvc_vui_parameters_present_flag;
    MVCVUI_t     MVCVUIParams;
} subset_seq_parameter_set_rbsp_t;
#endif



typedef struct pic_parameter_set_t {
    Boolean          Valid;                  // indicates the parameter set is valid

    unsigned int     pic_parameter_set_id;                             // ue(v)
    unsigned int     seq_parameter_set_id;                             // ue(v)
    Boolean          entropy_coding_mode_flag;                         // u(1)
    Boolean          bottom_field_pic_order_in_frame_present_flag;                           // u(1)

    unsigned int     num_slice_groups_minus1;                       // ue(v)
    unsigned int     slice_group_map_type;                          // ue(v)
    unsigned int     run_length_minus1[MAX_NUM_SLICE_GROUPS]; // ue(v)
    unsigned int     top_left         [MAX_NUM_SLICE_GROUPS];          // ue(v)
    unsigned int     bottom_right     [MAX_NUM_SLICE_GROUPS];      // ue(v)
    Boolean          slice_group_change_direction_flag;             // u(1)
    unsigned int     slice_group_change_rate_minus1;                // ue(v)
    unsigned int     pic_size_in_map_units_minus1;                  // ue(v)
    byte            *slice_group_id;                                // complete MBAmap u(v)

    int              num_ref_idx_l0_default_active_minus1;          // ue(v)
    int              num_ref_idx_l1_default_active_minus1;          // ue(v)
    Boolean          weighted_pred_flag;                            // u(1)
    unsigned int     weighted_bipred_idc;                           // u(2)
    int              pic_init_qp_minus26;                           // se(v)
    int              pic_init_qs_minus26;                           // se(v)
    int              chroma_qp_index_offset;                        // se(v)
    Boolean          deblocking_filter_control_present_flag;           // u(1)
    Boolean          constrained_intra_pred_flag;                      // u(1)
    Boolean          redundant_pic_cnt_present_flag;                   // u(1)

    Boolean          transform_8x8_mode_flag;                          // u(1)
    Boolean          pic_scaling_matrix_present_flag;                  // u(1)
    int              pic_scaling_list_present_flag[12];                // u(1)
    int              ScalingList4x4[6][16];                            // se(v)
    int              ScalingList8x8[6][64];                            // se(v)
    bool             UseDefaultScalingMatrix4x4Flag[6];
    bool             UseDefaultScalingMatrix8x8Flag[6];
    int              second_chroma_qp_index_offset;                    // se(v)
} pps_t;


// E.1.1 VUI parameter syntax
void vui_parameters(DataPartition *p, vui_t *vui);

// E.1.2 HRD parameters syntax
void hrd_parameters(DataPartition *p, hrd_t *hrd);

// 7.3.2.1 Sequence parameter set data syntax
void seq_parameter_set_rbsp(DataPartition *p, sps_t *sps);

// 7.3.2.1.1.1 Scaling list syntax
void scaling_list(int *scalingList, int sizeOfScalingList, bool *useDefaultScalingMatrixFlag, Bitstream *s);



pps_t *AllocPPS (void);
sps_t *AllocSPS (void);

void FreePPS (pps_t *pps);
void FreeSPS (sps_t *sps);

void MakePPSavailable (struct video_par *p_Vid, int id, pps_t *pps);

void ProcessSPS (struct video_par *p_Vid, NALU_t *nalu);
void ProcessPPS (struct video_par *p_Vid, NALU_t *nalu);

void CleanUpPPS(struct video_par *p_Vid);

void activate_sps (struct video_par *p_Vid, sps_t *sps);
void activate_pps (struct video_par *p_Vid, pps_t *pps);

void UseParameterSet (struct slice_t *currSlice);

#if (MVC_EXTENSION_ENABLE)
void ProcessSubsetSPS (struct video_par *p_Vid, NALU_t *nalu);
void init_subset_sps_list(subset_seq_parameter_set_rbsp_t *subset_sps_list, int iSize);
void reset_subset_sps(subset_seq_parameter_set_rbsp_t *subset_sps);
int  GetBaseViewId(struct video_par *p_Vid, subset_seq_parameter_set_rbsp_t **subset_sps);
#endif

#ifdef __cplusplus
}
#endif

#endif /* _PARSET_H_ */
