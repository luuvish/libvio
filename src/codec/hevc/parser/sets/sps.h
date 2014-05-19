#ifndef __SPS_H__
#define __SPS_H__

#define MAX_VPS_ID                  (16)
#define MAX_SPS_ID                  (32)
#define MAX_SUB_LAYERS              ( 7)
#define MAX_SHORT_TERM_REF_PIC_SETS (64)
#define MAX_LONG_TERM_REF_PICS      (32)

enum {
	CHROMA_FORMAT_IDC_MONO = 0,
	CHROMA_FORMAT_IDC_420  = 1,
	CHROMA_FORMAT_IDC_422  = 2,
	CHROMA_FORMAT_IDC_444  = 3
};

#include "scaling_list.h"
#include "vui.h"

typedef struct seq_parameter_set_t {
	uint32_t video_parameter_set_id                    : 4; // u(4)
	uint32_t sps_max_sub_layers_minus1                 : 3; // u(3)

	ptl_t    ptl;

	uint32_t seq_parameter_set_id                         ; // ue(v)
	uint32_t chroma_format_idc                            ; // ue(v)
	uint32_t separate_colour_plane_flag                : 1; // u(1)

	uint32_t pic_width_in_luma_samples                    ; // ue(v)
	uint32_t pic_height_in_luma_samples                   ; // ue(v)
	uint32_t pic_cropping_flag                         : 1; // u(1)
	uint32_t pic_crop_left_offset                         ; // ue(v)
	uint32_t pic_crop_right_offset                        ; // ue(v)
	uint32_t pic_crop_top_offset                          ; // ue(v)
	uint32_t pic_crop_bottom_offset                       ; // ue(v)
	uint32_t bit_depth_luma_minus8                        ; // ue(v)
	uint32_t bit_depth_chroma_minus8                      ; // ue(v)

	uint32_t log2_max_pic_order_cnt_lsb_minus4            ; // ue(v)
	uint32_t sps_max_dec_pic_buffering[MAX_SUB_LAYERS]    ; // ue(v)
	uint32_t sps_num_reorder_pics     [MAX_SUB_LAYERS]    ; // ue(v)
	uint32_t sps_max_latency_increase [MAX_SUB_LAYERS]    ; // ue(v)

	uint32_t restricted_ref_pic_lists_flag             : 1; // u(1)
	uint32_t lists_modification_present_flag           : 1; // u(1)
	uint32_t log2_min_luma_coding_block_size_minus3       ; // ue(v)
	uint32_t log2_diff_max_min_luma_coding_block_size     ; // ue(v)
	uint32_t log2_min_transform_block_size_minus2         ; // ue(v)
	uint32_t log2_diff_max_min_transform_block_size       ; // ue(v)
	uint32_t max_transform_hierarchy_depth_inter          ; // ue(v)
	uint32_t max_transform_hierarchy_depth_intra          ; // ue(v)

	uint32_t scaling_list_enable_flag                  : 1; // u(1)
	uint32_t sps_scaling_list_data_present_flag        : 1; // u(1)
	scaling_list_t sps_scaling_list;

	uint32_t amp_enabled_flag                          : 1; // u(1)
	uint32_t sample_adaptive_offset_enabled_flag       : 1; // u(1)
	uint32_t pcm_enabled_flag                          : 1; // u(1)
	uint32_t pcm_sample_bit_depth_luma_minus1          : 4; // u(4)
	uint32_t pcm_sample_bit_depth_chroma_minus1        : 4; // u(4)
	uint32_t log2_min_pcm_luma_coding_block_size_minus3   ; // ue(v)
	uint32_t log2_diff_max_min_pcm_luma_coding_block_size ; // ue(v)
	uint32_t pcm_loop_filter_disable_flag              : 1; // u(1)
	uint32_t sps_temporal_id_nesting_flag              : 1; // u(1)

	uint32_t num_short_term_ref_pic_sets                  ; // ue(v)
	rps_t    rps[MAX_SHORT_TERM_REF_PIC_SETS+1];
	uint32_t long_term_ref_pics_present_flag           : 1; // u(1)
	uint32_t num_long_term_ref_pics_sps                   ; // ue(v)
	uint32_t lt_ref_pic_poc_lsb_sps      [MAX_LONG_TERM_REF_PICS]; // u(v)
	uint32_t used_by_curr_pic_lt_sps_flag[MAX_LONG_TERM_REF_PICS]; // u(1)
	uint32_t sps_temporal_mvp_enable_flag              : 1; // u(1)
	uint32_t strong_intra_smoothing_enable_flag        : 1; // u(1)

	uint32_t vui_parameters_present_flag               : 1; // u(1)
	vui_t    vui;
} sps_t;

#endif /* __SPS_H__ */
