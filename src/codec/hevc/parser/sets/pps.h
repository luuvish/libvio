#ifndef __PPS_H__
#define __PPS_H__

#define MAX_SPS_ID          ( 32)
#define MAX_PPS_ID          (256)
#define MAX_REF_IDX         ( 16)
#define MAX_TILE_COLUMNS    (1920/8)
#define MAX_TILE_ROWS       (1088/8)

#include "scaling_list.h"

typedef struct pic_parameter_set_t {
	uint32_t pic_parameter_set_id                         ; // ue(v)
	uint32_t seq_parameter_set_id                         ; // ue(v)

	uint32_t sign_data_hiding_flag                     : 1; // u(1)
	uint32_t cabac_init_present_flag                   : 1; // u(1)
	uint32_t num_ref_idx_l0_default_active_minus1         ; // ue(v)
	uint32_t num_ref_idx_l1_default_active_minus1         ; // ue(v)
	int32_t  pic_init_qp_minus26                          ; // se(v)
	uint32_t constrained_intra_pred_flag               : 1; // u(1)
	uint32_t transform_skip_enabled_flag               : 1; // u(1)
	uint32_t cu_qp_delta_enabled_flag                  : 1; // u(1)
	uint32_t diff_cu_qp_delta_depth                       ; // ue(v)
	int32_t  pic_cb_qp_offset                             ; // se(v)
	int32_t  pic_cr_qp_offset                             ; // se(v)
	uint32_t pic_slice_chroma_qp_offsets_present_flag  : 1; // u(1)

	uint32_t weighted_pred_flag                        : 1; // u(1)
	uint32_t weighted_bipred_flag                      : 1; // u(1)
	uint32_t output_flag_present_flag                  : 1; // u(1)
	uint32_t transquant_bypass_enable_flag             : 1; // u(1)
	uint32_t dependent_slice_segments_enabled_flag     : 1; // u(1)

	uint32_t tiles_enabled_flag                        : 1; // u(1)
	uint32_t entropy_coding_sync_enabled_flag          : 1; // u(1)
	uint32_t num_tile_columns_minus1                      ; // ue(v)
	uint32_t num_tile_rows_minus1                         ; // ue(v)
	uint32_t uniform_spacing_flag                      : 1; // u(1)
	uint32_t column_width_minus1[MAX_TILE_COLUMNS]        ; // ue(v)
	uint32_t row_height_minus1  [MAX_TILE_ROWS]           ; // ue(v)
	uint32_t loop_filter_across_tiles_enabled_flag     : 1; // u(1)

	uint32_t loop_filter_across_slices_enabled_flag    : 1; // u(1)
	uint32_t deblocking_filter_control_present_flag    : 1; // u(1)
	uint32_t deblocking_filter_override_enabled_flag   : 1; // u(1)
	uint32_t pic_disable_deblocking_filter_flag        : 1; // u(1)
	int32_t  beta_offset_div2                             ; // se(v)
	int32_t  tc_offset_div2                               ; // se(v)

	uint32_t pic_scaling_list_data_present_flag        : 1; // u(1)
	scaling_list_t pic_scaling_list;

	uint32_t log2_parallel_merge_level_minus2             ; // ue(v)
	uint32_t slice_segment_header_extension_present_flag : 1; // u(1)
} pps_t;

#endif /* __PPS_H__ */
