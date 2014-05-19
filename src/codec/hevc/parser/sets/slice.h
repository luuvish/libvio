#ifndef __SLICE_H__
#define __SLICE_H__

#define MAX_ENTRY_POINT_OFFSETS (256)

enum {
	SLICE_TYPE_B = 0,
	SLICE_TYPE_P = 1,
	SLICE_TYPE_I = 2
};

typedef struct slice_header_t {
	uint32_t first_slice_segment_in_pic_flag           : 1; // u(1)
	uint32_t no_output_of_prior_pics_flag              : 1; // u(1)
	uint32_t pic_parameter_set_id                         ; // ue(v)

	uint32_t slice_segment_address                        ; // u(v)
	uint32_t dependent_slice_segment_flag              : 1; // u(1)
	uint32_t slice_type                                   ; // ue(v)
	uint32_t pic_output_flag                           : 1; // u(1)
	uint32_t colour_plane_id                           : 2; // u(2)

	uint32_t pic_order_cnt_lsb                            ; // u(v)
	uint32_t short_term_ref_pic_set_sps_flag           : 1; // u(1)
	uint32_t short_term_ref_pic_set_idx                   ; // u(v)
	uint32_t num_long_term_sps                            ; // ue(v)
	uint32_t num_long_term_pics                           ; // ue(v)
	uint32_t poc_lsb_lt              [MAX_LONG_TERM_PICS] ; // u(v)
	uint32_t used_by_curr_pic_lt_flag[MAX_LONG_TERM_PICS] ; // u(1)

	uint32_t slice_sao_luma_flag                       : 1; // u(1)
	uint32_t slice_sao_chroma_flag                     : 1; // u(1)
	uint32_t pic_temporal_mvp_enable_flag              : 1; // u(1)
	uint32_t num_ref_idx_active_override_flag          : 1; // u(1)
	uint32_t num_ref_idx_l0_active_minus1                 ; // ue(v)
	uint32_t num_ref_idx_l1_active_minus1                 ; // ue(v)
	list_modication_t list_modication;
	uint32_t mvd_l1_zero_flag                          : 1; // u(1)
	uint32_t cabac_init_flag                           : 1; // u(1)

	uint32_t collocated_from_l0_flag                   : 1; // u(1)
	uint32_t collocated_ref_idx                           ; // ue(v)
	pred_weight_t pred_weight;
	uint32_t five_minus_max_num_merge_cand                ; // ue(v)

	int32_t  slice_qp_delta                               ; // se(v)
	int32_t  slice_cb_qp_offset                           ; // se(v)
	int32_t  slice_cr_qp_offset                           ; // se(v)

	uint32_t deblocking_filter_override_flag           : 1; // u(1)
	uint32_t slice_disable_deblocking_filter_flag      : 1; // u(1)
	int32_t  beta_offset_div2                             ; // se(v)
	int32_t  tc_offset_div2                               ; // se(v)
	uint32_t slice_loop_filter_across_slices_enabled_flag : 1; // u(1)

	uint32_t num_entry_point_offsets                      ; // ue(v)
	uint32_t entry_point_offset[MAX_ENTRY_POINT_OFFSETS]  ; // u(v)
} slice_t;

#endif /* __SLICE_H__ */
