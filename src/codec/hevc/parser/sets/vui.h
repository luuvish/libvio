#ifndef __VUI_H__
#define __VUI_H__

#define MAX_SUB_LAYERS ( 8)
#define MAX_CPB        (10)

typedef struct vui_t {
	uint32_t aspect_ratio_info_present_flag;                // u(1)
	uint32_t aspect_ratio_idc;                              // u(8)
	uint32_t sar_width;                                     // u(16)
	uint32_t sar_height;                                    // u(16)

	uint32_t overscan_info_present_flag;                    // u(1)
	uint32_t overscan_appropriate_flag;                     // u(1)

	uint32_t video_signal_type_present_flag;                // u(1)
	uint32_t video_format;                                  // u(3)
	uint32_t video_full_range_flag;                         // u(1)
	uint32_t colour_description_present_flag;               // u(1)
	uint32_t colour_primaries;                              // u(8)
	uint32_t transfer_characteristics;                      // u(8)
	uint32_t matrix_coefficients;                           // u(8)

	uint32_t chroma_loc_info_present_flag;                  // u(1)
	uint32_t chroma_sample_loc_type_top_field;              // ue(v)
	uint32_t chroma_sample_loc_type_bottom_field;           // ue(v)

	uint32_t neutral_chroma_indication_flag;                // u(1)
	uint32_t field_seq_flag;                                // u(1)

	uint32_t hrd_parameters_present_flag;                   // u(1)
	uint32_t timing_info_present_flag;                      // u(1)
	uint32_t num_units_in_tick;                             // u(32)
	uint32_t time_scale;                                    // u(32)
	uint32_t nal_hrd_parameters_present_flag;               // u(1)
	uint32_t vcl_hrd_parameters_present_flag;               // u(1)
	uint32_t sub_pic_cpb_params_present_flag;               // u(1)
	uint32_t tick_divisor_minus2;                           // u(8)
	uint32_t du_cpb_removal_delay_length_minus1;            // u(5)
	uint32_t bit_rate_scale;                                // u(4)
	uint32_t cpb_size_scale;                                // u(4)
	uint32_t initial_cpb_removal_delay_length_minus1;       // u(5)
	uint32_t cpb_removal_delay_length_minus1;               // u(5)
	uint32_t dpb_output_delay_length_minus1;                // u(5)

	uint32_t fixed_pic_rate_flag      [MAX_SUB_LAYERS];     // u(1)
	uint32_t pic_duration_in_tc_minus1[MAX_SUB_LAYERS];     // ue(v)
	uint32_t low_delay_hrd_flag       [MAX_SUB_LAYERS];     // u(1)
	uint32_t cpb_cnt_minus1           [MAX_SUB_LAYERS];     // ue(v)
	uint32_t bit_size_value_minus1    [MAX_SUB_LAYERS][MAX_CPB][2]; // ue(v)
	uint32_t cpb_size_value_minus1    [MAX_SUB_LAYERS][MAX_CPB][2]; // ue(v)
	uint32_t cbr_flag                 [MAX_SUB_LAYERS][MAX_CPB][2]; // u(1)

	uint32_t bitstream_restriction_flag;                    // u(1)
	uint32_t tiles_fixed_structure_flag;                    // u(1)
	uint32_t motion_vectors_over_pic_boundaries_flag;       // u(1)
	uint32_t max_bytes_per_pic_denom;                       // ue(v)
	uint32_t max_bits_per_mincu_denom;                      // ue(v)
	uint32_t log2_max_mv_length_horizontal;                 // ue(v)
	uint32_t log2_max_mv_length_vertical;                   // ue(v)
} vui_t;

#endif /* __VUI_H__ */
