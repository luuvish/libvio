#include "rbsp.h"
#include "sps.h"
#include "vui.h"

void vui_parameters(vui_t *vui, sps_t *sps, rbsp_t *rbsp) {
	vui->aspect_ratio_info_present_flag = rbsp->u(1, "aspect_ratio_info_present_flag");
	if (vui->aspect_ratio_info_present_flag) {
		vui->aspect_ratio_idc = rbsp->u(8, "aspect_ratio_idc");
		if (vui->aspect_ratio_idc == 255) {
			vui->sar_width  = rbsp->u(16, "sar_width");
			vui->sar_height = rbsp->u(16, "sar_height");
		}
	}

	vui->overscan_info_present_flag = rbsp->u(1, "overscan_info_present_flag");
	if (vui->overscan_info_present_flag)
		vui->overscan_appropriate_flag = rbsp->u(1, "overscan_appropriate_flag");

	vui->video_signal_type_present_flag = rbsp->u(1, "video_signal_type_present_flag");
	if (vui->video_signal_type_present_flag) {
		vui->video_format                    = rbsp->u(3, "video_format");
		vui->video_full_range_flag           = rbsp->u(1, "video_full_range_flag");
		vui->colour_description_present_flag = rbsp->u(1, "colour_description_present_flag");
		if (vui->colour_description_present_flag) {
			vui->colour_primaries         = rbsp->u(8, "colour_primaries");
			vui->transfer_characteristics = rbsp->u(8, "transfer_characteristics");
			vui->matrix_coefficients      = rbsp->u(8, "matrix_coefficients");
		}
	}

	vui->chroma_loc_info_present_flag = rbsp->u(1, "chroma_loc_info_present_flag");
	if (vui->chroma_loc_info_present_flag) {
		vui->chroma_sample_loc_type_top_field    = rbsp->ue("chroma_sample_loc_type_top_field");
		vui->chroma_sample_loc_type_bottom_field = rbsp->ue("chroma_sample_loc_type_bottom_field");
	}

	vui->neutral_chroma_indication_flag = rbsp->u(1, "neutral_chroma_indication_flag");
	vui->field_seq_flag                 = rbsp->u(1, "field_seq_flag");

	vui->hrd_parameters_present_flag = rbsp->u(1, "hrd_parameters_present_flag");
	if (vui->hrd_parameters_present_flag) {
		vui->timing_info_present_flag = rbsp->u(1, "timing_info_present_flag");
		if (vui->timing_info_present_flag) {
			vui->num_units_in_tick = rbsp->u(32, "num_units_in_tick");
			vui->time_scale        = rbsp->u(32, "time_scale");
		}
		vui->nal_hrd_parameters_present_flag = rbsp->u(1, "nal_hrd_parameters_present_flag");
		vui->vcl_hrd_parameters_present_flag = rbsp->u(1, "vcl_hrd_parameters_present_flag");
		if (vui->nal_hrd_parameters_present_flag || vui->vcl_hrd_parameters_present_flag) {
			vui->sub_pic_cpb_params_present_flag = rbsp->u(1, "sub_pic_Cpb_params_present_flag");
			if (vui->sub_pic_cpb_params_present_flag) {
				vui->tick_divisor_minus2                = rbsp->u(8, "tick_divisor_minus2");
				vui->du_cpb_removal_delay_length_minus1 = rbsp->u(5, "du_cpb_removal_delay_length_minus1");
			}
			vui->bit_rate_scale                          = rbsp->u(4, "bit_rate_scale");
			vui->cpb_size_scale                          = rbsp->u(4, "cpb_size_scale");
			vui->initial_cpb_removal_delay_length_minus1 = rbsp->u(5, "initial_cpb_removal_delay_length_minus1");
			vui->cpb_removal_delay_length_minus1         = rbsp->u(5, "cpb_removal_delay_length_minus1");
			vui->dpb_output_delay_length_minus1          = rbsp->u(5, "dpb_output_delay_length_minus1");
		}
	}

	for (int i = 0; i <= sps->sps_max_sub_layers_minus1; i++) {
		vui->fixed_pic_rate_flag[i] = rbsp->u(1, "fixed_pic_rate_flag");
		if (vui->fixed_pic_rate_flag[i])
			vui->pic_duration_in_tc_minus1[i] = rbsp->ue("pic_duration_in_tc_minus1");
		vui->low_delay_hrd_flag[i] = rbsp->u(1, "low_delay_hrd_flag");
		vui->cpb_cnt_minus1    [i] = rbsp->ue("cpb_cnt_minus1");
		for (int nal = 0; nal < 2; nal++) {
			if (nal == 0 && vui->nal_hrd_parameters_present_flag ||
				nal == 1 && vui->vcl_hrd_parameters_present_flag) {
				for (int j = 0; j <= vui->cpb_cnt_minus1[i]; j++) {
					vui->bit_size_value_minus1[i][j][nal] = rbsp->ue("bit_size_value_minus1");
					vui->cpb_size_value_minus1[i][j][nal] = rbsp->ue("cpb_size_value_minus1");
					vui->cbr_flag             [i][j][nal] = rbsp->u(1, "cbr_flag");
				}
			}
		}
	}

	vui->bitstream_restriction_flag = rbsp->u(1, "bitstream_restriction_flag");
	if (vui->bitstream_restriction_flag) {
		vui->tiles_fixed_structure_flag              = rbsp->u(1, "tiles_fixed_structure_flag");
		vui->motion_vectors_over_pic_boundaries_flag = rbsp->u(1, "motion_vectors_over_pic_boundaries_flag");
		vui->max_bytes_per_pic_denom                 = rbsp->ue("max_bytes_per_pic_denom");
		vui->max_bits_per_mincu_denom                = rbsp->ue("max_bits_per_mincu_denom");
		vui->log2_max_mv_length_horizontal           = rbsp->ue("log2_max_mv_length_horizontal");
		vui->log2_max_mv_length_vertical             = rbps->ue("log2_max_mv_length_vertical");
	}
}
