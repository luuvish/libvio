#include "rbsp.h"
#include "slice.h"

enum {
	COLOUR_PLANE_Y   = 0,
	COLOUR_PLANE_CB  = 1,
	COLOUR_PLANE_CR  = 2,
	COLOUR_PLANE_NUM = 3
};

enum {
	SAO_TYPE_IDX_NONE_APPLIED = 0,
	SAO_TYPE_IDX_1D_0         = 1,
	SAO_TYPE_IDX_1D_90        = 2,
	SAO_TYPE_IDX_1D_135       = 3,
	SAO_TYPE_IDX_1D_45        = 4,
	SAO_TYPE_IDX_BAND_OFFSET  = 5
};

void slice_segment_layer_rbps(slice_t *slice, rbsp_t *rbsp) {
	slice_segment_header(slice, rbsp);
	slice_segment_data(slice, rbsp);

	rbsp_slice_segment_trailing_bits(rbsp);
}

void slice_segment_header(slice_t *slice, rbsp_t *rbsp) {
	slice->first_slice_segment_in_pic_flag = rbsp->u(1, "first_slice_in_pic_flag");
	if (nal->RapPicFlag)
		slice->no_output_of_prior_pics_flag = rbsp->u(1, "no_output_of_prior_pics_flag");
	slice->pic_parameter_set_id = rbsp->ue("pic_parameter_set_id");
	assert(slice->pic_parameter_set_id < MAX_PPS_ID);

	pps_t *pps = slice->pic_parameter_set_id;
	sps_t *sps = pps->seq_parameter_set_id;

	slice->slice_segment_address = 0;
	if (! slice->first_slice_segment_in_pic_flag) {
		slice->slice_segment_address = rbsp->u(ceil(log2(sps->PicSizeInCtbsY)), "slice_address");
		assert(slice->slice_segment_address > 0);
		assert(slice->slice_segment_address < sps->PicSizeInCtbsY);
	}

	slice->dependent_slice_segment_flag = 0;
	if (pps->dependent_slice_segments_enabled_flag && ! slice->first_slice_segment_in_pic_flag)
		slice->dependent_slice_segment_flag = rbsp->u(1, "dependent_slice_flag");

	if (slice->dependent_slice_segment_flag == 0)
		slice->SliceAddrRS = slice->slice_segment_address;
	else
		slice->SliceAddrRS = pps->CtbAddrTStoRS[pps->CtbAddrRStoTS[slice->slice_segment_address] - 1];

	if (! slice->dependent_slice_segment_flag) {
		slice->slice_type = rbsp->ue("slice_type");
		if (nal->nal_unit_type >= 7 && nal->nal_unit_type <= 12)
			assert(slice->slice_type == SLICE_TYPE_I);
		if (sps->sps_max_dec_pic_buffering[nal->TemporalId] == 0)
			assert(slice->slice_type == SLICE_TYPE_I);

		slice->pic_output_flag = 1;
		if (pps->output_flag_present_flag)
			slice->pic_output_flag = rbsp->u(1, "pic_output_flag");

		if (pps->separate_colour_plane_flag == 1) {
			slice->colour_plane_id = rbsp->u(2);
			assert(slice->colour_plane_id < COLOUR_PLANE_NUM);
		}

		slice->pic_order_cnt_lsb = 0;
		if (! nal->IdrPicFlag) {
			pic_order_count(slice, rbsp);

			slice->short_term_ref_pic_set_sps_flag = rbsp->u(1, "short_term_ref_pic_set_sps_flag");
			if (! slice->short_term_ref_pic_set_sps_flag)
				short_term_ref_pic_set(&pps->short_term_rps, pps->num_short_term_ref_pic_sets, rbsp);
			else {
				// ceil(log2(pps->num_short_term_ref_pic_sets))
				slice->short_term_ref_pic_set_idx = rbsp->ue("short_term_ref_pic_set_idx");
				assert(slice->short_term_ref_pic_set_idx < pps->num_short_term_ref_pic_sets);
			}

			if (slice->short_term_ref_pic_set_sps_flag)
				slice->StRpsIdx = slice->short_term_ref_pic_set_idx;
			else
				slice->StRpsIdx = pps->num_short_term_ref_pic_sets;

			long_term_pics(slice, rbsp);

			if (slice->nal_unit_type == NAL_UNIT_SLICE_BLA ||
				slice->nal_unit_type == NAL_UNIT_SLICE_BLANT ||
				slice->nal_unit_type == NAL_UNIT_SLICE_BLA_N_LP) {
				rps->numNegativePics = 0;
				rps->numPositivePics = 0;
				rps->numLongTermPics = 0;
			}
		}

		slice->slice_sao_luma_flag   = 0;
		slice->slice_sao_chroma_flag = 0;
		if (sps->sample_adaptive_offset_enabled_flag) {
			slice->slice_sao_luma_flag   = rbsp->u(1, "slice_sao_luma_flag");
			slice->slice_sao_chroma_flag = rbsp->u(1, "slice_sao_chroma_flag");
		}

		slice->slice_temporal_mvp_enable_flag = 0;
		slice->num_ref_idx_l0_active_minus1 = pps->num_ref_idx_l0_default_active_minus1;
		slice->num_ref_idx_l1_active_minus1 = pps->num_ref_idx_l1_default_active_minus1;
		if (slice->slice_type != SLICE_TYPE_I) {
			if (sps->sps_temporal_mvp_enable_flag)
				slice->slice_temporal_mvp_enable_flag = rbsp->u(1, "enable_temporal_mvp_flag");
			slice->num_ref_idx_active_override_flag = rbsp->u(1, "num_ref_idx_active_override_flag");
			if (slice->num_ref_idx_active_override_flag) {
				slice->num_ref_idx_l0_active_minus1 = rbsp->ue("num_ref_idx_l0_active_minus1");
				if (slice->slice_type == SLICE_TYPE_B)
					slice->num_ref_idx_l1_active_minus1 = rbsp->ue("num_ref_idx_l1_active_minus1");
				assert(slice->num_ref_idx_l0_active_minus1 < MAX_REF_IDX);
				assert(slice->num_ref_idx_l1_active_minus1 < MAX_REF_IDX);
			}
			if (sps->lists_modification_present_flag)
				ref_pic_list_modification(slice, rbsp);
		}
		slice->mvd_l1_zero_flag = 0;
		if (slice->slice_type == SLICE_TYPE_B)
			slice->mvd_l1_zero_flag = rbsp->u(1, "mvd_l1_zero_flag");
		slice->cabac_init_flag = 0;
		if (pps->cabac_init_present_flag && slice->slice_type != SLICE_TYPE_I)
			slice->cabac_init_flag = rbsp->u(1, "cabac_init_flag");

		slice->collocated_from_l0_flag = 1;
		slice->collocated_ref_idx      = 0;
		if (slice->slice_temporal_mvp_enable_flag) {
			if (slice->slice_type == SLICE_TYPE_B)
				slice->collocated_from_l0_flag = rbsp->u(1, "collocated_from_l0_flag");
			if (slice->slice_type != SLICE_TYPE_I &&
				((  slice->collocated_from_l0_flag && slice->num_ref_idx_l0_active_minus1 > 0) ||
				 (! slice->collocated_from_l0_flag && slice->num_ref_idx_l1_active_minus1 > 0)))
				slice->collocated_ref_idx = rbsp->ue("collocated_ref_idx");
		}

		if ((pps->weighted_pred_flag   && slice->slice_type == SLICE_TYPE_P) ||
			(pps->weighted_bipred_flag && slice->slice_type == SLICE_TYPE_B))
			pred_weight_table(slice, rbsp);

#if 1
		if (slice->slice_type == SLICE_TYPE_I) {
			slice->five_minus_max_num_merge_cand = rbsp->ue("five_minus_max_num_merge_cand");

			slice->MaxNumMergeCand = 5 - slice->five_minus_max_num_merge_cand;
		}
#else
		slice->five_minus_max_num_merge_cand = rbsp->ue("five_minus_max_num_merge_cand");

		slice->MaxNumMergeCand = 5 - slice->five_minus_max_num_merge_cand;
#endif

		slice->slice_qp_delta = rbsp->se("slice_qp_delta");

		slice->SliceQpY = 26 + pps->pic_init_qp_minus26 + slice->slice_qp_delta;
		assert(slice->SliceQpY >= -sps->QpBdOffsetY && slice->SliceQpY < 52);

		slice->slice_cb_qp_offset = 0;
		slice->slice_cr_qp_offset = 0;
		if (pps->pic_slice_chroma_qp_offsets_present_flag) {
			slice->slice_cb_qp_offset = rbsp->se("slice_qp_delta_cb");
			slice->slice_cr_qp_offset = rbsp->se("slice_qp_delta_cr");
			assert(slice->slice_cb_qp_offset >= -12 && slice->slice_cb_qp_offset <= 12);
			assert(slice->slice_cb_qp_offset >= -12 && slice->slice_cb_qp_offset <= 12);
			assert(pps->pic_cb_qp_offset + slice->slice_cb_qp_offset >= -12);
			assert(pps->pic_cb_qp_offset + slice->slice_cb_qp_offset <=  12);
			assert(pps->pic_cr_qp_offset + slice->slice_cr_qp_offset >= -12);
			assert(pps->pic_cr_qp_offset + slice->slice_cr_qp_offset <=  12);
		}

		slice->deblocking_filter_override_flag              = 0;
		slice->slice_disable_deblocking_filter_flag         = pps->pic_disable_deblocking_filter_flag;
		slice->beta_offset_div2                             = pps->beta_offset_div2;
		slice->tc_offset_div2                               = pps->tc_offset_div2;
		slice->slice_loop_filter_across_slices_enabled_flag = pps->loop_filter_across_slices_enabled_flag;
		if (pps->deblocking_filter_control_present_flag) {
			if (pps->deblocking_filter_override_enabled_flag)
				slice->deblocking_filter_override_flag = rbsp->u(1, "deblocking_filter_override_flag");
			if (slice->deblocking_filter_override_flag) {
				slice->slice_disable_deblocking_filter_flag = rbsp->u(1, "slice_disable_deblocking_filter_flag");
				if (! slice->slice_disable_deblocking_filter_flag) {
					slice->beta_offset_div2 = rbsp->se("beta_offset_div2");
					slice->tc_offset_div2   = rbsp->se("tc_offset_div2");
					assert(slice->beta_offset_div2 >= -6 && slice->beta_offset_div2 <= 6);
					assert(slice->tc_offset_div2   >= -6 && slice->tc_offset_div2   <= 6);
				}
			}
		}
		if (pps->loop_filter_across_slices_enabled_flag &&
			((sps->sample_adaptive_offset_enabled_flag &&
			  (slice->slice_sao_luma_flag || slice->slice_sao_chroma_flag)) ||
			 ! slice->slice_disable_deblocking_filter_flag))
			slice->slice_loop_filter_across_slices_enabled_flag = rbsp->u(1, "slice_loop_filter_across_slices_enabled_flag");
	}

	slice->num_entry_point_offsets = 0;
	if (pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag) {
		slice->num_entry_point_offsets = rbsp->ue("num_entry_point_offsets");
		if (pps->tiles_enabled_flag == 0 && pps->entropy_coding_sync_enabled_flag == 1)
			assert(slice->num_entry_point_offsets < sps->PicHeightInCtbsY);
		if (pps->tiles_enabled_flag == 1 && pps->entropy_coding_sync_enabled_flag == 0)
			assert(slice->num_entry_point_offsets <
			 	   (sps->num_tile_columns_minus1+1) * (sps->num_tile_rows_minus1+1));
		if (pps->tiles_enabled_flag == 1 && pps->entropy_coding_sync_enabled_flag == 1)
			assert(slice->num_entry_point_offsets <
				   (sps->num_tile_columns_minus1+1) * sps->PicHeightInCtbsY);

		if (slice->num_entry_point_offsets > 0) {
			uint32_t offset_len_minus1 = rbsp->ue("offset_len_minus1");
			assert(offset_len_minus1 < 32);
			for (int i = 0; i < slice->num_entry_point_offsets; i++)
				slice->entry_point_offset[i] = rbsp->u(offset_len_minus1+1, "entry_point_offset");
		}
	}

	if (pps->slice_segment_header_extension_present_flag) {
		uint32_t slice_segment_header_extenstion_length = rbsp->ue("slice_header_extenstion_length");
		assert(slice_segment_header_extenstion_length <= 256);
		for (int i = 0; i < slice_segment_header_extenstion_length; i++)
			uint32_t slice_segment_header_extenstion_data_byte = rbsp->u(8, "slice_header_extenstion_data_byte");
	}

	rbsp->byte_alignment();
}

void slice_segment_data(rbsp_t *rbsp) {
	do {
		coding_tree_unit(rbsp);
		uint32_t end_of_slice_segment_flag = rbsp->ae();

		CtbAddrTS++;
		CtbAddrRS = CtbAddrTStoRS[CtbAddrTS];

		if (! end_of_slice_segment_flag &&
			((tiles_enabled_flag && TileId[CtbAddrTS] != TileId[CtbAddrTS-1]) ||
			 (entropy_coding_sync_enabled_flag && CtbAddrTS % PicWidthInCtbs == 0))) {
			uint32_t end_of_sub_stream_one_bit = rbsp->ae();
			assert(end_of_sub_stream_one_bit == 1);

			byte_alignment();
		}
	} while (! end_of_slice_segment_flag);
}

void rbsp_slice_segment_trailing_bits(rbsp_t *rbsp) {
	rbsp_railing_bits(rbsp);

	while (more_rbsp_trailing_data()) {
		uint32_t cabac_zero_word = rbsp->f(16);
		assert(cabac_zero_word == 0x0000);
	}
}

void rbsp_railing_bits(rbsp_t *rbsp) {
	uint32_t rbsp_stop_one_bit = rbsp->f(1);
	assert(rbsp_stop_one_bit == 1);

	while (! byte_aligned(rbsp)) {
		uint32_t rbsp_alignment_zero_bit = rbsp->f(1);
		assert(rbsp_alignment_zero_bit == 0);
	}
}

void byte_alignment(rbsp_t *rbsp) {
	uint32_t bit_equal_to_one = rbsp->f(1);
	assert(bit_equal_to_one == 1);

	while (! byte_aligned(rbsp)) {
		uint32_t bit_equal_to_zero = rbsp->f(1);
		assert(bit_equal_to_zero == 0);
	}
}

void pic_order_count(slice_t *slice, rbsp_t *rbsp) {
	slice->pic_order_cnt_lsb = rbsp->u(sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "pic_order_cnt_lsb");
	assert(slice->pic_order_cnt_lsb < sps->MaxPicOrderCntLsb);

	int iMaxPOClsb  = 1 << (sps->log2_max_pic_order_cnt_lsb_minus4 + 4);
	int iPrevPOC    = slice->prevPOC;
	int iPrevPOClsb = iPrevPOC % iMaxPOClsb;
	int iPrevPOCmsb = iPrevPOC - iPrevPOClsb;

	int iPOClsb = slice->pic_order_cnt_lsb;
	int iPOCmsb = iPrevPOCmsb;
	if (iPOClsb < iPrevPOClsb && iPrevPOClsb - iPOClsb >= iMaxPOClsb/2)
		iPOCmsb += iMaxPOClsb;
	if (iPOClsb > iPrevPOClsb && iPOClsb - iPrevPOClsb > iMaxPOClsb/2)
		iPOCmsb -= iMaxPOClsb;
	if (slice->nal_unit_type == NAL_UNIT_SLICE_BLA ||
		slice->nal_unit_type == NAL_UNIT_SLICE_BLANT ||
		slice->nal_unit_type == NAL_UNIT_SLICE_BLA_N_LP)
		iPOCmsb = 0;
	slice->pic_order_count = iPOCmsb + iPOClsb;
}

void long_term_pics(slice_t *slice, rbsp_t *rbsp) {
	slice->num_long_term_sps  = 0;
	slice->num_long_term_pics = 0;

	if (sps->long_term_ref_pics_present_flag) {
		if (sps->num_long_term_ref_pics_sps > 0) {
			slice->num_long_term_sps = rbsp->ue("num_long_term_sps");
			assert(slice->num_long_term_sps <=
			   min(sps->num_long_term_ref_pics_sps,
			   	   sps->sps_max_dec_pic_buffering[sps->sps_max_temporal_layers_minus1]
			   	   - rps[slice->StRpsIdx].NumNegativePics
			   	   - rps[slice->StRpsIdx].NumPositivePics));
		}

		slice->num_long_term_pics = rbsp->ue("num_long_term_pics");
		assert(slice->num_long_term_pics <=
		 	   sps->sps_max_dec_pic_buffering[sps->sps_max_temporal_layers_minus1]
		 	   - rps[slice->StRpsIdx].NumNegativePics
		 	   - rps[slice->StRpsIdx].NumPositivePics
		 	   - slice->num_long_term_sps));

		for (int i = 0; i < slice->num_long_term_sps + slice->num_long_term_pics; i++) {
			if (i < slice->num_long_term_sps) {
				uint32_t lt_idx_sps = rbsp->u(ceil(log2(sps->num_long_term_ref_pics_sps)), "lt_idx_sps[i]");
				assert(lt_idx_sps < sps->num_long_term_ref_pics_sps);

				slice->PocLsbLt       [i] = sps->lt_ref_pic_poc_lsb_sps      [lt_idx_sps];
				slice->UsedByCurrPicLt[i] = sps->used_by_curr_pic_lt_sps_flag[lt_idx_sps];
			}
			else {
				uint32_t poc_lsb_lt               = rbsp->u(sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "poc_lsb_lt");
				uint32_t used_by_curr_pic_lt_flag = rbsp->u(1, "used_by_curr_pic_lt_flag");

				slice->PocLsbLt       [i] = poc_lsb_lt;
				slice->UsedByCurrPicLt[i] = used_by_curr_pic_lt_flag;
			}

			uint32_t delta_poc_msb_present_flag = rbsp->u(1, "delta_poc_msb_present_flag");
			if (delta_poc_msb_present_flag) {
				uint32_t delta_poc_msb_cycle_lt = rbsp->ue("delta_poc_msb_cycle_lt[i]");
				slice->DeltaPocMSBCycleLt[i] = delta_poc_msb_cycle_lt;
				if (! (i == 0 || i == slice->num_long_term_sps ||
				 	   slice->PocLsbLt[i-1] != slice->PocLsbLt[i]))
					slice->DeltaPocMSBCycleLt[i] += slice->DeltaPocMSBCycleLt[i-1];

				assert(slice->DeltaPocMSBCycleLt[i] * sps->MaxPicOrderCntLsb
					   + slice->pic_order_cnt_lsb - slice->PocLsbLt[i] >= 1);
				assert(slice->DeltaPocMSBCycleLt[i] * sps->MaxPicOrderCntLsb
					   + slice->pic_order_cnt_lsb - slice->PocLsbLt[i] < (1 << 24));
			}
		}
	}

	slice->NumPocTotalCurr = 0;
	for (int i = 0; i < rps[slice->StRpsIdx].NumNegativePics; i++) {
		if (rps[slice->StRpsIdx].UsedByCurrPicS0[i])
			slice->NumPocTotalCurr++;
	}
	for (int i = 0; i < rps[slice->StRpsIdx].NumPositivePics; i++) {
		if (rps[slice->StRpsIdx].UsedByCurrPicS1[i])
			slice->NumPocTotalCurr++;
	}
	for (int i = 0; i < slice->num_long_term_sps + slice->num_long_term_pics; i++) {
		if (rps[slice->StRpsIdx].UsedByCurrPicLt[i])
			slice->NumPocTotalCurr++;
	}
}
