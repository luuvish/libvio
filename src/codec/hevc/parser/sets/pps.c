#include "rbsp.h"
#include "pps.h"

void pic_parameter_set_rbsp(pps_t *pps, rbsp_t *rbsp) {
	pps->pic_parameter_set_id = rbsp->ue("pic_parameter_set_id");
	assert(pps->pic_parameter_set_id < MAX_PPS_ID);
	pps->seq_parameter_set_id = rbsp->ue("seq_parameter_set_id");
	assert(pps->seq_parameter_set_id < MAX_SPS_ID);

	pps->sign_data_hiding_flag   = rbsp->u(1, "sign_data_hiding_flag");
	pps->cabac_init_present_flag = rbsp->u(1, "cabac_init_present_flag");

	pps->num_ref_idx_l0_default_active_minus1 = rbsp->ue("num_ref_idx_l0_default_active_minus1");
	pps->num_ref_idx_l1_default_active_minus1 = rbsp->ue("num_ref_idx_l1_default_active_minus1");
	assert(pps->num_ref_idx_l0_default_active_minus1 < MAX_REF_IDX);
	assert(pps->num_ref_idx_l1_default_active_minus1 < MAX_REF_IDX);

	pps->pic_init_qp_minus26 = rbsp->se("pic_init_qp_minus26");
	assert(pps->pic_init_qp_minus26 >= -(26 + sps->QpBdOffsetY))
	assert(pps->pic_init_qp_minus26 < 26);

	pps->constrained_intra_pred_flag = rbsp->u(1, "constrained_intra_pred_flag");
	pps->transform_skip_enabled_flag = rbsp->u(1, "transform_skip_enabled_flag");
	pps->cu_qp_delta_enabled_flag    = rbsp->u(1, "cu_qp_delta_enabled_flag");
	if (pps->cu_qp_delta_enabled_flag)
		pps->diff_cu_qp_delta_depth  = rbsp->ue("diff_cu_qp_delta_depth");
	assert(pps->diff_cu_qp_delta_depth <= sps->log2_diff_max_min_luma_coding_block_size);

	pps->Log2MinCuQpDeltaSize = sps->Log2CtbSizeY - pps->diff_cu_qp_delta_depth;

	pps->pic_cb_qp_offset = rbsp->se("cb_qp_offset");
	pps->pic_cr_qp_offset = rbsp->se("cr_qp_offset");
	assert(pps->pic_cb_qp_offset >= -12 && pps->pic_cb_qp_offset <= 12);
	assert(pps->pic_cr_qp_offset >= -12 && pps->pic_cr_qp_offset <= 12);

	pps->pic_slice_chroma_qp_offsets_present_flag = rbsp->u(1, "slicelevel_chroma_qp_flag");

	pps->weighted_pred_flag                    = rbsp->u(1, "weighted_pred_flag");
	pps->weighted_bipred_flag                  = rbsp->u(1, "weighted_bipred_flag");
	pps->output_flag_present_flag              = rbsp->u(1, "output_flag_present_flag");
	pps->transquant_bypass_enable_flag         = rbsp->u(1, "transquant_bypass_enable_flag");
	pps->dependent_slice_segments_enabled_flag = rbsp->u(1, "dependent_slices_enabled_flag");

	pps->tiles_enabled_flag                    = rbsp->u(1, "tiles_enabled_flag");
	pps->entropy_coding_sync_enabled_flag      = rbsp->u(1, "entropy_coding_sync_enabled_flag");
#if 1
	pps->entropy_slice_enabled_flag            = rbsp->u(1, "entropy_slice_enabled_flag");
#endif
	pps->loop_filter_across_tiles_enabled_flag = 0;
	if (pps->tiles_enabled_flag) {
		pps->num_tile_columns_minus1 = rbsp->ue("num_tile_columns_minus1");
		pps->num_tile_rows_minus1    = rbsp->ue("num_tile_rows_minus1");
		assert(pps->column_width_minus1[i] < sps->picWidthInCtbsY);
		assert(pps->row_height_minus1  [i] < sps->picHeightInCtbsY);

		pps->uniform_spacing_flag = rbsp->u(1, "uniform_spacing_flag");
		if (! pps->uniform_spacing_flag) {
			for (int i = 0; i < pps->num_tile_columns_minus1; i++)
				pps->column_width_minus1[i] = rbsp->ue("column_width_minus1");
			for (int i = 0; i < pps->num_tile_rows_minus1; i++)
				pps->row_height_minus1[i]   = rbsp->ue("row_height_minus1");
		}

		if (pps->uniform_spacing_flag) {
			for (int i = 0; i <= pps->num_tile_columns_minus1; i++)
				pps->colWidth[i] = (i+1) * sps->PicWidthInCtbsY / (pps->num_tile_columns_minus1+1)
								 -  i    * sps->PicWidthInCtbsY / (pps->num_tile_columns_minus1+1);
			for (int j = 0; j <= pps->num_tile_rows_minus1; j++)
				pps->rowHeight[j] = (j+1) * sps->PicHeightInCtbsY / (pps->num_tile_rows_minus1+1)
								  -  j    * sps->PicHeightInCtbsY / (pps->num_tile_rows_minus1+1);
		}
		else {
			pps->colWidth[pps->num_tile_columns_minus1] = sps->PicWidthInCtbsY;
			for (int i = 0; i < pps->num_tile_columns_minus1; i++) {
				pps->colWidth[i] = pps->column_width_minus1[i] + 1;
				pps->colWidth[pps->num_tile_columns_minus1] -= pps->colWidth[i];
			}
			pps->rowHeight[pps->num_tile_rows_minus1] = sps->PicHeightInCtbsY;
			for (int j = 0; j < pps->num_tile_rows_minus1; j++) {
				pps->rowHeight[j] = pps->row_height_minus1[j] + 1;
				pps->rowHeight[pps->num_tile_rows_minus1] -= pps->rowHeight[j];
			}
		}

		for (pps->colBd[0] = 0, i = 0; i <= pps->num_tile_columns_minus1; i++)
			pps->colBd[i+1] = pps->colBd[i] + pps->colWidth[i];
		for (pps->rowBd[0] = 0, j = 0; j <= pps->num_tile_rows_minus1; j++)
			pps->rowBd[j+1] = pps->rowBd[j] + pps->rowHeight[j];

		for (int ctbAddrRS = 0; ctbAddrRS < sps->PicSizeInCtbsY; ctbAddrRS++) {
			int tbX = ctbAddrRS % sps->PicWidthInCtbsY;
			int tbY = ctbAddrRS / sps->PicWidthInCtbsY;
			int tileX = 0;
			int tileY = 0;
			for (int i = 0; i <= pps->num_tile_columns_minus1; i++) {
				if (tbX >= pps->colBd[i])
					tileX = i;
			}
			for (int j = 0; j <= pps->num_tile_rows_minus1; i++) {
				if (tbY >= pps->rowBd[j])
					tileY = j;
			}
			pps->CtbAddrRStoTS[ctbAddrRS] = 0;
			for (int i = 0; i < tileX; i++)
				pps->CtbAddrRStoTS[ctbAddrRS] += pps->rowHeight[tileY] * pps->colWidth[i];
			for (int j = 0; j < tileY; j++)
				pps->CtbAddrRStoTS[ctbAddrRS] += sps->PicWidthInCtbsY * pps->rowHeight[j];
			pps->CtbAddrRStoTS[ctbAddrRS] += (tbY - pps->rowBd[tileY]) * pps->colWidth[tileX]
			 							   + (tbX - pps->colBd[tileX]);
		}

		for (int ctbAddrRS = 0; ctbAddrRS < sps->PicSizeInCtbsY; ctbAddrRS++)
			pps->CtbAddrTStoRS[pps->CtbAddrRStoTS[ctbAddrRS]] = ctbAddrRS;

		for (int j = 0, iIdx = 0; j <= pps->num_tile_rows_minus1; j++) {
			for (int i = 0; i <= pps->num_tile_columns_minus1; i++, tIdx++)
				for (int y = pps->rowBd[j]; y < pps->rowBd[j+1]; y++)
					for (int x = pps->colBd[i]; x < pps->colBd[i+1]; x++)
						pps->TileId[pps->CtbAddrRStoTS[y * sps->PicWidthInCtbsY + x]] = tIdx;
		}

		for (int i = 0; i <= pps->num_tile_columns_minus1; i++)
			pps->ColumnWidthInLumaSamples[i] = pps->colWidth[i] << sps->Log2CtbSizeY;
		for (int j = 0; j <= pps->num_tile_rows_minus1; j++)
			pps->RowHeightInLumaSamples[i] = pps->rowHeight[j] << sps->Log2CtbSizeY;

#if 1
		if (pps->num_tile_columns_minus1 != 0 || pps->num_tile_rows_minus1 != 0)
#endif
		pps->loop_filter_across_tiles_enabled_flag = rbsp->u(1, "loop_filter_across_tiles_enabled_flag");
	}

	for (int y = 0; y < sps->PicHeightInCtbsY << (sps->Log2CtbSizeY - sps->Log2MinTrafoSize); y++) {
		for (int x = 0; x < sps->PicWidthInCtbsY << (sps->Log2CtbSizeY - sps->Log2MinTrafoSize); x++) {
			int tbX = (x << sps->Log2MinTrafoSize) >> sps->Log2CtbSizeY;
			int tbY = (y << sps->Log2MinTrafoSize) >> sps->Log2CtbSizeY;
			int ctbAddrRS = sps->PicWidthInCtbsY * tbY + tbX;
			pps->MinTbAddrZS[x][y] = pps->CtbAddrRStoTS[ctbAddrRS] << ((sps->Log2CtbSizeY - sps->Log2MinTrafoSize) * 2);
			for (int i = 0, p = 0; i < (sps->Log2CtbSizeY - sps->Log2MinTrafoSize); i++) {
				int m = 1 << i;
				p += (m & x ? m * m : 0) + (m & y ? 2 * m * m : 0);
			}
			pps->MinTbAddrZS[x][y] += p;
		}
	}

	pps->loop_filter_across_slices_enabled_flag  = rbsp->u(1, "loop_filter_across_slices_enabled_flag");
	pps->deblocking_filter_control_present_flag  = rbsp->u(1, "deblocking_filter_control_present_flag");
	pps->deblocking_filter_override_enabled_flag = 0;
	pps->pic_disable_deblocking_filter_flag      = 0;
	pps->beta_offset_div2                        = 0;
	pps->tc_offset_div2                          = 0;
	if (pps->deblocking_filter_control_present_flag) {
		pps->deblocking_filter_override_enabled_flag = rbsp->u(1, "deblocking_filter_override_enabled_flag");
		pps->pic_disable_deblocking_filter_flag      = rbsp->u(1, "pic_disable_deblocking_filter_flag");
		if (! pps->pic_disable_deblocking_filter_flag) {
			pps->beta_offset_div2 = rbsp->se("pps_beta_offset_div2");
			pps->tc_offset_div2   = rbsp->se("pps_tc_offset_div2");
			assert(pps->beta_offset_div2 >= -6 && pps->beta_offset_div2 <= 6);
			assert(pps->tc_offset_div2   >= -6 && pps->tc_offset_div2   <= 6);
		}
	}

	pps->pic_scaling_list_data_present_flag = rbsp->u(1, "pps_scaling_list_data_present_flag");
	if (sps->scaling_list_enable_flag == 0)
		assert(pps->pic_scaling_list_data_present_flag == 0);
	if (pps->pic_scaling_list_data_present_flag)
		scaling_list_data(&pps->pic_scaling_list, rbsp);

	pps->log2_parallel_merge_level_minus2            = rbsp->ue("log2_parallel_merge_level_minus2");
	assert(pps->log2_parallel_merge_level_minus2 <= sps->Log2CtbSizeY - 2);
	pps->slice_segment_header_extension_present_flag = rbsp->u(1, "slice_header_extension_present_flag");

	uint32_t pps_extension_flag = rbsp->u(1, "pps_extension_flag");
	if (pps_extension_flag) {
		while (rbsp->more_rbsp_data())
			uint32_t pps_extension_data_flag = rbsp->u(1, "pps_extension_data_flag");
	}

	rbsp->rbsp_trailing_bits();
}
