#include "rbsp.h"
#include "sps.h"

enum {
	SCAN_IDX_DIAG = 0,
	SCAN_IDX_HOR  = 1,
	SCAN_IDX_VER  = 2,
	SCAN_IDX_NUM  = 3
};

enum {
	SCAN_COMP_HOR = 0,
	SCAN_COMP_VER = 1,
	SCAN_COMP_NUM = 2
};

void seq_parameter_set_rbsp(sps_t *sps, rbsp_t *rbsp) {
	sps->video_parameter_set_id    = rbsp->u(4, "video_parameter_set_id");
	assert(sps->video_parameter_set_id < MAX_VPS_ID);
	sps->sps_max_sub_layers_minus1 = rbsp->u(3, "sps_max_sub_layers_minus1");
	assert(sps->sps_max_sub_layers_minus1 < MAX_SUB_LAYERS);
	uint32_t sps_reverved_zero_bit = rbsp->u(1, "sps_reverved_zero_bit");
	assert(sps_reverved_zero_bit == 0);

	profile_tier_level(&sps->ptl, 1, sps->sps_max_sub_layers_minus1, rbsp);

	sps->seq_parameter_set_id       = rbsp->ue("seq_parameter_set_id");
	assert(sps->seq_parameter_set_id < MAX_SPS_ID);
	sps->chroma_format_idc          = rbsp->ue("chroma_format_idc");
	assert(sps->chroma_format_idc <= CHROMA_FORMAT_IDC_444);
	sps->separate_colour_plane_flag = 0;
	if (sps->chroma_format_idc == CHROMA_FORMAT_IDC_444)
		sps->separate_colour_plane_flag = rbsp->u(1, "separate_colour_plane_flag");

	sps->ChromaArrayType = CHROMA_FORMAT_IDC_MONO;
	if (sps->separate_colour_plane_flag == 0)
		sps->ChromaArrayType = sps->chroma_format_idc;
	int SubWidthC [] = {1, 2, 2, 1, 1};
	int SubHeightC[] = {1, 2, 1, 1, 1};
	sps->SubWidthC  = SubWidthC [sps->chroma_format_idc];
	sps->SubHeightC = SubHeightC[sps->chroma_format_idc];

	sps->pic_width_in_luma_samples  = rbsp->ue("pic_width_in_luma_samples");
	sps->pic_height_in_luma_samples = rbsp->ue("pic_height_in_luma_samples");
	sps->pic_cropping_flag          = rbsp->u(1, "pic_cropping_flag");
	sps->pic_crop_left_offset       = 0;
	sps->pic_crop_right_offset      = 0;
	sps->pic_crop_top_offset        = 0;
	sps->pic_crop_bottom_offset     = 0;
	if (sps->pic_cropping_flag) {
		sps->pic_crop_left_offset   = rbsp->ue("pic_crop_left_offset");
		sps->pic_crop_right_offset  = rbsp->ue("pic_crop_right_offset");
		sps->pic_crop_top_offset    = rbsp->ue("pic_crop_top_offset");
		sps->pic_crop_bottom_offset = rbsp->ue("pic_crop_bottom_offset");

		int CropUnitX = 1;
		int CropUnitY = 1;
		if (sps->ChromaArrayType != CHROMA_FORMAT_IDC_MONO) {
			CropUnitX = sps->SubWidthC;
			CropUnitY = sps->SubHeightC;
		}
		assert(sps->pic_crop_left_offset + sps->pic_crop_right_offset  < sps->pic_width_in_luma_samples  / CropUnitX);
		assert(sps->pic_crop_top_offset  + sps->pic_crop_bottom_offset < sps->pic_height_in_luma_samples / CropUnitY);
	}

	sps->bit_depth_luma_minus8   = rbsp->ue("bit_depth_luma_minus8");
	sps->bit_depth_chroma_minus8 = rbsp->ue("bit_depth_chroma_minus8");
	assert(sps->bit_depth_luma_minus8   <= 6);
	assert(sps->bit_depth_chroma_minus8 <= 6);

	sps->BitDepthY   = 8 + sps->bit_depth_luma_minus8;
	sps->QpBdOffsetY = 6 * sps->bit_depth_luma_minus8;
	sps->BitDepthC   = 8 + sps->bit_depth_chroma_minus8;
	sps->QpBdOffsetC = 6 * sps->bit_depth_chroma_minus8;

#if 1
	sps->pcm_enabled_flag                    = rbsp->u(1, "pcm_enabled_flag");
	sps->pcm_loop_filter_disable_flag        = 0;
	if (sps->pcm_enabled_flag) {
		sps->pcm_sample_bit_depth_luma_minus1   = rbsp->u(4, "pcm_bit_depth_luma_minus1");
		sps->pcm_sample_bit_depth_chroma_minus1 = rbsp->u(4, "pcm_bit_depth_chroma_minus1");
	}
#endif

	sps->log2_max_pic_order_cnt_lsb_minus4 = rbsp->ue("log2_max_pic_order_cnt_lsb_minus4");
	assert(sps->log2_max_pic_order_cnt_lsb_minus4 <= 12);

	sps->MaxPicOrderCntLsb = 1 << (sps->log2_max_pic_order_cnt_lsb_minus4 + 4);

	for (int i = 0; i <= sps->sps_max_sub_layers_minus1; i++) {
		sps->sps_max_dec_pic_buffering[i] = rbsp->ue("max_dec_pic_buffering");
		sps->sps_num_reorder_pics     [i] = rbsp->ue("num_reorder_pics");
		sps->sps_max_latency_increase [i] = rbsp->ue("max_latency_increase");
		assert(sps->sps_max_dec_pic_buffering[i] <= MaxDpbSize);
		assert(sps->sps_num_reorder_pics     [i] <= sps->sps_max_dec_pic_buffering[i]);
		assert(sps->sps_max_latency_increase [i] <= ((1<<32)-2));
		if (i > 0) {
			assert(sps->sps_max_dec_pic_buffering[i-1] <= sps->sps_max_dec_pic_buffering[i]);
			assert(sps->sps_num_reorder_pics     [i-1] <= sps->sps_num_reorder_pics     [i]);
		}
	}

	sps->restricted_ref_pic_lists_flag       = rbsp->u(1, "restricted_ref_pic_lists_flag");
	sps->lists_modification_present_flag     = 0;
	if (sps->restricted_ref_pic_lists_flag)
		sps->lists_modification_present_flag = rbsp->u(1, "lists_modification_present_flag");

	sps->log2_min_luma_coding_block_size_minus3   = rbsp->ue("log2_min_coding_block_size_minus3");
	sps->log2_diff_max_min_luma_coding_block_size = rbsp->ue("log2_diff_max_min_coding_block_size");

	sps->Log2MinCbSizeY     = sps->log2_min_luma_coding_block_size_minus3 + 3;
	sps->Log2CtbSizeY       = sps->Log2MinCbSizeY + sps->log2_diff_max_min_luma_coding_block_size;
	sps->MinCbSizeY         = 1 << sps->Log2MinCbSizeY;
	sps->CtbSizeY           = 1 << sps->Log2CtbSizeY;
	sps->PicWidthInMinCbsY  = sps->pic_width_in_luma_samples  / sps->MinCbSizeY;
	sps->PicWidthInCtbsY    = (sps->pic_width_in_luma_samples  + sps->CtbSizeY) / sps->CtbSizeY;
	sps->PicHeightInMinCbsY = sps->pic_height_in_luma_samples / sps->MinCbSizeY;
	sps->PicHeightInCtbsY   = (sps->pic_height_in_luma_samples + sps->CtbSizeY) / sps->CtbSizeY;
	sps->PicSizeInMinCbsY   = sps->PicWidthInMinCbsY * sps->PicHeightInMinCbsY;
	sps->PicSizeInCtbsY     = sps->PicWidthInCtbsY   * sps->PicHeightInCtbsY;
	sps->PicSizeInSamplesY  = sps->pic_width_in_luma_samples * sps->pic_height_in_luma_samples;
	sps->CtbWidthC          = 0;
	sps->CtbHeightC         = 0;
	if (sps->chroma_format_idc != CHROMA_FORMAT_IDC_MONO && sps->separate_colour_plane_flag == 0) {
		sps->CtbWidthC      = sps->CtbSizeY / sps->SubWidthC;
		sps->CtbHeightC     = sps->CtbSizeY / sps->SubHeightC;
	}
	assert(sps->pic_width_in_luma_samples  % sps->MinCbSizeY == 0);
	assert(sps->pic_height_in_luma_samples % sps->MinCbSizeY == 0);

	sps->log2_min_transform_block_size_minus2   = rbsp->ue("log2_min_transform_block_size_minus2");
	sps->log2_diff_max_min_transform_block_size = rbsp->ue("log2_diff_max_min_transform_block_size");

	sps->Log2MinTrafoSize = sps->log2_min_transform_block_size_minus2 + 2;
	sps->Log2MaxTrafoSize = sps->Log2MinTrafoSize + sps->log2_diff_max_min_transform_block_size;
	assert(sps->Log2MinTrafoSize < sps->Log2MinCbSizeY);
	assert(sps->Log2MaxTrafoSize <= min(sps->Log2CtbSizeY, 5));
	// ScanOrder[log2BlockSize][scanIdx][sPos][sComp]
	// log2BlockSize ranging from min(2, Log2MinTrafoSize-2) to max(2, Log2MaxTrafoSize-2)
	for (int log2BlockSize = min(2, sps->Log2MinTrafoSize-2); log2BlockSize <= 3; log2BlockSize++)
		scan_order(sps->ScanOrder[log2BlockSize], 1 << log2BlockSize);

#if 1
	if (sps->pcm_enabled_flag) {
		sps->log2_min_pcm_luma_coding_block_size_minus3   = rbsp->ue("log2_min_pcm_coding_block_size_minus3");
		sps->log2_diff_max_min_pcm_luma_coding_block_size = rbsp->ue("log2_diff_max_min_pcm_coding_block_size");
	}
#endif

	sps->max_transform_hierarchy_depth_inter = rbsp->ue("max_transform_hierarchy_depth_inter");
	sps->max_transform_hierarchy_depth_intra = rbsp->ue("max_transform_hierarchy_depth_intra");
	assert(sps->max_transform_hierarchy_depth_inter <= sps->Log2CtbSizeY - sps->Log2MinTrafoSize);
	assert(sps->max_transform_hierarchy_depth_intra <= sps->Log2CtbSizeY - sps->Log2MinTrafoSize);

	sps->scaling_list_enable_flag               = rbsp->u(1, "scaling_list_enable_flag");
	sps->sps_scaling_list_data_present_flag     = 0;
	if (sps->scaling_list_enable_flag) {
		sps->sps_scaling_list_data_present_flag = rbsp->u(1, "sps_scaling_list_data_present_flag");
		if (sps->sps_scaling_list_data_present_flag)
			scaling_list_data(&sps->sps_scaling_list, rbsp);
	}

	sps->amp_enabled_flag                    = rbsp->u(1, "asymmetric_motion_partitions_enabled_flag");
	sps->sample_adaptive_offset_enabled_flag = rbsp->u(1, "sample_adaptive_offset_enabled_flag");
#if 1
	sps->pcm_loop_filter_disable_flag     = 0;
	if (sps->pcm_enabled_flag)
		sps->pcm_loop_filter_disable_flag = rbsp->u(1, "pcm_loop_filter_disable_flag");
#endif
#if 0
	sps->pcm_enabled_flag                    = rbsp->u(1, "pcm_enabled_flag");
	sps->pcm_loop_filter_disable_flag        = 0;
	if (sps->pcm_enabled_flag) {
		sps->pcm_sample_bit_depth_luma_minus1   = rbsp->u(4, "pcm_sample_bit_depth_luma_minus1");
		sps->pcm_sample_bit_depth_chroma_minus1 = rbsp->u(4, "pcm_sample_bit_depth_chroma_minus1");

		sps->PCMBitDepthY = 1 + sps->pcm_sample_bit_depth_luma_minus1;
		sps->PCMBitDepthC = 1 + sps->pcm_sample_bit_depth_chroma_minus1;
		assert(sps->PCMBitDepthY <= sps->BitDepthY);
		assert(sps->PCMBitDepthC <= sps->BitDepthC);

		sps->log2_min_pcm_luma_coding_block_size_minus3   = rbsp->ue("log2_min_pcm_coding_block_size_minus3");
		sps->log2_diff_max_min_pcm_luma_coding_block_size = rbsp->ue("log2_diff_max_min_pcm_coding_block_size");
		sps->pcm_loop_filter_disable_flag                 = rbsp->u(1, "pcm_loop_filter_disable_flag");

		sps->Log2MinIpcmCbSizeY = sps->log2_min_pcm_luma_coding_block_size_minus3 + 3;
		sps->Log2MaxIpcmCbSizeY = sps->Log2MaxIpcmCbSizeY + sps->log2_diff_max_min_pcm_luma_coding_block_size;
		assert(sps->Log2MinIpcmCbSizeY >= sps->Log2MinCbSizeY);
		assert(sps->Log2MinIpcmCbSizeY <= min(sps->Log2CtbSizeY, 5));
		assert(sps->Log2MaxIpcmCbSizeY <= min(sps->Log2CtbSizeY, 5));
	}
#endif
	sps->sps_temporal_id_nesting_flag = rbsp->u(1, "temporal_id_nesting_flag");
	assert(sps->sps_temporal_id_nesting_flag == vps->vps_temporal_id_nesting_flag);

	sps->num_short_term_ref_pic_sets = rbsp->ue("num_short_term_ref_pic_sets");
	assert(sps->num_short_term_ref_pic_sets <= MAX_SHORT_TERM_REF_PIC_SETS);
	for (int i = 0; i < sps->num_short_term_ref_pic_sets; i++)
		short_term_ref_pic_set(&sps->rps[i], i, rbsp);
	sps->long_term_ref_pics_present_flag = rbsp->u(1, "long_term_ref_pics_present_flag");
	if (sps->long_term_ref_pics_present_flag) {
		sps->num_long_term_ref_pics_sps = rbsp->ue("num_long_term_ref_pic_sps");
		for (int i = 0; i < sps->num_long_term_ref_pics_sps; i++) {
			sps->lt_ref_pic_poc_lsb_sps[i]       = rbsp->u(sps->log2_max_pic_order_cnt_lsb_minus4 + 4, "lt_ref_pic_poc_lsb_sps");
			sps->used_by_curr_pic_lt_sps_flag[i] = rbsp->u(1, "used_by_curr_pic_lt_sps_flag[i]");
		}
	}
	sps->sps_temporal_mvp_enable_flag       = rbsp->u(1, "sps_temporal_mvp_enable_flag");
#if 0
	sps->strong_intra_smoothing_enable_flag = rbsp->u(1, "strong_intra_smoothing_enable_flag");
#endif

	sps->vui_parameters_present_flag = rbsp->u(1, "vui_parameters_present_flag");
	if (sps->vui_parameters_present_flag)
		vui_parameters(&sps->vui, sps, rbsp);

	uint32_t sps_extension_flag = rbsp->u(1, "sps_extension_flag");
	if (sps_extension_flag) {
		while (rbsp->more_rbsp_data())
			uint32_t sps_extension_data_flag = rbsp->u(1, "sps_extension_data_flag");
	}

	rbsp->rbsp_trailing_bits();
}

void scan_order(int ScanOrder[SCAN_IDX_NUM][64][SCAN_COMP_NUM], int blkSize) {
	bool stopLoop;
	int x, y, i;

	i = 0;
	x = 0;
	y = 0;
	stopLoop = false;
	while (! stopLoop) {
		while (y >= 0) {
			if (x < blkSize && y < blkSize) {
				ScanOrder[SCAN_IDX_DIAG][i][SCAN_COMP_HOR] = x;
				ScanOrder[SCAN_IDX_DIAG][i][SCAN_COMP_VER] = y;
				i++;
			}
			y--;
			x++;
		}
		y = x;
		x = 0;
		if (i >= blkSize * blkSize)
			stopLoop = true;
	}

	i = 0;
	y = 0;
	while (y < blkSize) {
		x = 0;
		while (x < blkSize) {
			ScanOrder[SCAN_IDX_HOR][i][SCAN_COMP_HOR] = x;
			ScanOrder[SCAN_IDX_HOR][i][SCAN_COMP_VER] = y;
			x++;
			i++;
		}
		y++;
	}

	i = 0;
	x = 0;
	while (x < blkSize) {
		y = 0;
		while (y < blkSize) {
			ScanOrder[SCAN_IDX_VER][i][SCAN_COMP_HOR] = x;
			ScanOrder[SCAN_IDX_VER][i][SCAN_COMP_VER] = y;
			y++;
			i++;
		}
		x++;
	}
}
