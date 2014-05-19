
enum {
	MODE_INTER = 0,
	MODE_INTRA = 1
};

enum {
	PART_2Nx2N = 0,
	PART_2NxN  = 1,
	PART_Nx2N  = 2,
	PART_NxN   = 3,
	PART_2NxnU = 4,
	PART_2NxnD = 5,
	PART_nLx2N = 6,
	PART_nRx2N = 7
};

void coding_tree_unit(xCtb, yCtb) {
	xCtb = (CtbAddrRS % sps->PicWidthInCtbsY ) << sps->Log2CtbSizeY;
	yCtb = (CtbAddrRS % sps->PicHeightInCtbsY) << sps->Log2CtbSizeY;
	CtbAddrInSliceSeg = CtbAddrRS - slice->slice_segment_address;

	if (slice->slice_sao_luma_flag || slice->slice_sao_chroma_flag)
		sao(xCtb >> sps->Log2CtbSizeY, yCtb >> sps->Log2CtbSizeY, 0);
	coding_quadtree(xCtb, yCtb, sps->Log2CtbSizeY, 0);
}

void coding_quadtree(x0, y0, log2CbSize, ctDepth, rbsp_t *rbsp) {
	if (x0 + (1 << log2CbSize) <= sps->pic_width_in_luma_samples &&
		y0 + (1 << log2CbSize) <= sps->pic_height_in_luma_samples &&
		log2CbSize > sps->Log2MinCbSizeY)
		split_cu_flag[x0][y0] = rbsp->ae();

	if (pps->cu_qp_delta_enabled_flag && log2CbSize >= pps->Log2MinCuQpDeltaSize) {
		IsCuQpDeltaCoded = 0;
		CuQpDelta = 0;
	}

	if (split_cu_flag[x0][y0]) {
		x1 = x0 + ((1 << log2CbSize) >> 1);
		y1 = y0 + ((1 << log2CbSize) >> 1);
		coding_quadtree(x0, y0, log2CbSize-1, ctDepth+1, rbsp);
		if (x1 < sps->pic_width_in_luma_samples)
			coding_quadtree(x1, y0, log2CbSize-1, ctDepth+1, rbsp);
		if (y1 < sps->pic_height_in_luma_samples)
			coding_quadtree(x0, y1, log2CbSize-1, ctDepth+1, rbsp);
		if (x1 < sps->pic_width_in_luma_samples &&
		 	y1 < sps->pic_height_in_luma_samples)
			coding_quadtree(x1, y1, log2CbSize-1, ctDepth+1, rbsp);
	}
	else
		coding_unit(x0, y0, log2CbSize, rbsp);
}

void coding_unit(x0, y0, log2CbSize, rbsp_t *rbsp) {
	cu->cu_transquant_bypass_flag = 0;
	if (pps->transquant_bypass_enable_flag)
		cu->cu_transquant_bypass_flag = rbsp->ae();

	cu->skip_flag[x0][y0] = 0;
	if (slice_type != SLICE_TYPE_I)
		cu->skip_flag[x0][y0] = rbsp->ae();

	int nCbS = 1 << log2CbSize;

	if (slice->slice_type == SLICE_TYPE_I)
		cu->PredMode[x0][y0] = MODE_INTER;
	else if (cu->skip_flag[x0][y0])
		cu->PredMode[x0][y0] = MODE_SKIP;
	else {
		pred_mode_flag = rbsp->ae(CTX_PRED_MODE, "pred_mode_flag");
		cu->PredMode[x0][y0] = pred_mode_flag == 0 ? MODE_INTER : MODE_INTRA;
	}
	for (int y = y0; y < y0 + nCbS; y++) {
		for (int x = x0; x < x0 + nCbS; x++)
			cu->PredMode[x][y] = cu->PredMode[x0][y0];
	}

	cu->PartMode       = PART_2Nx2N;
	cu->IntraSplitFlag = 0;
	if (! cu->skip_flag[x0][y0] &&
		(cu->PredMode[x0][y0] != MODE_INTRA || log2CbSize == sps->Log2MinCbSizeY)) {
		part_mode = rbsp->ae();
		if (cu->PredMode[x0][y0] == MODE_INTRA)
			assert(part_mode == 0 || part_mode == 1);
		if (cu->PredMode[x0][y0] == MODE_INTER) {
			if (log2CbSize > sps->Log2MinCbSizeY) {
				if (amp_enabled_flag == 1)
					assert((part_mode >= 0 && part_mode <= 2) ||
						   (part_mode >= 4 && part_mode <= 7));
				else
					assert(part_mode >= 0 && part_mode <= 2);
			}
			if (log2CbSize == 3)
				assert(part_mode >= 0 && part_mode <= 2);
			if (log2CbSize > 3)
				assert(part_mode >= 0 && part_mode <= 3);
		}

		cu->PartMode       = part_mode;
		cu->IntraSplitFlag = 0;
		if (cu->PredMode[x0][y0] == MODE_INTRA) {
			cu->PartMode       = part_mode == 0 ? PART_2Nx2N : PART_NxN;
			cu->IntraSplitFlag = part_mode == 0 ? 0 : 1;
		}
	}

	if (cu->PredMode[x0][y0] == MODE_INTRA) {
		cu->pcm_flag[x0][y0] = 0;
		if (cu->PartMode == PART_2Nx2N && sps->pcm_enabled_flag &&
			log2CbSize >= sps->Log2MinIpcmCbSizeY &&
			log2CbSize <= sps->Log2MaxIpcmCbSizeY)
			cu->pcm_flag[x0][y0] = rbsp->ae(CTX_PCM_FLAG, "pcm_flag");
		for (int y = y0; y < y0 + nCbS; y++) {
			for (int x = x0; x < x0 + nCbS; x++)
				cu->pcm_flag[x][y] = cu->pcm_flag[x0][y0];
		}

		if (cu->pcm_flag[x0][y0]) {
			while (! byte_aligned())
				pcm_alignment_zero_bit = rbsp->ae();
			pcm_sample(x0, y0, log2CbSize, rbsp);
		}
	}

	if (cu->PredMode[x0][y0] == MODE_INTRA && ! cu->pcm_flag[x0][y0]) {
		pbOffset = (cu->PartMode == PART_NxN) ? nCbS/2 : nCbS;
		for (int j = 0; j < nCbS; j += pbOffset) {
			for (int i = 0; i < nCbS; i += pbOffset)
				cu->prev_intra_luma_pred_flag[x0 + i][y0 + j] = rbsp->ae(CTX_INTRA_PRED_Y);
		}
		for (int j = 0; j < nCbS; j += pbOffset) {
			for (int i = 0; i < nCbS; i += pbOffset) {
				if (cu->prev_intra_luma_pred_flag[x0 + i][y0 + j])
					cu->mpm_idx[x0 + i][y0 + j] = rbsp->ae();
				else
					cu->rem_intra_luma_pred_mode[x0 + i][y0 + j] = rbsp->ae();
			}
		}
		cu->intra_chroma_pred_mode[x0][y0] = rbsp->ae(CTX_INTRA_PRED_C);
	}

	if (cu->PredMode[x0][y0] != MODE_INTRA) {
		switch (cu->PartMode) {
		case PART_2Nx2N:
			prediction_unit(x0, y0, nCbS, nCbS);
			break;
		case PART_2NxN:
			prediction_unit(x0, y0, nCbS, nCbS/2);
			prediction_unit(x0, y0 + nCbS/2, nCbS, nCbS/2);
			break;
		case PART_Nx2N:
			prediction_unit(x0, y0, nCbS/2, nCbS);
			prediction_unit(x0 + nCbS/2, y0, nCbS/2, nCbS);
			break;
		case PART_2NxnU:
			prediction_unit(x0, y0, nCbS, nCbS/4);
			prediction_unit(x0, y0 + nCbS/4, nCbS, nCbS*3/4);
			break;
		case PART_2NxnD:
			prediction_unit(x0, y0, nCbS, nCbS*3/4);
			prediction_unit(x0, y0 + nCbS*3/4, nCbS, nCbS/4);
			break;
		case PART_nLx2N:
			prediction_unit(x0, y0, nCbS/4, nCbS);
			prediction_unit(x0 + nCbS/4, y0, nCbS*3/4, nCbS);
			break;
		case PART_nRx2N:
			prediction_unit(x0, y0, nCbS*3/4, nCbS);
			prediction_unit(x0 + nCbS*3/4, y0, nCbS/4, nCbS);
			break;
		case PART_NxN:
			prediction_unit(x0, y0, nCbS/2, nCbS/2);
			prediction_unit(x0 + nCbS/2, y0, nCbS/2, nCbS/2);
			prediction_unit(x0, y0 + nCbS/2, nCbS/2, nCbS/2);
			prediction_unit(x0 + nCbS/2, y0 + nCbS/2, nCbS/2, nCbS/2);
			break;
		}
	}

	if (! cu->skip_flag[x0][y0] && ! cu->pcm_flag[x0][y0]) {
		if (cu->PredMode[x0][y0] != PRED_INTRA && ! (cu->PartMode == PART_2Nx2N && merge_flag[x0][y0]))
			no_residual_syntax_flag = rbsp->ae();
		if (! no_residual_syntax_flag) {
			cu->MaxTrafoDepth = cu->PredMode[x0][y0] == PRED_INTRA ?
				max_transform_hierarchy_depth_intra + cu->IntraSplitFlag :
				max_transform_hierarchy_depth_inter;
			transform_tree(x0, y0, x0, y0, log2CbSize, 0, 0);
		}
	}
}

void pcm_sample(x0, y0, log2CbSize, rbsp_t *rbsp) {
	for (int i = 0; i < (1 << (log2CbSize << 1)); i++)
		pcm_sample_luma[i] = rbsp->u(sps->PCMBitDepthY);
	for (int i = 0; i < (1 << (log2CbSize << 1)) >> 1; i++)
		pcm_sample_chroma[i] = rbsp->u(sps->PCMBitDepthC);
}
