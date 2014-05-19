
void transform_tree(x0, y0, xBase, yBase, log2TrafoSize, trafoDepth, blkIdx) {
	interSplitFlag = 0;
	if (max_transform_hierarchy_depth_inter == 0 &&
		PredMode[x0][y0] == MODE_INTER && PartMode != PART_2Nx2N && trafoDepth == 0)
		interSplitFlag = 1;

	split_transform_flag[x0][y0][trafoDepth] = 0;
	if (log2TrafoSize <= Log2MaxTrafoSize &&
	 	log2TrafoSize > Log2MinTrafoSize &&
		trafoDepth < MaxTrafoDepth && !(IntraSplitFlag && trafoDepth == 0))
		split_transform_flag[x0][y0][trafoDepth] = rbsp->ae();
	else if (log2TrafoSize > Log2MaxTrafoSize ||
			 (IntraSplitFlag == 1 && trafoDepth == 0) ||
			 interSplitFlag == 1)
		split_transform_flag[x0][y0][trafoDepth] = 1;

	cbf_cb[x0][y0][trafoDepth] = 0;
	cbf_cr[x0][y0][trafoDepth] = 0;
	if (trafoDepth > 0 && log2TrafoSize == 2) {
		cbf_cb[x0][y0][trafoDepth] = cbf_cb[xBase][yBase][trafoDepth-1];
		cbf_cr[x0][y0][trafoDepth] = cbf_cr[xBase][yBase][trafoDepth-1];
	}
	if (log2TrafoSize > 2) {
		if (trafoDepth == 0 || cbf_cb[xBase][yBase][trafoDepth-1])
			cbf_cb[x0][y0][trafoDepth] = rbsp->ae();
		if (trafoDepth == 0 || cbf_cr[xBase][yBase][trafoDepth-1])
			cbf_cr[x0][y0][trafoDepth] = rbsp->ae();
	}

	if (split_transform_flag[x0][y0][trafoDepth]) {
		x1 = x0 + ((1 << log2TrafoSize) >> 1);
		y1 = y0 + ((1 << log2TrafoSize) >> 1);
		transform_tree(x0, y0, x0, y0, log2TrafoSize-1, trafoDepth+1, 0);
		transform_tree(x1, y0, x0, y0, log2TrafoSize-1, trafoDepth+1, 1);
		transform_tree(x0, y1, x0, y0, log2TrafoSize-1, trafoDepth+1, 2);
		transform_tree(x1, y1, x0, y0, log2TrafoSize-1, trafoDepth+1, 3);
	}
	else {
		cbf_luma[x0][y0][trafoDepth] = 0;
		if (PredMode[x0][y0] == MODE_INTRA || trafoDepth != 0 ||
			cbf_cb[x0][y0][trafoDepth] || cbf_cr[x0][y0][trafoDepth])
			cbf_luma[x0][y0][trafoDepth] = rbsp->ae();
		transform_unit(x0, y0, xBase, yBase, log2TrafoSize, trafoDepth, blkIdx);
	}
}

void transform_unit(x0, y0, xBase, yBase, log2TrafoSize, trafoDepth, blkIdx) {
	if (cbf_luma[x0][y0][trafoDepth] || cbf_cb[x0][y0][trafoDepth] || cbf_cr[x0][y0][trafoDepth]) {
		cu_qp_delta_sign = 0;
		if (cu_qp_delta_enabled_flag && ! IsCuQpDeltaCoded) {
			cu_qp_delta_abs = rbsp->ae();
			if (cu_qp_delta_abs) {
				cu_qp_delta_sign = rbsp->ae();

				IsCuQpDeltaCoded = 1;
				CuQpDelta = cu_qp_delta_abs * (1 - 2 * cu_qp_delta_sign);
				assert(CuQpDelta >= -(26 + QpBdOffsetY/2) && CuQpDelta < (26 + QpBdOffsetY/2));
			}
		}

		if (cbf_luma[x0][y0][trafoDepth])
			residual_coding(x0, y0, log2TrafoSize, 0);
		if (log2TrafoSize > 2) {
			if (cbf_cb[x0][y0][trafoDepth])
				residual_coding(x0, y0, log2TrafoSize, 1);
			if (cbf_cr[x0][y0][trafoDepth])
				residual_coding(x0, y0, log2TrafoSize, 2);
		}
		else if (blkIdx == 3) {
			if (cbf_cb[xBase][yBase][trafoDepth])
				residual_coding(xBase, yBase, log2TrafoSize, 1);
			if (cbf_cr[xBase][yBase][trafoDepth])
				residual_coding(xBase, yBase, log2TrafoSize, 2);
		}
	}
}

void residual_coding(x0, y0, log2TrafoSize, cIdx) {
	scanIdx = 0;
	if (PredMode[x0][y0] == MODE_INTRA) {
		assert((log2TrafoSize == 2) || (log2TrafoSize == 3 && cIdx == 0))

		if (IntraPredMode >= 6 && IntraPredMode <= 14)
			scanIdx = 2;
		if (IntraPredMode >= 22 && IntraPredMode <= 30)
			scanIdx = 1;
	}

	transform_skip_flag[x0][y0][cIdx] = 0;
	if (transform_skip_enabled && ! cu_transquant_bypass_flag && log2TrafoSize == 2)
		transform_skip_flag[x0][y0][cIdx] = rbsp->ae();

	last_significant_coeff_x_prefix = rbsp->ae();
	last_significant_coeff_y_prefix = rbsp->ae();
	assert(last_significant_coeff_x_prefix < (log2TrafoSize << 1));
	assert(last_significant_coeff_y_prefix < (log2TrafoSize << 1));

	LastSignificantCoeffX = last_significant_coeff_x_prefix;
	if (last_significant_coeff_x_prefix > 3) {
		last_significant_coeff_x_suffix = rbsp->ae();
		assert(last_significant_coeff_x_suffix < (1 << ((last_significant_coeff_x_prefix>>1)-1)));

		LastSignificantCoeffX = (1 << ((last_significant_coeff_x_prefix>>1)-1))
							  * (2 + (last_significant_coeff_x_prefix & 1))
							  + last_significant_coeff_x_suffix;
	}
	LastSignificantCoeffY = last_significant_coeff_y_prefix;
	if (last_significant_coeff_y_prefix > 3) {
		last_significant_coeff_y_suffix = rbsp->ae();
		assert(last_significant_coeff_y_suffix < (1 << ((last_significant_coeff_y_prefix>>1)-1)));

		LastSignificantCoeffY = (1 << ((last_significant_coeff_y_prefix>>1)-1))
							  * (2 + (last_significant_coeff_y_prefix & 1))
							  + last_significant_coeff_y_suffix;
	}
	if (scanIdx == 2) {
		temp                  = LastSignificantCoeffX;
		LastSignificantCoeffX = LastSignificantCoeffY;
		LastSignificantCoeffY = temp;
	}

	lastScanPos  = 16;
	lastSubBlock = (1 << (log2TrafoSize-2)) * (1 << (log2TrafoSize-2)) - 1;
	do {
		if (lastScanPos == 0) {
			lastScanPos = 16;
			lastSubBlock--;
		}
		lastScanPos--;

		xS = ScanOrder[log2TrafoSize-2][scanIdx][lastSubBlock][0];
		yS = ScanOrder[log2TrafoSize-2][scanIdx][lastSubBlock][1];
		xC = (xS << 2) + ScanOrder[2][scanIdx][lastScanPos][0];
		yC = (yS << 2) + ScanOrder[2][scanIdx][lastScanPos][1];
	} while (xC != LastSignificantCoeffX || yC != LastSignificantCoeffY);

	for (int i = lastSubBlock; i >= 0; i--) {
		xS = ScanOrder[log2TrafoSize-2][scanIdx][i][0];
		yS = ScanOrder[log2TrafoSize-2][scanIdx][i][1];

		coded_sub_block_flag[xS][yS] = 0;
		if ((xS == 0 && yS == 0) ||
			(xS == (LastSignificantCoeffX>>2) && yS == (LastSignificantCoeffY>>2)))
			coded_sub_block_flag[xS][yS] = 1;
		inferSigCoeffFlag = 0;
		if (i < lastSubBlock && i > 0) {
			coded_sub_block_flag[xS][yS] = rbsp->ae();
			inferSigCoeffFlag = 1;
		}

		for (int n = (i == lastSubBlock ? lastScanPos-1 : 15); n >= 0; n--) {
			xC = (xS << 2) + ScanOrder[2][scanIdx][n][0];
			yC = (yS << 2) + ScanOrder[2][scanIdx][n][1];

			significant_coeff_flag[xC][yC] = 0;
			if ((xC == LastSignificantCoeffX && yC == LastSignificantCoeffY) ||
				((xC & 3) == 0 && (yC & 3) == 0 && inferSigCoeffFlag == 1))
				significant_coeff_flag[xC][yC] = 1;
			if (coded_sub_block_flag[xS][yS] && (n > 0 || ! inferSigCoeffFlag)) {
				significant_coeff_flag[xC][yC] = rbsp->ae();
				if (significant_coeff_flag[xC][yC])
					inferSigCoeffFlag = 0;
			}
		}

		firstSigScanPos = 16;
		lastSigScanPos = -1;
		numGreater1Flag = 0;
		lastGreater1ScanPos = -1;
		for (int n = 15; n >= 0; n--) {
			xC = (xS << 2) + ScanOrder[2][scanIdx][n][0];
			yC = (yS << 2) + ScanOrder[2][scanIdx][n][1];

			coeff_abs_level_greater1_flag[n] = 0;
			coeff_abs_level_greater2_flag[n] = 0;
			if (significant_coeff_flag[xC][yC]) {
				if (numGreater1Flag < 8) {
					coeff_abs_level_greater1_flag[n] = rbsp->ae();
					numGreater1Flag++;
					if (coeff_abs_level_greater1_flag[n] && lastGreater1ScanPos == -1)
						lastGreater1ScanPos = n;
				}
				if (lastSigScanPos == -1)
					lastSigScanPos = n;
				firstSigScanPos = n;
			}
		}

		signHidden = lastSigScanPos - firstSigScanPos > 3 && ! cu_transquant_bypass_flag;
		if (lastGreater1ScanPos != -1)
			coeff_abs_level_greater2_flag[lastGreater1ScanPos] = rbsp->ae();

		for (int n = 15; n >= 0; n--) {
			xC = (xS << 2) + ScanOrder[2][scanIdx][n][0];
			yC = (yS << 2) + ScanOrder[2][scanIdx][n][1];

			coeff_sign_flag[n] = 0;
			if (significant_coeff_flag[xC][yC] &&
				(! sign_data_hiding_flag || ! signHidden || n != firstSigScanPos))
				coeff_sign_flag[n] = rbsp->ae();
		}

		numSigCoeff = 0;
		sumAbsLevel = 0;
		for (int n = 15; n >= 0; n--) {
			xC = (xS << 2) + ScanOrder[2][scanIdx][n][0];
			yC = (yS << 2) + ScanOrder[2][scanIdx][n][1];

			coeff_abs_level_remaining[n] = 0;
			if (significant_coeff_flag[xC][yC]) {
				baseLevel = 1 + coeff_abs_level_greater1_flag[n] + coeff_abs_level_greater2_flag[n];
				if (baseLevel == numSigCoeff < 8 ? (n == lastGreater1ScanPos ? 3 : 2) : 1)
					coeff_abs_level_remaining[n] = rbsp->ae();
				TransCoeffLevel[x0][y0][cIdx][xC][yC] =
					(coeff_abs_level_remaining[n] + baseLevel) * (1 - 2 * coeff_sign_flag[n]);
				assert(TransCoeffLevel[x0][y0][cIdx][xC][yC] >= -32768);
				assert(TransCoeffLevel[x0][y0][cIdx][xC][yC] < 32768);

				if (sign_data_hiding_flag && signHidden) {
					sumAbsLevel += (coeff_abs_level_remaining[n] + baseLevel);
					if (n == firstSigScanPos && sumAbsLevel % 2 == 1)
						TransCoeffLevel[x0][y0][cIdx][xC][yC] = -TransCoeffLevel[x0][y0][cIdx][xC][yC];
				}
				numSigCoeff++;
			}
		}
	}
}
