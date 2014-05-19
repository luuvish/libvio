
enum {
	PRED_L0 = 0,
	PRED_L1 = 1,
	PRED_BI = 2
};

void prediction_unit(x0, y0, nPbW, nPbH, rbsp_t *rbsp) {
	mvp_l0_flag   [x0][y0] = 0;
	mvp_l1_flag   [x0][y0] = 0;
	merge_idx     [x0][y0] = 0;
	inter_pred_idc[x0][y0] = PRED_L0;
	ref_idx_l0    [x0][y0] = 0;
	ref_idx_l1    [x0][y0] = 0;

	if (skip_flag[x0][y0]) {
		if (MaxNumMergeCand > 1)
			merge_idx[x0][y0] = rbsp->ae();
	}
	else {
		merge_flag[x0][y0] = rbsp->ae();
		if (merge_flag[x0][y0]) {
			if (MaxNumMergeCand > 1)
				merge_idx[x0][y0] = rbsp->ae();
		}
		else {
			if (slice_type == SLICE_B) {
				inter_pred_idc[x0][y0] = rbsp->ae();
				if (inter_pred_idc[x0][y0])
					assert(nPbW + nPbH != 12);
			}

			if (inter_pred_idc[x0][y0] != PRED_L1) {
				if (num_ref_idx_l0_active_minus1 > 0)
					ref_idx_l0[x0][y0] = rbsp->ae();
				mvd_coding(x0, y0, 0);
				mvp_l0_flag[x0][y0] = rbsp->ae();
			}

			if (inter_pred_idc[x0][y0] != PRED_L0) {
				if (num_ref_idx_l1_active_minus1 > 0)
					ref_idx_l1[x0][y0] = rbsp->ae();
				if (mvd_l1_zero_flag && inter_pred_idc[x0][y0] == PRED_BI) {
					mvd_l1[x0][y0][0] = 0;
					mvd_l1[x0][y0][1] = 0;
				}
				else
					mvd_coding(x0, y0, 1);
				mvp_l1_flag[x0][y0] = rbsp->ae();
			}
		}
	}
}

void mvd_coding(x0, y0, refList) {
	abs_mvd_greater0_flag[0] = rbsp->ae();
	abs_mvd_greater0_flag[1] = rbsp->ae();
	abs_mvd_greater1_flag[0] = 0;
	if (abs_mvd_greater0_flag[0])
		abs_mvd_greater1_flag[0] = rbsp->ae();
	abs_mvd_greater1_flag[1] = 0;
	if (abs_mvd_greater0_flag[1])
		abs_mvd_greater1_flag[1] = rbsp->ae();

	abs_mvd_minus2[0] = 0;
	mvd_sign_flag [0] = 0;
	if (abs_mvd_greater0_flag[0]) {
		abs_mvd_minus2[0] = -1;
		if (abs_mvd_greater1_flag[0])
			abs_mvd_minus2[0] = rbsp->ae();
		mvd_sign_flag[0] = rbsp->ae();
	}
	abs_mvd_minus2[1] = 0;
	mvd_sign_flag [1] = 0;
	if (abs_mvd_greater0_flag[1]) {
		abs_mvd_minus2[1] = -1;
		if (abs_mvd_greater1_flag[1])
			abs_mvd_minus2[1] = rbsp->ae();
		mvd_sign_flag[1] = rbsp->ae();
	}

	Mvd[0] = abs_mvd_greater0_flag[0] * (abs_mvd_minus2[0] + 2) * (1 - 2 * mvd_sign_flag[0]);
	Mvd[1] = abs_mvd_greater0_flag[1] * (abs_mvd_minus2[1] + 2) * (1 - 2 * mvd_sign_flag[1]);
	assert(Mvd[0] >= -(1<<15) && Mvd[0] < (1<<15));
	assert(Mvd[1] >= -(1<<15) && Mvd[1] < (1<<15));

	if (refList == 0) {
		MvdL0[x0][y0][0] = Mvd[0];
		MvdL0[x0][y0][1] = Mvd[1];
	}
	else {
		MvdL1[x0][y0][0] = Mvd[0];
		MvdL1[x0][y0][1] = Mvd[1];
	}
}
