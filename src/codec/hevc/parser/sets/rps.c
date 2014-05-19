#include "rbsp.h"
#include "rps.h"

void short_term_ref_pic_set(rps_t *rps, int idx, rbsp_t *rbsp) {
	uint32_t inter_ref_pic_set_prediction_flag = rbsp->u(1, "inter_ref_pic_set_prediction_flag");
	if (idx == 0)
		assert(inter_ref_pic_set_prediction_flag == 0);

	if (inter_ref_pic_set_prediction_flag) {
		uint32_t delta_idx_minus1 = 0;
		if (idx == sps->num_short_term_ref_pic_sets)
			delta_idx_minus1 = rbsp->ue("delta_idx_minus1");
		assert(delta_idx_minus1 < idx);

		int RIdx = idx - (delta_idx_minus1 + 1);

		uint32_t delta_rps_sign       = rbsp->u(1, "delta_rps_sign");
		uint32_t abs_delta_rps_minus1 = rbsp->ue("abs_delta_rps_minus1");
		int DeltaRPS = (1 - 2 * delta_rps_sign) * (abs_delta_rps_minus1 + 1);

		uint32_t used_by_curr_pic_flag[MAX_SHORT_TERM_REF_PIC_SETS+1];
		uint32_t use_delta_flag       [MAX_SHORT_TERM_REF_PIC_SETS+1];
		for (int j = 0; j <= rps[RIdx].NumDeltaPocs; j++) {
			used_by_curr_pic_flag[j] = rbsp->u(1, "used_by_curr_pic_flag[j]");
			use_delta_flag       [j] = 1;
			if (! used_by_curr_pic_flag[j])
				use_delta_flag[j] = rbsp->u(1, "use_delta_flag[j]");
		}

		int i = 0;
		for (int j = rps[RIdx].NumPositivePics-1; j >= 0; j--) {
			int dPoc = rps[RIdx].DeltaPocS1[j] + DeltaRPS;
			if (dPoc < 0 && use_delta_flag[rps[RIdx].NumNegativePics + j]) {
				rps[idx].DeltaPocS0     [i  ] = dPoc;
				rps[idx].UsedByCurrPicS0[i++] = used_by_curr_pic_flag[rps[RIdx].NumNegativePics + j];
			}
		}
		if (DeltaRPS < 0 && use_delta_flag[rps[RIdx].NumDeltaPocs]) {
			rps[idx].DeltaPocS0     [i  ] = DeltaRPS;
			rps[idx].UsedByCurrPicS0[i++] = used_by_curr_pic_flag[rps[RIdx].NumDeltaPocs];
		}
		for (int j = 0; j < rps[RIdx].NumNegativePics; j++) {
			int dPoc = rps[RIdx].DeltaPocS0[j] + DeltaRPS;
			if (dPoc < 0 && use_delta_flag[j]) {
				rps[idx].DeltaPocS0     [i  ] = dPoc;
				rps[idx].UsedByCurrPicS0[i++] = used_by_curr_pic_flag[j];
			}
		}
		rps[idx].NumNegativePics = i;

		int i = 0;
		for (int j = rps[RIdx].NumNegativePics-1; j >= 0; j--) {
			int dPoc = rps[RIdx].DeltaPocS0[j] + DeltaRPS;
			if (dPoc > 0 && use_delta_flag[j]) {
				rps[idx].DeltaPocS1     [i  ] = dPoc;
				rps[idx].UsedByCurrPicS1[i++] = used_by_curr_pic_flag[j];
			}
		}
		if (DeltaRPS > 0 && use_delta_flag[rps[RIdx].NumDeltaPocs]) {
			rps[idx].DeltaPocS1     [i  ] = DeltaRPS;
			rps[idx].UsedByCurrPicS1[i++] = used_by_curr_pic_flag[rps[RIdx].NumDeltaPocs];
		}
		for (int j = 0; j < rps[RIdx].NumPositivePics; j++) {
			int dPoc = rps[RIdx].DeltaPocS1[j] + DeltaRPS;
			if (dPoc > 0 && use_delta_flag[rps[RIdx].NumNegativePics + j]) {
				rps[idx].DeltaPocS1     [i  ] = dPoc;
				rps[idx].UsedByCurrPicS1[i++] = used_by_curr_pic_flag[rps[RIdx].NumNegativePics + j];
			}
		}
		rps[idx].NumPositivePics = i;
	}
	else {
		uint32_t num_negative_pics = rbsp->ue("num_negative_pics");
		uint32_t num_positive_pics = rbsp->ue("num_positive_pics");
		assert(num_negative_pics <=
		 	   sps->sps_max_dec_pic_buffering[sps->sps_max_temporal_layers_minus1]);
		assert(num_positive_pics <=
		 	   sps->sps_max_dec_pic_buffering[sps->sps_max_temporal_layers_minus1]
		 	   - num_negative_pics);

		rps[idx].NumNegativePics = num_negative_pics;
		rps[idx].NumPositivePics = num_positive_pics;

		int prevDeltaPocS0 = 0;
		for (int i = 0; i < num_negative_pics; i++) {
			uint32_t delta_poc_s0_minus1      = rbsp->ue("delta_poc_s0_minus1[i]");
			uint32_t used_by_curr_pic_s0_flag = rbsp->u(1, "used_by_curr_pic_s0_flag[i]");
			assert(delta_poc_s0_minus1 < (1 << 15));

			prevDeltaPocS0 -= (delta_poc_s0_minus1 + 1);
			rps[idx].DeltaPocS0     [i] = prevDeltaPocS0;
			rps[idx].UsedByCurrPicS0[i] = used_by_curr_pic_s0_flag;
		}
		int prevDeltaPocS1 = 0;
		for (int i = 0; i < num_positive_pics; i++) {
			uint32_t delta_poc_s1_minus1      = rbsp->ue("delta_poc_s1_minus1[i]");
			uint32_t used_by_curr_pic_s1_flag = rbsp->u(1, "used_by_curr_pic_s1_flag[i]");
			assert(delta_poc_s1_minus1 < (1 << 15));

			prevDeltaPocS1 += (delta_poc_s1_minus1 + 1);
			rps[idx].DeltaPocS1     [i] = prevDeltaPocS1;
			rps[idx].UsedByCurrPicS1[i] = used_by_curr_pic_s1_flag;
		}
	}

	rps[idx].NumDeltaPocs = rps[idx].NumNegativePics + rps[idx].NumPositivePics;
}

void ref_pic_list_modification(slice_t *slice, rbsp_t *rbsp) {
	if (slice->slice_type != SLICE_TYPE_I) {
		rbsp->u(1, "ref_pic_list_modification_flag_l0");
		if (ref_pic_list_modification_flag_l0 && slice->NumPocTotalCurr > 1) {
			for (int i = 0; i <= slice->num_ref_idx_l0_active_minus1; i++) {
				uint32_t list_entry_l0 = rbsp->u(ceil(log2(slice->NumPocTotalCurr)), "list_entry_l0[i]");
				assert(list_entry_l0 < slice->NumPocTotalCurr);
			}
		}
	}
	if (slice->slice_type == SLICE_TYPE_B) {
		rbsp->u(1, "ref_pic_list_modification_flag_l1");
		if (ref_pic_list_modification_flag_l1 && slice->NumPocTotalCurr > 1) {
			for (int i = 0; i <= slice->num_ref_idx_l1_active_minus1; i++) {
				uint32_t list_entry_l1 = rbsp->u(ceil(log2(slice->NumPocTotalCurr)), "list_entry_l1[i]");
				assert(list_entry_l1 < slice->NumPocTotalCurr);
			}
		}
	}
}
