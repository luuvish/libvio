#include "rbsp.h"
#include "pred_weight.h"

void pred_weight_table(pred_weight_t *weight, rbsp_t *rbsp) {
	weight->luma_log2_weight_denom = rbsp->ue("luma_log2_weight_denom");
	assert(weight->luma_log2_weight_denom < 8);
	if (sps->chroma_format_idc != 0) {
		int32_t delta_chroma_log2_weight_denom = rbsp->se("delta_chroma_log2_weight_denom");
		weight->ChromaLog2WeightDenom = weight->luma_log2_weight_denom + delta_chroma_log2_weight_denom;
		assert(weight->ChromaLog2WeightDenom >= 0 && weight->ChromaLog2WeightDenom < 8);
	}

	if (slice->slice_type != SLICE_TYPE_I) {
		for (int i = 0; i <= slice->num_ref_idx_l0_active_minus1; i++)
			weight->luma_weight_l0_flag[i] = rbsp->u(1, "luma_weight_lX_flag");
		if (sps->chroma_format_idc != 0) {
			for (int i = 0; i <= slice->num_ref_idx_l0_active_minus1; i++)
				weight->chroma_weight_l0_flag[i] = rbsp->u(1, "chroma_weight_lX_flag");
		}

		for (int i = 0; i <= slice->num_ref_idx_l0_active_minus1; i++) {
			weight->LumaWeightL0[i] = (1 << weight->luma_log2_weight_denom);
			weight->LumaOffsetL0[i] = 0;

			if (weight->luma_weight_l0_flag[i]) {
				int32_t delta_luma_weight_l0 = rbsp->se("delta_luma_weight_lX");
				int32_t luma_offset_l0       = rbsp->se("luma_offset_lX");
				assert(delta_luma_weight_l0 >= -128 && delta_luma_weight_l0 < 128);
				assert(luma_offset_l0       >= -128 && luma_offset_l0       < 128);

				weight->LumaWeightL0[i] += delta_luma_weight_l0;
				weight->LumaOffsetL0[i] = luma_offset_l0;
			}

			if (sps->chroma_format_idc != 0) {
				for (int j = 0; j < 2; j++) {
					weight->ChromaWeightL0[i][j] = (1 << weight->ChromaLog2WeightDenom);
					weight->ChromaOffsetL0[i][j] = 0;

					if (weight->chroma_weight_l0_flag[i]) {
						int32_t delta_chroma_weight_l0 = rbsp->se("delta_chroma_weight_lX");
						int32_t delta_chroma_offset_l0 = rbps->se("delta_chroma_offset_lX");
						assert(delta_chroma_weight_l0 >= -128 && delta_chroma_weight_l0 < 128);
						assert(delta_chroma_offset_l0 >= -512 && delta_chroma_offset_l0 < 512);

						weight->ChromaWeightL0[i][j] += delta_chroma_weight_l0;
						weight->ChromaOffsetL0[i][j] =
						 	clip(-128, 127,
								 delta_chroma_offset_l0
								 - (((128 * weight->ChromaWeightL0[i][j]) >> weight->ChromaLog2WeightDenom) - 128));
					}
				}
			}
		}
	}

	if (slice->slice_type == SLICE_TYPE_B) {
		for (int i = 0; i <= slice->num_ref_idx_l1_active_minus1; i++)
			weight->luma_weight_l1_flag[i] = rbsp->u(1, "luma_weight_lX_flag");
		if (sps->chroma_format_idc != 0) {
			for (int i = 0; i <= slice->num_ref_idx_l1_active_minus1; i++)
				weight->chroma_weight_l1_flag[i] = rbsp->u(1, "chroma_weight_lX_flag");
		}

		for (int i = 0; i <= slice->num_ref_idx_l1_active_minus1; i++) {
			weight->LumaWeightL1[i] = (1 << weight->luma_log2_weight_denom);
			weight->LumaOffsetL1[i] = 0;
			if (weight->luma_weight_l1_flag[i]) {
				int32_t delta_luma_weight_l1 = rbsp->se("delta_luma_weight_lX");
				int32_t luma_offset_l1       = rbsp->se("luma_offset_lX");
				assert(delta_luma_weight_l1 >= -128 && delta_luma_weight_l1 < 128);
				assert(luma_offset_l1       >= -128 && luma_offset_l1       < 128);

				weight->LumaWeightL1[i] += delta_luma_weight_l1;
				weight->LumaOffsetL1[i] = luma_offset_l1;
			}

			if (sps->chroma_format_idc != 0) {
				for (int j = 0; j < 2; j++) {
					weight->ChromaWeightL1[i][j] = (1 << weight->chroma_log2_weight_denom);
					weight->ChromaOffsetL1[i][j] = 0;

					if (weight->chroma_weight_l1_flag[i]) {
						int32_t delta_chroma_weight_l1 = rbsp->se("delta_chroma_weight_lX");
						int32_t delta_chroma_offset_l1 = rbps->se("delta_chroma_offset_lX");
						assert(delta_chroma_weight_l1 >= -128 && delta_chroma_weight_l1 < 128);
						assert(delta_chroma_offset_l1 >= -512 && delta_chroma_offset_l1 < 512);

						weight->ChromaWeightL1[i][j] += delta_chroma_weight_l1;
						weight->ChromaOffsetL1[i][j] =
						 	clip(-128, 127,
						 		 delta_chroma_offset_l1
								 - (((128 * weight->ChromaWeightL1[i][j]) >> weight->ChromaLog2WeightDenom) - 128));
					}
				}
			}
		}
	}

	int sumWeightL0Flags = 0;
	if (slice->slice_type != SLICE_TYPE_I) {
		for (int i = 0; i <= slice->num_ref_idx_l0_active_minus1; i++) {
			sumWeightL0Flags += weight->luma_weight_l0_flag[i];
			if (sps->chroma_format_idc != 0)
				sumWeightL0Flags += 2 * weight->chroma_weight_l0_flag[i];
		}
	}
	int sumWeightL1Flags = 0;
	if (slice->slice_type == SLICE_TYPE_B) {
		for (int i = 0; i <= slice->num_ref_idx_l1_active_minus1; i++) {
			sumWeightL1Flags += weight->luma_weight_l1_flag[i];
			if (sps->chroma_format_idc != 0)
				sumWeightL1Flags += 2 * weight->chroma_weight_l1_flag[i];
		}
	}
	if (slice->slice_type == SLICE_TYPE_P)
		assert(sumWeightL0Flags <= 24);
	if (slice->slice_type == SLICE_TYPE_B)
		assert(sumWeightL0Flags + sumWeightL1Flags <= 24);
}
