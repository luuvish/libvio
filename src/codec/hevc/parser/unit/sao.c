#include "rbsp.h"
#include "slice.h"

enum {
	SAO_TYPE_IDX_NONE_APPLIED = 0,
	SAO_TYPE_IDX_BAND_OFFSET  = 1,
	SAO_TYPE_IDX_EDGE_OFFSET  = 2
};

enum {
	SAO_EO_CLASS_1D_0   = 0,
	SAO_EO_CLASS_1D_90  = 1,
	SAO_EO_CLASS_1D_135 = 2,
	SAO_EO_CLASS_1D_45  = 3
}

void sao(rx, ry, rbsp_t *rbsp) {
	sao_merge_left_flag = 0;
	if (rx > 0) {
		leftCtbInSliceSeg = CtbAddrInSliceSeg > 0;
		leftCtbInTile = TileId[CtbAddrTS] == TileId[CtbAddrRStoTS[CtbAddrRS-1]];
		if (leftCtbInSliceSeg && leftCtbInTile)
			sao_merge_left_flag = rbsp->ae(AE_SAO_MERGE_FLAG, "sao_merge_left_flag");
	}

	sao_merge_up_flag = 0;
	if (ry > 0 && ! sao_merge_left_flag) {
		upCtbInSliceSeg = (CtbAddrRS - PicWidthInCtbsY) >= slice->slice_segment_address;
		upCtbInTile = TileId[CtbAddrTS] == TileId[CtbAddrRStoTS[CtbAddrRS-PicWidthInCtbsY]];
		if (upCtbInSliceSeg && upCtbInTile)
			sao_merge_up_flag = rbsp->ae(AE_SAO_MERGE_FLAG, "sao_merge_up_flag");
	}

	for (int cIdx = 0; cIdx < 3; cIdx++) {
		if ((slice->slice_sao_luma_flag   && cIdx == 0) ||
			(slice->slice_sao_chroma_flag && cIdx > 0)) {
			if (sao_merge_left_flag) {
				SaoTypeIdx[cIdx][rx][ry] = SaoTypeIdx[cIdx][rx-1][ry];

				for (int i = 0; i < 4; i++)
					sao_offset_abs[cIdx][rx][ry][i] = sao_offset_abs[cIdx][rx-1][ry][i];
				if (SaoTypeIdx[cIdx][rx][ry] == SAO_TYPE_IDX_BAND_OFFSET) {
					for (int i = 0; i < 4; i++)
						sao_offset_sign[cIdx][rx][ry][i] = sao_offset_sign[cIdx][rx-1][ry][i];
					sao_band_position[cIdx][rx][ry] = sao_band_position[cIdx][rx-1][ry];
				}
				if (SaoTypeIdx[cIdx][rx][ry] == SAO_TYPE_IDX_EDGE_OFFSET)
					SaoEoClass[cIdx][rx][ry] = SaoEoClass[cIdx][rx-1][ry];
			}
			else if (sao_merge_up_flag) {
				SaoTypeIdx[cIdx][rx][ry] = SaoTypeIdx[cIdx][rx][ry-1];

				for (int i = 0; i < 4; i++)
					sao_offset_abs[cIdx][rx][ry][i] = sao_offset_abs[cIdx][rx][ry-1][i];
				if (SaoTypeIdx[cIdx][rx][ry] == SAO_TYPE_IDX_BAND_OFFSET) {
					for (int i = 0; i < 4; i++)
						sao_offset_sign[cIdx][rx][ry][i] = sao_offset_sign[cIdx][rx][ry-1][i];
					sao_band_position[cIdx][rx][ry] = sao_band_position[cIdx][rx][ry-1];
				}
				if (SaoTypeIdx[cIdx][rx][ry] == SAO_TYPE_IDX_EDGE_OFFSET)
					SaoEoClass[cIdx][rx][ry] = SaoEoClass[cIdx][rx][ry-1];
			}
			else {
				if (cIdx == 0) {
					sao_type_idx_luma = rbsp->ae();
					SaoTypeIdx[0][rx][ry] = sao_type_idx_luma;
				}
				if (cIdx == 1) {
					sao_type_idx_chroma = rbsp->ae();
					SaoTypeIdx[1][rx][ry] = sao_type_idx_chroma;
					SaoTypeIdx[2][rx][ry] = sao_type_idx_chroma;
				}

				if (SaoTypeIdx[cIdx][rx][ry] != SAO_TYPE_IDX_NONE_APPLIED) {
					for (int i = 0; i < 4; i++)
						sao_offset_abs[cIdx][rx][ry][i] = rbsp->ae();
				}
				if (SaoTypeIdx[cIdx][rx][ry] == SAO_TYPE_IDX_BAND_OFFSET) {
					for (int i = 0; i < 4; i++) {
						if (sao_offset_abs[cIdx][rx][ry][i] != 0)
							sao_offset_sign[cIdx][rx][ry][i] = rbsp->ae();
					}
					sao_band_position[cIdx][rx][ry] = rbsp->ae();
				}
				if (SaoTypeIdx[cIdx][rx][ry] == SAO_TYPE_IDX_EDGE_OFFSET) {
					if (cIdx == 0) {
						sao_eo_class_luma = rbsp->ae();
						SaoEoClass[0][rx][ry] = sao_eo_class_luma;
					}
					if (cIdx == 1) {
						sao_eo_class_chroma = rbsp->ae();
						SaoEoClass[1][rx][ry] = sao_eo_class_chroma;
						SaoEoClass[2][rx][ry] = sao_eo_class_chroma;
					}
				}
			}

			int bitDepth = cIdx == 0 ? sps->BitDepthY : sps->BitDepthC;

			offsetSign = -1;
			if (SaoTypeIdx[cIdx][rx][ry] == 2)
				offsetSign = i > 1 ? -1 : 1;
			if (SaoTypeIdx[cIdx][rx][ry] == 0)
				offsetSign = 1;

			SaoOffsetVal[cIdx][rx][ry][0] = 0;
			for (int i = 0; i < 4; i++) {
				SaoOffsetVal[cIdx][rx][ry][i+1] =
					(offsetSign * sao_offset_abs[cIdx][rx][ry][i])
					 << (bitDepth - min(bitDepth,10));
			}
		}
	}
}
