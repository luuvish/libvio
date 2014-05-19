#include <stdint.h>

enum {
	SCALING_LIST_4x4   = 0,
	SCALING_LIST_8x8   = 1,
	SCALING_LIST_16x4  = 1,
	SCALING_LIST_4x16  = 1,
	SCALING_LIST_16x16 = 2,
	SCALING_LIST_32x8  = 2,
	SCALING_LIST_8x32  = 2,
	SCALING_LIST_32x32 = 3,
	SCALING_LIST_NUM   = 4
};

#define SCALING_LIST_32x32_MATRIX ( 2)
#define SCALING_MATRIX_NUM        ( 6)
#define SCALING_COEF_NUM          (64)
#define SCALING_COEF_START        ( 8)

typedef struct scaling_list_t {
	uint32_t dc_coef[SCALING_LIST_NUM-SCALING_LIST_16x16][SCALING_MATRIX_NUM];
	uint32_t coef[SCALING_LIST_NUM][SCALING_MATRIX_NUM][SCALING_COEF_NUM];
} scaling_list_t;


const SCALING_LIST_INTRA_4x4[16] = {
	16, 16, 16, 16,
	16, 16, 16, 16,
	16, 16, 16, 16,
	16, 16, 16, 16
};

const SCALING_LIST_INTER_4x4[16] = {
	16, 16, 16, 16,
	16, 16, 16, 16,
	16, 16, 16, 16,
	16, 16, 16, 16
};

const SCALING_LIST_INTRA_8x8[64] = {
	16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 17, 16, 17, 16, 17, 18,
	17, 18, 18, 17, 18, 21, 19, 20,
	21, 20, 19, 21, 24, 22, 22, 24,
	24, 22, 22, 24, 25, 25, 27, 30,
	27, 25, 25, 29, 31, 35, 35, 31,
	29, 36, 41, 44, 41, 36, 47, 54,
	54, 47, 65, 70, 65, 88, 88, 115
};

const SCALING_LIST_INTER_8x8[64] = {
	16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 17, 17, 17, 17, 17, 18,
	18, 18, 18, 18, 18, 20, 20, 20,
	20, 20, 20, 20, 24, 24, 24, 24,
	24, 24, 24, 24, 25, 25, 25, 25,
	25, 25, 25, 28, 28, 28, 28, 28,
	28, 33, 33, 33, 33, 33, 41, 41,
	41, 41, 54, 54, 54, 71, 71, 91
};

const SCALING_LIST[3][6][64] = {
	{
	 	SCALING_LIST_INTRA_4x4, SCALING_LIST_INTRA_4x4, SCALING_LIST_INTRA_4x4,
		SCALING_LIST_INTER_4x4, SCALING_LIST_INTER_4x4, SCALING_LIST_INTER_4x4
	},
	{
		SCALING_LIST_INTRA_8x8, SCALING_LIST_INTRA_8x8, SCALING_LIST_INTRA_8x8,
		SCALING_LIST_INTER_8x8, SCALING_LIST_INTER_8x8, SCALING_LIST_INTER_8x8
	},
	{
		SCALING_LIST_INTRA_8x8, SCALING_LIST_INTER_8x8,
	}
};


void scaling_list_data(scaling_list_t *scaling, rbsp_t *rbsp) {
	for (int sizeId = 0; sizeId < SCALING_LIST_NUM; sizeId++) {
		for (int matrixId = 0; matrixId < (sizeId == SCALING_LIST_32x32 ? SCALING_LIST_32x32_MATRIX : SCALING_MATRIX_NUM); matrixId++) {
			uint32_t scaling_list_pred_mode_flag = rbsp->u(1, "scaling_list_pred_mode_flag");
			if (! scaling_list_pred_mode_flag) {
				uint32_t scaling_list_pred_matrix_id_delta = rbsp->ue("scaling_list_pred_matrix_id_delta");
				assert(scaling_list_pred_matrix_id_delta <= matrixId);
				int refMatrixId = matrixId - scaling_list_pred_matrix_id_delta;
				for (int i = 0; i < SCALING_COEF_NUM; i++)
					scaling->coef[sizeId][matrixId][i] = scaling->coef[sizeId][refMatrixId][i];
			}
			else {
				int nextCoef = SCALING_COEF_START;
				int coefNum = min(SCALING_COEF_NUM, (1 << (4 + (sizeId<<1))));

				if (sizeId > SCALING_LIST_8x8) {
					int32_t scaling_list_dc_coef_minus8 = rbsp->se("scaling_list_dc_coef_minus8");
					assert(scaling_list_dc_coef_minus8 >= -7 && scaling_list_dc_coef_minus8 < 248);
					nextCoef = scaling->dc_coef[sizeId-SCALING_LIST_16x16][matrixId] = scaling_list_dc_coef_minus8 + 8;
				}

				for (int i = 0; i < coefNum; i++) {
					int32_t scaling_list_delta_coef = rbsp->se("scaling_list_delta_coef");
					assert(scaling_list_delta_coef >= -128 && scaling_list_delta_coef < 128);
					nextCoef = (nextCoef + scaling_list_delta_coef + 256) % 256;
					assert(nextCoef > 0);
					scaling->coef[sizeId][matrixId][i] = nextCoef;
				}
			}
		}
	}
}

void scaling_list_parameter(scaling_list *scaling) {
	for (int matrixId = 0; matrixId < 6; matrixId++) {
		for (int i = 0; i < 16; i++) {
			int x = sps->ScanOrder[2][SCAN_IDX_DIAG][i][SCAN_COMP_HOR];
			int y = sps->ScanOrder[2][SCAN_IDX_DIAG][i][SCAN_COMP_VER];
			pps->ScalingFactor[SCALING_LIST_4x4][matrixId][x][y] = scaling->coef[SCALING_LIST_4x4][matrixId][i];
		}
	}

	for (int matrixId = 0; matrixId < 6; matrixId++) {
		for (int i = 0; i < 64; i++) {
			int x = sps->ScanOrder[3][SCAN_IDX_DIAG][i][SCAN_COMP_HOR];
			int y = sps->ScanOrder[3][SCAN_IDX_DIAG][i][SCAN_COMP_VER];
			pps->ScalingFactor[SCALING_LIST_8x8][matrixId][x][y] = scaling->coef[SCALING_LIST_8x8][matrixId][i];
		}
	}

	for (int matrixId = 0; matrixId < 6; matrixId++) {
		for (int i = 0; i < 64; i++) {
			int x = sps->ScanOrder[3][SCAN_IDX_DIAG][i][SCAN_COMP_HOR];
			int y = sps->ScanOrder[3][SCAN_IDX_DIAG][i][SCAN_COMP_VER];
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 2; k++)
					pps->ScalingFactor[SCALING_LIST_16x16][matrixId][x*2+k][y*2+j] = scaling->coef[SCALING_LIST_16x16][matrixId][i];
			}
		}
		pps->ScalingFactor[SCALING_LIST_16x16][matrixId][0][0] = scaling->dc_coef[SCALING_LIST_16x16-SCALING_LIST_16x16][matrixId];
	}

	for (int matrixId = 0; matrixId < 2; matrixId++) {
		for (int i = 0; i < 64; i++) {
			int x = sps->ScanOrder[3][SCAN_IDX_DIAG][i][SCAN_COMP_HOR];
			int y = sps->ScanOrder[3][SCAN_IDX_DIAG][i][SCAN_COMP_VER];
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 4; k++)
					pps->ScalingFactor[SCALING_LIST_32x32][matrixId][x*4+k][y*4+j] = scaling->coef[SCALING_LIST_32x32][matrixId][i];
			}
		}
		pps->ScalingFactor[SCALING_LIST_32x32][matrixId][0][0] = scaling->dc_coef[SCALING_LIST_32x32-SCALING_LIST_16x16][matrixId];
	}
}
