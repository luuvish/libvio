#ifndef __PRED_WEIGHT_H__
#define __PRED_WEIGHT_H__

typedef struct pred_weight_t {
	uint32_t luma_log2_weight_denom;                        // ue(v)
	uint32_t chroma_log2_weight_denom;                      // se(v)

	uint32_t luma_weight_l0_flag  [MAX_REF_IDX];            // u(1)
	uint32_t luma_weight_l0       [MAX_REF_IDX];            // se(v)
	uint32_t luma_offset_l0       [MAX_REF_IDX];            // se(v)
	uint32_t chroma_weight_l0_flag[MAX_REF_IDX];            // u(1)
	uint32_t chroma_weight_l0     [MAX_REF_IDX][2];         // se(v)
	uint32_t chroma_offset_l0     [MAX_REF_IDX][2];         // se(v)

	uint32_t luma_weight_l1_flag  [MAX_REF_IDX];            // u(1)
	uint32_t luma_weight_l1       [MAX_REF_IDX];            // se(v)
	uint32_t luma_offset_l1       [MAX_REF_IDX];            // se(v)
	uint32_t chroma_weight_l1_flag[MAX_REF_IDX];            // u(1)
	uint32_t chroma_weight_l1     [MAX_REF_IDX][2];         // se(v)
	uint32_t chroma_offset_l1     [MAX_REF_IDX][2];         // se(v)
} pred_weight_t;

#endif /* __PRED_WEIGHT_H__ */
