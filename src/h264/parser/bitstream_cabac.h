#ifndef _CABAC_H_
#define _CABAC_H_


struct cabac_context_t {
    uint8_t pStateIdx; // index into state-table CP
    uint8_t valMPS;    // Least Probable Symbol 0/1 CP

    void init(int8_t m, int8_t n, uint8_t SliceQpY);
};


#define NUM_MB_TYPE_CTX        11
#define NUM_B8_TYPE_CTX         9
#define NUM_MV_RES_CTX         10
#define NUM_REF_NO_CTX          6
#define NUM_DELTA_QP_CTX        4
#define NUM_MB_AFF_CTX          4

#define NUM_BLOCK_TYPES        22

#define NUM_TRANSFORM_SIZE_CTX  3
#define NUM_IPR_CTX             2
#define NUM_CIPR_CTX            4
#define NUM_CBP_CTX             4
#define NUM_BCBP_CTX            4
#define NUM_MAP_CTX            15
#define NUM_LAST_CTX           15
#define NUM_ONE_CTX             5
#define NUM_ABS_CTX             5

struct cabac_contexts_t {
    cabac_context_t mb_type_contexts       [3][NUM_MB_TYPE_CTX];
    cabac_context_t b8_type_contexts       [2][NUM_B8_TYPE_CTX];
    cabac_context_t mv_res_contexts        [2][NUM_MV_RES_CTX];
    cabac_context_t ref_no_contexts        [2][NUM_REF_NO_CTX];
    cabac_context_t delta_qp_contexts         [NUM_DELTA_QP_CTX];
    cabac_context_t mb_aff_contexts           [NUM_MB_AFF_CTX];

    cabac_context_t transform_size_contexts   [NUM_TRANSFORM_SIZE_CTX];
    cabac_context_t ipr_contexts              [NUM_IPR_CTX];
    cabac_context_t cipr_contexts             [NUM_CIPR_CTX];
    cabac_context_t cbp_contexts           [3][NUM_CBP_CTX];
    cabac_context_t bcbp_contexts             [NUM_BLOCK_TYPES][NUM_BCBP_CTX];
    cabac_context_t map_contexts           [2][NUM_BLOCK_TYPES][NUM_MAP_CTX];
    cabac_context_t last_contexts          [2][NUM_BLOCK_TYPES][NUM_LAST_CTX];
    cabac_context_t one_contexts              [NUM_BLOCK_TYPES][NUM_ONE_CTX];
    cabac_context_t abs_contexts              [NUM_BLOCK_TYPES][NUM_ABS_CTX];

    void init(uint8_t slice_type, uint8_t cabac_init_idc, uint8_t SliceQpY);
};


#endif /* _CABAC_H_ */
