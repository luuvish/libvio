/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : cabac.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _CABAC_H_
#define _CABAC_H_


namespace arrow {
namespace video {
namespace h264  {


struct cabac_context_t {
    uint8_t pStateIdx; // index into state-table CP
    uint8_t valMPS;    // Least Probable Symbol 0/1 CP

    void init(int8_t m, int8_t n, uint8_t SliceQpY);
};


#define NUM_SKIP_CTX            3
#define NUM_MB_AFF_CTX          3

#define NUM_MB_TYPE_CTX        11
#define NUM_B8_TYPE_CTX         4
#define NUM_TRANSFORM_SIZE_CTX  3
#define NUM_CBP_L_CTX           4
#define NUM_CBP_C_CTX           8
#define NUM_DELTA_QP_CTX        4
#define NUM_IPR_CTX             2
#define NUM_CIPR_CTX            4
#define NUM_REF_NO_CTX          6
#define NUM_MVD_CTX             7

#define NUM_BLOCK_TYPES        22
#define NUM_BCBP_CTX            4
#define NUM_MAP_CTX            15
#define NUM_LAST_CTX           15
#define NUM_ONE_CTX             5
#define NUM_ABS_CTX             5

struct cabac_contexts_t {
    cabac_context_t skip_contexts             [NUM_SKIP_CTX];
    cabac_context_t mb_aff_contexts           [NUM_MB_AFF_CTX];

    cabac_context_t mb_type_contexts          [NUM_MB_TYPE_CTX];
    cabac_context_t b8_type_contexts          [NUM_B8_TYPE_CTX];
    cabac_context_t transform_size_contexts   [NUM_TRANSFORM_SIZE_CTX];
    cabac_context_t cbp_l_contexts            [NUM_CBP_L_CTX];
    cabac_context_t cbp_c_contexts            [NUM_CBP_C_CTX];
    cabac_context_t delta_qp_contexts         [NUM_DELTA_QP_CTX];
    cabac_context_t ipr_contexts              [NUM_IPR_CTX];
    cabac_context_t cipr_contexts             [NUM_CIPR_CTX];
    cabac_context_t ref_no_contexts           [NUM_REF_NO_CTX];
    cabac_context_t mvd_x_contexts            [NUM_MVD_CTX];
    cabac_context_t mvd_y_contexts            [NUM_MVD_CTX];

    cabac_context_t bcbp_contexts             [NUM_BLOCK_TYPES][NUM_BCBP_CTX];
    cabac_context_t map_contexts           [2][NUM_BLOCK_TYPES][NUM_MAP_CTX];
    cabac_context_t last_contexts          [2][NUM_BLOCK_TYPES][NUM_LAST_CTX];
    cabac_context_t one_contexts              [NUM_BLOCK_TYPES][NUM_ONE_CTX];
    cabac_context_t abs_contexts              [NUM_BLOCK_TYPES][NUM_ABS_CTX];

    void init(uint8_t slice_type, uint8_t cabac_init_idc, uint8_t SliceQpY);
};


};
};
};


#endif /* _CABAC_H_ */
