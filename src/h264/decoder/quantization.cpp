#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "memalloc.h"
#include "image.h"
#include "neighbour.h"
#include "transform.h"
#include "quantization.h"


quantization_t quantization;


static int Flat_4x4_16[16] = {
    16, 16, 16, 16,
    16, 16, 16, 16,
    16, 16, 16, 16,
    16, 16, 16, 16
};

static int Flat_8x8_16[64] = {
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16,
    16, 16, 16, 16, 16, 16, 16, 16
};

// Table 7-3 Specification of default scaling lists Default_4x4_Intra and Default_4x4_Inter

static int Default_4x4_Intra[16] = {
     6, 13, 20, 28,
    13, 20, 28, 32,
    20, 28, 32, 37,
    28, 32, 37, 42
};

static int Default_4x4_Inter[16] = {
    10, 14, 20, 24,
    14, 20, 24, 27,
    20, 24, 27, 30,
    24, 27, 30, 34
};

// Table 7-4 Specification of default scaling lists Default_8x8_Intra and Default_8x8_Inter

static int Default_8x8_Intra[64] = {
     6, 10, 13, 16, 18, 23, 25, 27,
    10, 11, 16, 18, 23, 25, 27, 29,
    13, 16, 18, 23, 25, 27, 29, 31,
    16, 18, 23, 25, 27, 29, 31, 33,
    18, 23, 25, 27, 29, 31, 33, 36,
    23, 25, 27, 29, 31, 33, 36, 38,
    25, 27, 29, 31, 33, 36, 38, 40,
    27, 29, 31, 33, 36, 38, 40, 42
};

static int Default_8x8_Inter[64] = {
     9, 13, 15, 17, 19, 21, 22, 24,
    13, 13, 17, 19, 21, 22, 24, 25,
    15, 17, 19, 21, 22, 24, 25, 27,
    17, 19, 21, 22, 24, 25, 27, 28,
    19, 21, 22, 24, 25, 27, 28, 30,
    21, 22, 24, 25, 27, 28, 30, 32,
    22, 24, 25, 27, 28, 30, 32, 33,
    24, 25, 27, 28, 30, 32, 33, 35
};


// exported variables
static const int dequant_coef8[6][8][8] = {
    {{ 20, 19, 25, 19, 20, 19, 25, 19 },
     { 19, 18, 24, 18, 19, 18, 24, 18 },
     { 25, 24, 32, 24, 25, 24, 32, 24 },
     { 19, 18, 24, 18, 19, 18, 24, 18 },
     { 20, 19, 25, 19, 20, 19, 25, 19 },
     { 19, 18, 24, 18, 19, 18, 24, 18 },
     { 25, 24, 32, 24, 25, 24, 32, 24 },
     { 19, 18, 24, 18, 19, 18, 24, 18 }},
    {{ 22, 21, 28, 21, 22, 21, 28, 21 },
     { 21, 19, 26, 19, 21, 19, 26, 19 },
     { 28, 26, 35, 26, 28, 26, 35, 26 },
     { 21, 19, 26, 19, 21, 19, 26, 19 },
     { 22, 21, 28, 21, 22, 21, 28, 21 },
     { 21, 19, 26, 19, 21, 19, 26, 19 },
     { 28, 26, 35, 26, 28, 26, 35, 26 },
     { 21, 19, 26, 19, 21, 19, 26, 19 }},
    {{ 26, 24, 33, 24, 26, 24, 33, 24 },
     { 24, 23, 31, 23, 24, 23, 31, 23 },
     { 33, 31, 42, 31, 33, 31, 42, 31 },
     { 24, 23, 31, 23, 24, 23, 31, 23 },
     { 26, 24, 33, 24, 26, 24, 33, 24 },
     { 24, 23, 31, 23, 24, 23, 31, 23 },
     { 33, 31, 42, 31, 33, 31, 42, 31 },
     { 24, 23, 31, 23, 24, 23, 31, 23 }},
    {{ 28, 26, 35, 26, 28, 26, 35, 26 },
     { 26, 25, 33, 25, 26, 25, 33, 25 },
     { 35, 33, 45, 33, 35, 33, 45, 33 },
     { 26, 25, 33, 25, 26, 25, 33, 25 },
     { 28, 26, 35, 26, 28, 26, 35, 26 },
     { 26, 25, 33, 25, 26, 25, 33, 25 },
     { 35, 33, 45, 33, 35, 33, 45, 33 },
     { 26, 25, 33, 25, 26, 25, 33, 25 }},
    {{ 32, 30, 40, 30, 32, 30, 40, 30 },
     { 30, 28, 38, 28, 30, 28, 38, 28 },
     { 40, 38, 51, 38, 40, 38, 51, 38 },
     { 30, 28, 38, 28, 30, 28, 38, 28 },
     { 32, 30, 40, 30, 32, 30, 40, 30 },
     { 30, 28, 38, 28, 30, 28, 38, 28 },
     { 40, 38, 51, 38, 40, 38, 51, 38 },
     { 30, 28, 38, 28, 30, 28, 38, 28 }},
    {{ 36, 34, 46, 34, 36, 34, 46, 34 },
     { 34, 32, 43, 32, 34, 32, 43, 32 },
     { 46, 43, 58, 43, 46, 43, 58, 43 },
     { 34, 32, 43, 32, 34, 32, 43, 32 },
     { 36, 34, 46, 34, 36, 34, 46, 34 },
     { 34, 32, 43, 32, 34, 32, 43, 32 },
     { 46, 43, 58, 43, 46, 43, 58, 43 },
     { 34, 32, 43, 32, 34, 32, 43, 32 }}
};

static void set_dequant4x4(int (*InvLevelScale4x4)[4],  const int (*dequant)[4], int *qmatrix)
{
    for (int j = 0; j < 4; j++) {
        *(*InvLevelScale4x4      ) = *(*dequant      ) * *qmatrix++;
        *(*InvLevelScale4x4   + 1) = *(*dequant   + 1) * *qmatrix++;
        *(*InvLevelScale4x4   + 2) = *(*dequant   + 2) * *qmatrix++;
        *(*InvLevelScale4x4++ + 3) = *(*dequant++ + 3) * *qmatrix++;
    }
}

static void set_dequant8x8(int (*InvLevelScale8x8)[8],  const int (*dequant)[8], int *qmatrix)
{
    for (int j = 0; j < 8; j++) {
        *(*InvLevelScale8x8      ) = *(*dequant      ) * *qmatrix++;
        *(*InvLevelScale8x8   + 1) = *(*dequant   + 1) * *qmatrix++;
        *(*InvLevelScale8x8   + 2) = *(*dequant   + 2) * *qmatrix++;
        *(*InvLevelScale8x8   + 3) = *(*dequant   + 3) * *qmatrix++;
        *(*InvLevelScale8x8   + 4) = *(*dequant   + 4) * *qmatrix++;
        *(*InvLevelScale8x8   + 5) = *(*dequant   + 5) * *qmatrix++;
        *(*InvLevelScale8x8   + 6) = *(*dequant   + 6) * *qmatrix++;
        *(*InvLevelScale8x8++ + 7) = *(*dequant++ + 7) * *qmatrix++;
    }
}

static void CalculateQuant4x4Param(slice_t *currSlice)
{
    const int (*p_dequant_coef)[4][4] = dequant_coef;
    int  (*InvLevelScale4x4_Intra_0)[4][4] = currSlice->InvLevelScale4x4_Intra[0];
    int  (*InvLevelScale4x4_Intra_1)[4][4] = currSlice->InvLevelScale4x4_Intra[1];
    int  (*InvLevelScale4x4_Intra_2)[4][4] = currSlice->InvLevelScale4x4_Intra[2];
    int  (*InvLevelScale4x4_Inter_0)[4][4] = currSlice->InvLevelScale4x4_Inter[0];
    int  (*InvLevelScale4x4_Inter_1)[4][4] = currSlice->InvLevelScale4x4_Inter[1];
    int  (*InvLevelScale4x4_Inter_2)[4][4] = currSlice->InvLevelScale4x4_Inter[2];

    for (int k = 0; k < 6; k++) {
        set_dequant4x4(*InvLevelScale4x4_Intra_0++, *p_dequant_coef  , currSlice->qmatrix[0]);
        set_dequant4x4(*InvLevelScale4x4_Intra_1++, *p_dequant_coef  , currSlice->qmatrix[1]);
        set_dequant4x4(*InvLevelScale4x4_Intra_2++, *p_dequant_coef  , currSlice->qmatrix[2]);
        set_dequant4x4(*InvLevelScale4x4_Inter_0++, *p_dequant_coef  , currSlice->qmatrix[3]);
        set_dequant4x4(*InvLevelScale4x4_Inter_1++, *p_dequant_coef  , currSlice->qmatrix[4]);
        set_dequant4x4(*InvLevelScale4x4_Inter_2++, *p_dequant_coef++, currSlice->qmatrix[5]);
    }
}

static void CalculateQuant8x8Param(slice_t *currSlice)
{
    const int (*p_dequant_coef)[8][8] = dequant_coef8;
    int  (*InvLevelScale8x8_Intra_0)[8][8] = currSlice->InvLevelScale8x8_Intra[0];
    int  (*InvLevelScale8x8_Intra_1)[8][8] = currSlice->InvLevelScale8x8_Intra[1];
    int  (*InvLevelScale8x8_Intra_2)[8][8] = currSlice->InvLevelScale8x8_Intra[2];
    int  (*InvLevelScale8x8_Inter_0)[8][8] = currSlice->InvLevelScale8x8_Inter[0];
    int  (*InvLevelScale8x8_Inter_1)[8][8] = currSlice->InvLevelScale8x8_Inter[1];
    int  (*InvLevelScale8x8_Inter_2)[8][8] = currSlice->InvLevelScale8x8_Inter[2];

    for (int k = 0; k < 6; k++) {
        set_dequant8x8(*InvLevelScale8x8_Intra_0++, *p_dequant_coef  , currSlice->qmatrix[6]);
        set_dequant8x8(*InvLevelScale8x8_Inter_0++, *p_dequant_coef++, currSlice->qmatrix[7]);
    }

    p_dequant_coef = dequant_coef8;
    if (currSlice->active_sps->chroma_format_idc == 3) { // 4:4:4
        for (int k = 0; k < 6; k++) {
            set_dequant8x8(*InvLevelScale8x8_Intra_1++, *p_dequant_coef  , currSlice->qmatrix[8]);
            set_dequant8x8(*InvLevelScale8x8_Inter_1++, *p_dequant_coef  , currSlice->qmatrix[9]);
            set_dequant8x8(*InvLevelScale8x8_Intra_2++, *p_dequant_coef  , currSlice->qmatrix[10]);
            set_dequant8x8(*InvLevelScale8x8_Inter_2++, *p_dequant_coef++, currSlice->qmatrix[11]);
        }
    }
}


void assign_quant_params(slice_t *currSlice)
{
    sps_t *sps = currSlice->active_sps;
    pps_t *pps = currSlice->active_pps;

    if (!pps->pic_scaling_matrix_present_flag &&
        !sps->seq_scaling_matrix_present_flag) {
        for (int i = 0; i < 12; i++)
            currSlice->qmatrix[i] = (i < 6) ? Flat_4x4_16 : Flat_8x8_16;
    } else {
        int n_ScalingList = (sps->chroma_format_idc != YUV444) ? 8 : 12;
        if (sps->seq_scaling_matrix_present_flag) { // check sps first
            for (int i = 0; i < n_ScalingList; i++) {
                if (i < 6) {
                    if (!sps->seq_scaling_list_present_flag[i]) { // fall-back rule A
                        if (i == 0)
                            currSlice->qmatrix[i] = Default_4x4_Intra;
                        else if (i == 3)
                            currSlice->qmatrix[i] = Default_4x4_Inter;
                        else
                            currSlice->qmatrix[i] = currSlice->qmatrix[i-1];
                    } else {
                        if (sps->UseDefaultScalingMatrix4x4Flag[i])
                            currSlice->qmatrix[i] = (i < 3) ? Default_4x4_Intra : Default_4x4_Inter;
                        else
                            currSlice->qmatrix[i] = sps->ScalingList4x4[i];
                    }
                } else {
                    if (!sps->seq_scaling_list_present_flag[i]) { // fall-back rule A
                        if (i == 6)
                            currSlice->qmatrix[i] = Default_8x8_Intra;
                        else if (i == 7)
                            currSlice->qmatrix[i] = Default_8x8_Inter;
                        else
                            currSlice->qmatrix[i] = currSlice->qmatrix[i-2];
                    } else {
                        if (sps->UseDefaultScalingMatrix8x8Flag[i-6])
                            currSlice->qmatrix[i] = (i == 6 || i == 8 || i == 10) ? Default_8x8_Intra : Default_8x8_Inter;
                        else
                            currSlice->qmatrix[i] = sps->ScalingList8x8[i-6];
                    }
                }
            }
        }

        if (pps->pic_scaling_matrix_present_flag) { // then check pps
            for (int i = 0; i < n_ScalingList; i++) {
                if (i < 6) {
                    if (!pps->pic_scaling_list_present_flag[i]) { // fall-back rule B
                        if (i == 0) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                currSlice->qmatrix[i] = Default_4x4_Intra;
                        } else if (i == 3) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                currSlice->qmatrix[i] = Default_4x4_Inter;
                        } else
                            currSlice->qmatrix[i] = currSlice->qmatrix[i-1];
                    } else {
                        if (pps->UseDefaultScalingMatrix4x4Flag[i])
                            currSlice->qmatrix[i] = (i < 3) ? Default_4x4_Intra : Default_4x4_Inter;
                        else
                            currSlice->qmatrix[i] = pps->ScalingList4x4[i];
                    }
                } else {
                    if (!pps->pic_scaling_list_present_flag[i]) { // fall-back rule B
                        if (i == 6) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                currSlice->qmatrix[i] = Default_8x8_Intra;
                        } else if (i == 7) {
                            if (!sps->seq_scaling_matrix_present_flag)
                                currSlice->qmatrix[i] = Default_8x8_Inter;
                        } else  
                            currSlice->qmatrix[i] = currSlice->qmatrix[i-2];
                    } else {
                        if (pps->UseDefaultScalingMatrix8x8Flag[i-6])
                            currSlice->qmatrix[i] = (i == 6 || i == 8 || i == 10) ? Default_8x8_Intra : Default_8x8_Inter;
                        else
                            currSlice->qmatrix[i] = pps->ScalingList8x8[i-6];
                    }
                }
            }
        }
    }

    CalculateQuant4x4Param(currSlice);
    if (pps->transform_8x8_mode_flag)
        CalculateQuant8x8Param(currSlice);
}


static inline int rshift_rnd_sf(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}

int quantization_t::inverse_quantize(mb_t* mb, bool uv, ColorPlane pl, int i0, int j0, int levarr)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int qp_per = mb->qp_scaled[pl] / 6;
    int qp_rem = mb->qp_scaled[pl] % 6;
    int transform_pl = sps->separate_colour_plane_flag ? slice->colour_plane_id : pl;

    if (uv) {
        int (*InvLevelScale4x4)[4] = mb->is_intra_block ?
            slice->InvLevelScale4x4_Intra[pl][qp_rem] :
            slice->InvLevelScale4x4_Inter[pl][qp_rem];
        levarr = rshift_rnd_sf((levarr * InvLevelScale4x4[j0][i0]) << qp_per, 4);
    } else if (!mb->transform_size_8x8_flag) {
        int (*InvLevelScale4x4)[4] = mb->is_intra_block ?
            slice->InvLevelScale4x4_Intra[transform_pl][qp_rem] :
            slice->InvLevelScale4x4_Inter[transform_pl][qp_rem];
        levarr = rshift_rnd_sf((levarr * InvLevelScale4x4[j0][i0]) << qp_per, 4);
    } else {
        int (*InvLevelScale8x8)[8] = mb->is_intra_block ?
            slice->InvLevelScale8x8_Intra[transform_pl][qp_rem] :
            slice->InvLevelScale8x8_Inter[transform_pl][qp_rem];
        levarr = rshift_rnd_sf((levarr * InvLevelScale8x8[j0][i0]) << qp_per, 6);
    }

    return levarr;
}


// Table 8-13 Specification of mapping of idx to Cij for zig-zag and field scan

static const uint8_t ZIGZAG_SCAN_4x4[2][16][2] = {
    {{ 0, 0 }, { 1, 0 }, { 0, 1 }, { 0, 2 },
     { 1, 1 }, { 2, 0 }, { 3, 0 }, { 2, 1 },
     { 1, 2 }, { 0, 3 }, { 1, 3 }, { 2, 2 },
     { 3, 1 }, { 3, 2 }, { 2, 3 }, { 3, 3 }},
    {{ 0, 0 }, { 0, 1 }, { 1, 0 }, { 0, 2 },
     { 0, 3 }, { 1, 1 }, { 1, 2 }, { 1, 3 },
     { 2, 0 }, { 2, 1 }, { 2, 2 }, { 2, 3 },
     { 3, 0 }, { 3, 1 }, { 3, 2 }, { 3, 3 }}
};

// Table 8-14 Specification of mapping of idx to Cij for 8x8 zig-zag and 8x8 field scan

static const uint8_t ZIGZAG_SCAN_8x8[2][64][2] = {
    {{ 0, 0 }, { 1, 0 }, { 0, 1 }, { 0, 2 }, { 1, 1 }, { 2, 0 }, { 3, 0 }, { 2, 1 },
     { 1, 2 }, { 0, 3 }, { 0, 4 }, { 1, 3 }, { 2, 2 }, { 3, 1 }, { 4, 0 }, { 5, 0 },
     { 4, 1 }, { 3, 2 }, { 2, 3 }, { 1, 4 }, { 0, 5 }, { 0, 6 }, { 1, 5 }, { 2, 4 },
     { 3, 3 }, { 4, 2 }, { 5, 1 }, { 6, 0 }, { 7, 0 }, { 6, 1 }, { 5, 2 }, { 4, 3 },
     { 3, 4 }, { 2, 5 }, { 1, 6 }, { 0, 7 }, { 1, 7 }, { 2, 6 }, { 3, 5 }, { 4, 4 },
     { 5, 3 }, { 6, 2 }, { 7, 1 }, { 7, 2 }, { 6, 3 }, { 5, 4 }, { 4, 5 }, { 3, 6 },
     { 2, 7 }, { 3, 7 }, { 4, 6 }, { 5, 5 }, { 6, 4 }, { 7, 3 }, { 7, 4 }, { 6, 5 },
     { 5, 6 }, { 4, 7 }, { 5, 7 }, { 6, 6 }, { 7, 5 }, { 7, 6 }, { 6, 7 }, { 7, 7 }},
    {{ 0, 0 }, { 0, 1 }, { 0, 2 }, { 1, 0 }, { 1, 1 }, { 0, 3 }, { 0, 4 }, { 1, 2 },
     { 2, 0 }, { 1, 3 }, { 0, 5 }, { 0, 6 }, { 0, 7 }, { 1, 4 }, { 2, 1 }, { 3, 0 },
     { 2, 2 }, { 1, 5 }, { 1, 6 }, { 1, 7 }, { 2, 3 }, { 3, 1 }, { 4, 0 }, { 3, 2 },
     { 2, 4 }, { 2, 5 }, { 2, 6 }, { 2, 7 }, { 3, 3 }, { 4, 1 }, { 5, 0 }, { 4, 2 },
     { 3, 4 }, { 3, 5 }, { 3, 6 }, { 3, 7 }, { 4, 3 }, { 5, 1 }, { 6, 0 }, { 5, 2 },
     { 4, 4 }, { 4, 5 }, { 4, 6 }, { 4, 7 }, { 5, 3 }, { 6, 1 }, { 6, 2 }, { 5, 4 },
     { 5, 5 }, { 5, 6 }, { 5, 7 }, { 6, 3 }, { 7, 0 }, { 7, 1 }, { 6, 4 }, { 6, 5 },
     { 6, 6 }, { 6, 7 }, { 7, 2 }, { 7, 3 }, { 7, 4 }, { 7, 5 }, { 7, 6 }, { 7, 7 }}
};

void quantization_t::coeff_luma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    slice_t *slice = mb->p_Slice;

    bool field = slice->field_pic_flag || mb->mb_field_decoding_flag;
    const uint8_t (*zigzag_scan_4x4)[2] = ZIGZAG_SCAN_4x4[field];

    int i0 = zigzag_scan_4x4[runarr][0];
    int j0 = zigzag_scan_4x4[runarr][1];

    slice->cof[pl][j0 * 4][i0 * 4] = levarr;
}

void quantization_t::coeff_luma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    slice_t* slice = mb->p_Slice;

    bool field = slice->field_pic_flag || mb->mb_field_decoding_flag;
    const uint8_t (*zigzag_scan_4x4)[2] = ZIGZAG_SCAN_4x4[field];
    const uint8_t (*zigzag_scan_8x8)[2] = ZIGZAG_SCAN_8x8[field];

    int i0, j0;

    if (!mb->transform_size_8x8_flag) {
        i0 = zigzag_scan_4x4[runarr][0];
        j0 = zigzag_scan_4x4[runarr][1];
        mb->s_cbp[pl].blk |= ((int64_t)0x01 << (y0 * 4 + x0));
    } else {
        i0 = zigzag_scan_8x8[runarr][0];
        j0 = zigzag_scan_8x8[runarr][1];
        mb->s_cbp[pl].blk |= ((int64_t)0x33 << (y0 * 4 + x0));
    }

    if (!mb->TransformBypassModeFlag)
        levarr = this->inverse_quantize(mb, false, pl, i0, j0, levarr);
    slice->cof[pl][y0 * 4 + j0][x0 * 4 + i0] = levarr;
}

void quantization_t::coeff_chroma_dc(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int i0, j0;

    if (sps->ChromaArrayType == 1) {
        i0 = runarr % 2;
        j0 = runarr / 2;
    }
    if (sps->ChromaArrayType == 2) {
        i0 = ZIGZAG_SCAN_4x4[1][runarr][0];
        j0 = ZIGZAG_SCAN_4x4[1][runarr][1];
    }

    slice->cof[pl][j0 * 4][i0 * 4] = levarr;
}

void quantization_t::coeff_chroma_ac(mb_t* mb, ColorPlane pl, int x0, int y0, int runarr, int levarr)
{
    slice_t* slice = mb->p_Slice;

    bool field = slice->field_pic_flag || mb->mb_field_decoding_flag;
    const uint8_t (*zigzag_scan_4x4)[2] = ZIGZAG_SCAN_4x4[field];

    int i0 = zigzag_scan_4x4[runarr][0];
    int j0 = zigzag_scan_4x4[runarr][1];

    if (!mb->TransformBypassModeFlag)
        levarr = this->inverse_quantize(mb, true, pl, i0, j0, levarr);
    slice->cof[pl][y0 * 4 + j0][x0 * 4 + i0] = levarr;
}
