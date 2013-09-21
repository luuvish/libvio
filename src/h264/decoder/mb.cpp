#include "global.h"
#include "dpb.h"
#include "slice.h"
#include "macroblock.h"
#include "transform.h"
#include "intra_prediction.h"
#include "inter_prediction.h"


using namespace vio::h264;

intra_prediction_t intra_prediction;

//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

static void mb_pred_skip(mb_t* mb, ColorPlane curr_plane)
{
    slice_t* slice = mb->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;

    perform_mc(mb, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    transform.inverse_transform_inter(mb, curr_plane, slice->slice_type == SP_slice);
}

static void mb_pred_p_inter(mb_t* mb, ColorPlane curr_plane)
{
    slice_t* slice = mb->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;

    int step_h0 = BLOCK_STEP[mb->mb_type][0];
    int step_v0 = BLOCK_STEP[mb->mb_type][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = mb->SubMbType    [block8x8];
            int pred_dir = mb->SubMbPredMode[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    perform_mc(mb, curr_plane, dec_picture, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    transform.inverse_transform_inter(mb, curr_plane, slice->slice_type == SP_slice);
    if (mb->CodedBlockPatternLuma != 0 || mb->CodedBlockPatternChroma != 0)
        slice->is_reset_coeff = false;
}

static void mb_pred_b_direct(mb_t* mb, ColorPlane curr_plane)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;

    if (slice->direct_spatial_mv_pred_flag) {
        if (sps->direct_8x8_inference_flag)
            get_direct8x8spatial(mb, dec_picture);
        else
            get_direct4x4spatial(mb, dec_picture);
    } else {
        for (int block8x8 = 0; block8x8 < 4; block8x8++) {
            if (sps->direct_8x8_inference_flag)
                get_direct8x8temporal(mb, dec_picture, block8x8);
            else
                get_direct4x4temporal(mb, dec_picture, block8x8);
        }
    }

    int step_h0 = 2; //BLOCK_STEP[mb->mb_type][0];
    int step_v0 = 2; //BLOCK_STEP[mb->mb_type][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = mb->SubMbType    [block8x8];
            int pred_dir = mb->SubMbPredMode[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];
            if (mv_mode == 0) {
                step_h4 = sps->direct_8x8_inference_flag ? 2 : 1;
                step_v4 = sps->direct_8x8_inference_flag ? 2 : 1;
            }

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    perform_mc(mb, curr_plane, dec_picture, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    transform.inverse_transform_inter(mb, curr_plane, slice->slice_type == SP_slice); 
    if (mb->CodedBlockPatternLuma != 0 || mb->CodedBlockPatternChroma != 0)
        slice->is_reset_coeff = false;
}

static void mb_pred_b_inter8x8(mb_t* mb, ColorPlane curr_plane)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;

    int step_h0 = BLOCK_STEP[mb->mb_type][0];
    int step_v0 = BLOCK_STEP[mb->mb_type][1];

    int pred_dirs[4];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            pred_dirs[block8x8] = get_inter8x8(mb, dec_picture, block8x8);
        }
    }

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = mb->SubMbType[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];
            int pred_dir = pred_dirs[block8x8];
            if (mv_mode == 0) {
                step_h4 = sps->direct_8x8_inference_flag ? 2 : 1;
                step_v4 = sps->direct_8x8_inference_flag ? 2 : 1;
            }

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    perform_mc(mb, curr_plane, dec_picture, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    transform.inverse_transform_inter(mb, curr_plane, slice->slice_type == SP_slice);
    if (mb->CodedBlockPatternLuma != 0 || mb->CodedBlockPatternChroma != 0)
        slice->is_reset_coeff = false;
}


static void mb_pred_ipcm(mb_t* mb)
{
    int i, j, k;
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;

    for (i = 0; i < 16; ++i) {
        for (j = 0; j < 16 ; ++j)
            dec_picture->imgY[mb->mb.y * 16 + i][mb->mb.x * 16 + j] = (imgpel) slice->cof[0][i][j];
    }

    if (sps->chroma_format_idc != YUV400 && sps->separate_colour_plane_flag == 0) {
        for (k = 0; k < 2; ++k) {
            for (i = 0; i < sps->MbHeightC; ++i) {
                for (j = 0;j < sps->MbWidthC; ++j)
                    dec_picture->imgUV[k][mb->mb.y * sps->MbHeightC + i][mb->mb.x * sps->MbWidthC + j] = (imgpel) slice->cof[k + 1][i][j];
            }
        }
    }

    // for deblocking filter
    mb->update_qp(0);

    memset(mb->nz_coeff, 16, 3 * BLOCK_PIXELS * sizeof(byte));

    // for CABAC decoding of MB skip flag
    mb->mb_skip_flag = 0;

    //for deblocking filter CABAC
    mb->cbp_blks[0] = 0xFFFF;

    //For CABAC decoding of Dquant
    slice->last_dquant = 0;
    slice->is_reset_coeff = false;
    slice->is_reset_coeff_cr = false;
}

static void mb_pred_intra(mb_t* mb, ColorPlane curr_plane)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int blockoffset = mb->mb_type == I4MB ? 1 : mb->mb_type == I8MB ? 4 : 16;

    for (int block4x4 = 0; block4x4 < 16; block4x4 += blockoffset) {
        int ioff = ((block4x4 / 4) % 2) * 8 + ((block4x4 % 4) % 2) * 4;
        int joff = ((block4x4 / 4) / 2) * 8 + ((block4x4 % 4) / 2) * 4;

        if (mb->mb_type == I4MB)
            intra_prediction.intra_pred_4x4(mb, curr_plane, ioff, joff);
        else if (mb->mb_type == I8MB)
            intra_prediction.intra_pred_8x8(mb, curr_plane, ioff, joff);
        else if (mb->mb_type == I16MB)
            intra_prediction.intra_pred_16x16(mb, curr_plane, ioff, joff);

        if (mb->mb_type == I4MB)
            transform.inverse_transform_4x4(mb, curr_plane, ioff, joff);
        else if (mb->mb_type == I8MB)
            transform.inverse_transform_8x8(mb, curr_plane, ioff, joff);
        else if (mb->mb_type == I16MB)
            transform.inverse_transform_16x16(mb, curr_plane, ioff, joff);
    }

    if (mb->mb_type == I16MB || mb->CodedBlockPatternLuma != 0 || mb->CodedBlockPatternChroma != 0)
        slice->is_reset_coeff = false;

    if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
        intra_prediction.intra_pred_chroma(mb);

        for (int uv = 0; uv < 2; uv++)
            transform.inverse_transform_chroma(mb, (ColorPlane)(uv + 1));
        if (mb->CodedBlockPatternChroma)
            slice->is_reset_coeff_cr = false;
    }
}


static void set_chroma_vector(mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    VideoParameters* p_Vid = mb->p_Vid;

    if (!slice->MbaffFrameFlag) {
        if (!slice->field_pic_flag) {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < slice->listXsize[l]; k++)
                    slice->chroma_vector_adjustment[l][k] = 0; 
            }
        } else if (!slice->bottom_field_flag) {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < slice->listXsize[l]; k++) {
                    if (p_Vid->structure != slice->listX[l][k]->structure)
                        slice->chroma_vector_adjustment[l][k] = -2; 
                    else
                        slice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        } else {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < slice->listXsize[l]; k++) {
                    if (p_Vid->structure != slice->listX[l][k]->structure)
                        slice->chroma_vector_adjustment[l][k] = 2; 
                    else
                        slice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        }
    } else {
        int mb_nr = (mb->mbAddrX & 0x01);
        int k,l;  

        //////////////////////////
        // find out the correct list offsets
        if (mb->mb_field_decoding_flag) {
            int list_offset = slice->MbaffFrameFlag && mb->mb_field_decoding_flag ?
                              mb->mbAddrX % 2 ? 4 : 2 : 0;

            for (l = LIST_0 + list_offset; l <= (LIST_1 + list_offset); l++) {
                for (k = 0; k < slice->listXsize[l]; k++) {
                    if (mb_nr == 0 && slice->listX[l][k]->structure == BOTTOM_FIELD)
                        slice->chroma_vector_adjustment[l][k] = -2; 
                    else if (mb_nr == 1 && slice->listX[l][k]->structure == TOP_FIELD)
                        slice->chroma_vector_adjustment[l][k] = 2; 
                    else
                        slice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        } else {
            for (l = LIST_0; l <= (LIST_1); l++) {
                for(k = 0; k < slice->listXsize[l]; k++)
                    slice->chroma_vector_adjustment[l][k] = 0; 
            }
        }
    }

    int max_vmv_r;
    sps_t* sps = slice->active_sps;
    if (sps->level_idc <= 10)
        max_vmv_r = 64 * 4;
    else if (sps->level_idc <= 20)
        max_vmv_r = 128 * 4;
    else if (sps->level_idc <= 30)
        max_vmv_r = 256 * 4;
    else
        max_vmv_r = 512 * 4; // 512 pixels in quarter pixels

    slice->max_mb_vmv_r = slice->field_pic_flag || mb->mb_field_decoding_flag ? max_vmv_r >> 1 : max_vmv_r;
}

static void decode_one_component(mb_t* mb, ColorPlane curr_plane)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    VideoParameters* p_Vid = slice->p_Vid;

    if (sps->ChromaArrayType == 3) {
        if (!mb->is_intra_block) {
            storable_picture *vidref = p_Vid->no_reference_picture;
            int noref = slice->framepoc < p_Vid->recovery_poc;
            for (int j = 0; j < 6; j++) {
                for (int i = 0; i < slice->listXsize[j] ; i++) {
                    storable_picture *curr_ref = slice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = (curr_plane == PLANE_Y) ? curr_ref->imgY : curr_ref->imgUV[curr_plane-1];
                    }
                }
            }
        }
    }

    if (mb->mb_type == IPCM) {
        mb_pred_ipcm(mb);
        return;
    }

    if (mb->mb_type == I16MB || mb->mb_type == I4MB || mb->mb_type == I8MB) {
        mb_pred_intra(mb, curr_plane);
        return;
    }

    if (mb->mb_type == PSKIP) {
        if (slice->slice_type == B_SLICE)
            mb_pred_b_direct(mb, curr_plane);
        else
            mb_pred_skip(mb, curr_plane);
        return;
    }

    if (slice->slice_type == P_SLICE || slice->slice_type == SP_SLICE ||
        mb->mb_type == P16x16 || mb->mb_type == P16x8 || mb->mb_type == P8x16) {
        mb_pred_p_inter(mb, curr_plane);
        return;
    }

    mb_pred_b_inter8x8(mb, curr_plane);
}

void macroblock_t::decode()
{
    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;

    set_chroma_vector(this);

    decode_one_component(this, PLANE_Y);

    if (sps->ChromaArrayType == 3) {
        decode_one_component(this, PLANE_U);
        decode_one_component(this, PLANE_V);

        slice->is_reset_coeff    = false;
        slice->is_reset_coeff_cr = false;
    }
}
