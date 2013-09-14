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

static void mb_pred_skip(mb_t *currMB, ColorPlane curr_plane)
{
    slice_t *currSlice = currMB->p_Slice;
    StorablePicture *dec_picture = currSlice->dec_picture;

    perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    transform.inverse_transform_inter(currMB, curr_plane, currSlice->slice_type == SP_slice);
}

static void mb_pred_p_inter(mb_t *currMB, ColorPlane curr_plane)
{
    slice_t *currSlice = currMB->p_Slice;
    StorablePicture *dec_picture = currSlice->dec_picture;

    int partmode = (currMB->mb_type == P8x8 ? 4 : currMB->mb_type);
    int step_h0  = BLOCK_STEP[partmode][0];
    int step_v0  = BLOCK_STEP[partmode][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = currMB->b8mode[block8x8];
            int pred_dir = currMB->b8pdir[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    transform.inverse_transform_inter(currMB, curr_plane, currSlice->slice_type == SP_slice);
    if (currMB->CodedBlockPatternLuma != 0 || currMB->CodedBlockPatternChroma != 0)
        currSlice->is_reset_coeff = false;
}

static void mb_pred_b_direct(mb_t *currMB, ColorPlane curr_plane)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    StorablePicture *dec_picture = currSlice->dec_picture;

    if (currSlice->direct_spatial_mv_pred_flag) {
        if (sps->direct_8x8_inference_flag)
            get_direct8x8spatial(currMB, dec_picture);
        else
            get_direct4x4spatial(currMB, dec_picture);
    } else {
        for (int block8x8 = 0; block8x8 < 4; block8x8++) {
            if (sps->direct_8x8_inference_flag)
                get_direct8x8temporal(currMB, dec_picture, block8x8);
            else
                get_direct4x4temporal(currMB, dec_picture, block8x8);
        }
    }

    //int partmode = (currMB->mb_type == P8x8 ? 4 : currMB->mb_type);
    int step_h0 = 2; //BLOCK_STEP[partmode][0];
    int step_v0 = 2; //BLOCK_STEP[partmode][1];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = currMB->b8mode[block8x8];
            int pred_dir = currMB->b8pdir[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];
            if (mv_mode == 0) {
                step_h4 = sps->direct_8x8_inference_flag ? 2 : 1;
                step_v4 = sps->direct_8x8_inference_flag ? 2 : 1;
            }

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    transform.inverse_transform_inter(currMB, curr_plane, currSlice->slice_type == SP_slice); 
    if (currMB->CodedBlockPatternLuma != 0 || currMB->CodedBlockPatternChroma != 0)
        currSlice->is_reset_coeff = false;
}

static void mb_pred_b_inter8x8(mb_t *currMB, ColorPlane curr_plane)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    StorablePicture *dec_picture = currSlice->dec_picture;

    int partmode = (currMB->mb_type == P8x8 ? 4 : currMB->mb_type);
    int step_h0  = BLOCK_STEP[partmode][0];
    int step_v0  = BLOCK_STEP[partmode][1];

    int pred_dirs[4];

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            pred_dirs[block8x8] = get_inter8x8(currMB, dec_picture, block8x8);
        }
    }

    for (int j0 = 0; j0 < 4; j0 += step_v0) {
        for (int i0 = 0; i0 < 4; i0 += step_h0) {
            int block8x8 = 2 * (j0 >> 1) + (i0 >> 1);
            int mv_mode  = currMB->b8mode[block8x8];
            int step_h4  = BLOCK_STEP[mv_mode][0];
            int step_v4  = BLOCK_STEP[mv_mode][1];
            int pred_dir = pred_dirs[block8x8];
            if (mv_mode == 0) {
                step_h4 = sps->direct_8x8_inference_flag ? 2 : 1;
                step_v4 = sps->direct_8x8_inference_flag ? 2 : 1;
            }

            for (int j = j0; j < j0 + step_v0; j += step_v4) {
                for (int i = i0; i < i0 + step_h0; i += step_h4)
                    perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, step_h4 * 4, step_v4 * 4);
            }
        }
    }

    transform.inverse_transform_inter(currMB, curr_plane, currSlice->slice_type == SP_slice);
    if (currMB->CodedBlockPatternLuma != 0 || currMB->CodedBlockPatternChroma != 0)
        currSlice->is_reset_coeff = false;
}


static void mb_pred_ipcm(mb_t *currMB)
{
    int i, j, k;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    StorablePicture *dec_picture = currSlice->dec_picture;

    for (i = 0; i < 16; ++i) {
        for (j = 0; j < 16 ; ++j)
            dec_picture->imgY[currMB->mb.y * 16 + i][currMB->mb.x * 16 + j] = (imgpel) currSlice->cof[0][i][j];
    }

    if (sps->chroma_format_idc != YUV400 && sps->separate_colour_plane_flag == 0) {
        for (k = 0; k < 2; ++k) {
            for (i = 0; i < sps->MbHeightC; ++i) {
                for (j = 0;j < sps->MbWidthC; ++j)
                    dec_picture->imgUV[k][currMB->mb.y * sps->MbHeightC + i][currMB->mb.x * sps->MbWidthC + j] = (imgpel) currSlice->cof[k + 1][i][j];
            }
        }
    }

    // for deblocking filter
    currMB->update_qp(0);

    memset(currMB->nz_coeff, 16, 3 * BLOCK_PIXELS * sizeof(byte));

    // for CABAC decoding of MB skip flag
    currMB->mb_skip_flag = 0;

    //for deblocking filter CABAC
    currMB->cbp_blks[0] = 0xFFFF;

    //For CABAC decoding of Dquant
    currSlice->last_dquant = 0;
    currSlice->is_reset_coeff = false;
    currSlice->is_reset_coeff_cr = false;
}

static void mb_pred_intra(mb_t *currMB, ColorPlane curr_plane)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int blockoffset = currMB->mb_type == I4MB ? 1
                    : currMB->mb_type == I8MB ? 4 : 16;

    for (int block4x4 = 0; block4x4 < 16; block4x4 += blockoffset) {
        int ioff = ((block4x4 / 4) % 2) * 8 + ((block4x4 % 4) % 2) * 4;
        int joff = ((block4x4 / 4) / 2) * 8 + ((block4x4 % 4) / 2) * 4;

        if (currMB->mb_type == I4MB)
            intra_prediction.intra_pred_4x4(currMB, curr_plane, ioff, joff);
        else if (currMB->mb_type == I8MB)
            intra_prediction.intra_pred_8x8(currMB, curr_plane, ioff, joff);
        else if (currMB->mb_type == I16MB)
            intra_prediction.intra_pred_16x16(currMB, curr_plane, ioff, joff);

        if (currMB->mb_type == I4MB)
            transform.inverse_transform_4x4(currMB, curr_plane, ioff, joff);
        else if (currMB->mb_type == I8MB)
            transform.inverse_transform_8x8(currMB, curr_plane, ioff, joff);
        else if (currMB->mb_type == I16MB)
            transform.inverse_transform_16x16(currMB, curr_plane, ioff, joff);
    }

    if (currMB->mb_type == I16MB || currMB->CodedBlockPatternLuma != 0 || currMB->CodedBlockPatternChroma != 0)
        currSlice->is_reset_coeff = false;

    if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
        intra_prediction.intra_pred_chroma(currMB);

        for (int uv = 0; uv < 2; uv++)
            transform.inverse_transform_chroma(currMB, (ColorPlane)(uv + 1));
        if (currMB->CodedBlockPatternChroma)
            currSlice->is_reset_coeff_cr = false;
    }
}


static void set_chroma_vector(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    if (!currSlice->MbaffFrameFlag) {
        if (!currSlice->field_pic_flag) {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++)
                    currSlice->chroma_vector_adjustment[l][k] = 0; 
            }
        } else if (!currSlice->bottom_field_flag) {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++) {
                    if (p_Vid->structure != currSlice->listX[l][k]->structure)
                        currSlice->chroma_vector_adjustment[l][k] = -2; 
                    else
                        currSlice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        } else {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++) {
                    if (p_Vid->structure != currSlice->listX[l][k]->structure)
                        currSlice->chroma_vector_adjustment[l][k] = 2; 
                    else
                        currSlice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        }
    } else {
        int mb_nr = (currMB->mbAddrX & 0x01);
        int k,l;  

        //////////////////////////
        // find out the correct list offsets
        if (currMB->mb_field_decoding_flag) {
            int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                              currMB->mbAddrX % 2 ? 4 : 2 : 0;

            for (l = LIST_0 + list_offset; l <= (LIST_1 + list_offset); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++) {
                    if (mb_nr == 0 && currSlice->listX[l][k]->structure == BOTTOM_FIELD)
                        currSlice->chroma_vector_adjustment[l][k] = -2; 
                    else if (mb_nr == 1 && currSlice->listX[l][k]->structure == TOP_FIELD)
                        currSlice->chroma_vector_adjustment[l][k] = 2; 
                    else
                        currSlice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        } else {
            for (l = LIST_0; l <= (LIST_1); l++) {
                for(k = 0; k < currSlice->listXsize[l]; k++)
                    currSlice->chroma_vector_adjustment[l][k] = 0; 
            }
        }
    }

    int max_vmv_r;
    sps_t *sps = currSlice->active_sps;
    if (sps->level_idc <= 10)
        max_vmv_r = 64 * 4;
    else if (sps->level_idc <= 20)
        max_vmv_r = 128 * 4;
    else if (sps->level_idc <= 30)
        max_vmv_r = 256 * 4;
    else
        max_vmv_r = 512 * 4; // 512 pixels in quarter pixels

    currSlice->max_mb_vmv_r = currSlice->field_pic_flag || currMB->mb_field_decoding_flag ? max_vmv_r >> 1 : max_vmv_r;
}

static void decode_one_component(mb_t *currMB, ColorPlane curr_plane)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currSlice->p_Vid;

    if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
        if (!currMB->is_intra_block) {
            StorablePicture *vidref = p_Vid->no_reference_picture;
            int noref = currSlice->framepoc < p_Vid->recovery_poc;
            for (int j = 0; j < 6; j++) {
                for (int i = 0; i < currSlice->listXsize[j] ; i++) {
                    StorablePicture *curr_ref = currSlice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = (curr_plane == PLANE_Y) ? curr_ref->imgY : curr_ref->imgUV[curr_plane-1];
                    }
                }
            }
        }
    }

    if (currMB->mb_type == IPCM) {
        mb_pred_ipcm(currMB);
        return;
    }

    if (currMB->mb_type == I16MB || currMB->mb_type == I4MB || currMB->mb_type == I8MB) {
        mb_pred_intra(currMB, curr_plane);
        return;
    }

    if (currMB->mb_type == PSKIP) {
        if (currSlice->slice_type == B_SLICE)
            mb_pred_b_direct(currMB, curr_plane);
        else
            mb_pred_skip(currMB, curr_plane);
        return;
    }

    if (currSlice->slice_type == P_SLICE || currSlice->slice_type == SP_SLICE ||
        currMB->mb_type == P16x16 || currMB->mb_type == P16x8 || currMB->mb_type == P8x16) {
        mb_pred_p_inter(currMB, curr_plane);
        return;
    }

    mb_pred_b_inter8x8(currMB, curr_plane);
}

void macroblock_t::decode()
{
    slice_t *slice = this->p_Slice;
    sps_t *sps = slice->active_sps;

    set_chroma_vector(this);

    decode_one_component(this, PLANE_Y);

    if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
        decode_one_component(this, PLANE_U);
        decode_one_component(this, PLANE_V);

        slice->is_reset_coeff    = false;
        slice->is_reset_coeff_cr = false;
    }
}
