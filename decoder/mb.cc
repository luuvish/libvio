
#include "global.h"
#include "dpb.h"
#include "slice.h"
#include "macroblock.h"
#include "transform.h"
#include "intra_prediction.h"
#include "inter_prediction.h"
#include "mb.h"


// number of intra prediction modes
#define NO_INTRA_PMODE  9


static void set_chroma_vector(Macroblock *currMB)
{
    Slice *currSlice = currMB->p_Slice;
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


//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

static int mb_pred_skip(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
        copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
        copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
    }
    return 1;
}

static int mb_pred_sp_skip(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    iTransform(currMB, curr_plane, 1);
    return 1;
}

static int mb_pred_p_inter(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;

    set_chroma_vector(currMB);

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

    iTransform(currMB, curr_plane, currSlice->slice_type == SP_SLICE);

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

static int mb_pred_b_d8x8temporal(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    set_chroma_vector(currMB);

    for (int block8x8 = 0; block8x8 < 4; block8x8++) {
        if (sps->direct_8x8_inference_flag)
            get_direct8x8temporal(currMB, dec_picture, block8x8);
        else
            get_direct4x4temporal(currMB, dec_picture, block8x8);
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

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }
    return 1;
}

static int mb_pred_b_d8x8spatial(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    set_chroma_vector(currMB);

    if (sps->direct_8x8_inference_flag)
        get_direct8x8spatial(currMB, dec_picture);
    else
        get_direct4x4spatial(currMB, dec_picture);

    //int partmode = (currMB->mb_type == P8x8 ? 4 : currMB->mb_type);
    int step_h0  = 2; //BLOCK_STEP[partmode][0];
    int step_v0  = 2; //BLOCK_STEP[partmode][1];

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

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if (sps->chroma_format_idc != YUV400 && sps->chroma_format_idc != YUV444) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, sps->MbWidthC, sps->MbHeightC);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }

    return 1;
}

static int mb_pred_b_inter8x8(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    set_chroma_vector(currMB);

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

    iTransform(currMB, curr_plane, 0);
    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}


static int mb_pred_ipcm(Macroblock *currMB)
{
    int i, j, k;
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;
    StorablePicture *dec_picture = currSlice->dec_picture;

    for (i = 0; i < 16; ++i) {
        for (j = 0;j < 16 ; ++j)
            dec_picture->imgY[currMB->pix_y + i][currMB->pix_x + j] = (imgpel) currSlice->cof[0][i][j];
    }

    if (sps->chroma_format_idc != YUV400 && sps->separate_colour_plane_flag == 0) {
        for (k = 0; k < 2; ++k) {
            for (i = 0; i < sps->MbHeightC; ++i) {
                for (j = 0;j < sps->MbWidthC; ++j)
                    dec_picture->imgUV[k][currMB->pix_c_y + i][currMB->pix_c_x + j] = (imgpel) currSlice->cof[k + 1][i][j];
            }
        }
    }

    // for deblocking filter
    update_qp(currMB, 0);

    // for CAVLC: Set the nz_coeff to 16.
    // These parameters are to be used in CAVLC decoding of neighbour blocks  
    memset(p_Vid->nz_coeff[currMB->mbAddrX][0][0], 16, 3 * BLOCK_PIXELS * sizeof(byte));

    // for CABAC decoding of MB skip flag
    currMB->mb_skip_flag = 0;

    //for deblocking filter CABAC
    currMB->s_cbp[0].blk = 0xFFFF;

    //For CABAC decoding of Dquant
    currSlice->last_dquant = 0;
    currSlice->is_reset_coeff = FALSE;
    currSlice->is_reset_coeff_cr = FALSE;
    return 1;
}

static void intra_cr_decoding(Macroblock *currMB, int yuv)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    StorablePicture *dec_picture = currSlice->dec_picture;
    imgpel **curUV;
    int uv;
    int b8,b4;
    int ioff, joff;
    int i,j;

    intra_pred_chroma(currMB);// last argument is ignored, computes needed data for both uv channels

    void (*itrans_4x4)(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
    itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : itrans4x4_ls;

    int num_uv_blocks;
    if (sps->chroma_format_idc != YUV400)
        num_uv_blocks = (((1 << sps->chroma_format_idc) & (~0x1)) >> 1);
    else
        num_uv_blocks = 0;

    for (uv = 0; uv < 2; uv++) {

        curUV = dec_picture->imgUV[uv];

        if (currMB->is_lossless) {
            if (currMB->intra_chroma_pred_mode == Intra_Chroma_Vertical || currMB->intra_chroma_pred_mode == Intra_Chroma_Horizontal)
                Inv_Residual_trans_Chroma(currMB, uv) ;
            else {
                for (j = 0; j < sps->MbHeightC; j++)
                    for (i = 0; i < sps->MbWidthC; i++)
                        currSlice->mb_rres[uv + 1][j][i] = currSlice->cof[uv + 1][j][i];
            }
        }

        if (currMB->mb_type != SI4MB && currMB->cbp >> 4) {
            for (b8 = 0; b8 < num_uv_blocks; b8++) {
                for (b4 = 0; b4 < 4; b4++) {
                    joff = subblk_offset_y[yuv][b8][b4];          
                    ioff = subblk_offset_x[yuv][b8][b4];          

                    itrans_4x4(currMB, (ColorPlane)(uv + 1), ioff, joff);

                    copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_rec[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
                }
            }
            currSlice->is_reset_coeff_cr = FALSE;
        } else if (currMB->mb_type == SI4MB) {
            itrans_sp_cr(currMB, uv);

            for (joff  = 0; joff < 8; joff += 4) {
                for (ioff = 0; ioff < 8; ioff += 4) {
                    itrans_4x4(currMB, (ColorPlane)(uv + 1), ioff, joff);

                    copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_rec[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
                }
            }
            currSlice->is_reset_coeff_cr = FALSE;
        } else { 
            for (b8 = 0; b8 < num_uv_blocks; b8++) {
                for (b4 = 0; b4 < 4; b4++) {
                    joff = subblk_offset_y[yuv][b8][b4];
                    ioff = subblk_offset_x[yuv][b8][b4];          

                    copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_pred[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
                }
            }
        }
    }
}

static int mb_pred_intra(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    int yuv = dec_picture->chroma_format_idc - 1;
    int block4x4;
    int blockoffset = currMB->mb_type == I4MB ? 1
                    : currMB->mb_type == I8MB ? 4 : 16;

    for (block4x4 = 0; block4x4 < 16; block4x4 += blockoffset) {
        int ioff = ((block4x4 / 4) % 2) * 8 + ((block4x4 % 4) % 2) * 4;
        int joff = ((block4x4 / 4) / 2) * 8 + ((block4x4 % 4) / 2) * 4;

        if (currMB->mb_type == I4MB) {
            intra_pred_4x4(currMB, curr_plane, ioff, joff);
            if (currMB->is_lossless)
                Inv_Residual_trans_4x4(currMB, curr_plane, ioff, joff);
            else
                itrans4x4(currMB, curr_plane, ioff, joff);
            copy_image_data_4x4(&currImg[currMB->pix_y + joff], &currSlice->mb_rec[curr_plane][joff], currMB->pix_x + ioff, ioff);
        } else if (currMB->mb_type == I8MB) {
            intra_pred_8x8(currMB, curr_plane, ioff, joff);
            if (currMB->is_lossless)
                Inv_Residual_trans_8x8(currMB, curr_plane, ioff, joff);
            else
                itrans8x8(currMB, curr_plane, ioff, joff);
            copy_image_data_8x8(&currImg[currMB->pix_y + joff], &currSlice->mb_rec[curr_plane][joff], currMB->pix_x + ioff, ioff);
        } else if (currMB->mb_type == I16MB) {
            intra_pred_16x16(currMB, curr_plane, 0, 0);
            if (currMB->is_lossless) {
                Inv_Residual_trans_16x16(currMB, curr_plane);
                copy_image_data_16x16(&currImg[currMB->pix_y + joff], &currSlice->mb_rec[curr_plane][joff], currMB->pix_x + ioff, ioff);
            } else
                iMBtrans4x4(currMB, curr_plane, 0);
        }
    }

    if (dec_picture->chroma_format_idc != YUV400 && dec_picture->chroma_format_idc != YUV444)
        intra_cr_decoding(currMB, yuv);

    if (currMB->mb_type == I16MB || currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}


static int decode_one_component(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;

    currMB->ipmode_DPCM = NO_INTRA_PMODE; 

    if (currMB->mb_type == IPCM)
        mb_pred_ipcm(currMB);
    else if (currMB->mb_type == I16MB || currMB->mb_type == I4MB || currMB->mb_type == I8MB)
        mb_pred_intra(currMB, curr_plane, currImg, dec_picture);
    else if (currMB->mb_type == P16x16 || currMB->mb_type == P16x8 || currMB->mb_type == P8x16)
        mb_pred_p_inter(currMB, curr_plane, dec_picture);
    else if (currMB->mb_type == PSKIP) {
        if (currSlice->slice_type == P_SLICE)
            mb_pred_skip(currMB, curr_plane, currImg, dec_picture);
        else if (currSlice->slice_type == SP_SLICE)
            mb_pred_sp_skip(currMB, curr_plane, dec_picture);
        else if (currSlice->slice_type == B_SLICE) {
            if (currSlice->direct_spatial_mv_pred_flag)
                mb_pred_b_d8x8spatial(currMB, curr_plane, currImg, dec_picture);
            else
                mb_pred_b_d8x8temporal(currMB, curr_plane, currImg, dec_picture);
        }
    } else if (currSlice->slice_type == P_SLICE || currSlice->slice_type == SP_SLICE)
        mb_pred_p_inter(currMB, curr_plane, dec_picture);
    else
        mb_pred_b_inter8x8(currMB, curr_plane, dec_picture);

    return 1;
}

// probably a better way (or place) to do this, but I'm not sure what (where) it is [CJV]
// this is intended to make get_block_luma faster, but I'm still performing
// this at the MB level, and it really should be done at the slice level
static void init_cur_imgy(VideoParameters *p_Vid, Slice *currSlice, int pl)
{
    int i, j;
    if (p_Vid->active_sps->separate_colour_plane_flag == 0) {
        StorablePicture *vidref = p_Vid->no_reference_picture;
        int noref = currSlice->framepoc < p_Vid->recovery_poc;
        if (pl == PLANE_Y) {
            for (j = 0; j < 6; j++) {
                for (i = 0; i < currSlice->listXsize[j] ; i++) {
                    StorablePicture *curr_ref = currSlice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = curr_ref->imgY;
                    }
                }
            }
        } else {
            for (j = 0; j < 6; j++) {
                for (i = 0; i < currSlice->listXsize[j]; i++) {
                    StorablePicture *curr_ref = currSlice->listX[j][i];
                    if (curr_ref) {
                        curr_ref->no_ref = noref && (curr_ref == vidref);
                        curr_ref->cur_imgY = curr_ref->imgUV[pl-1]; 
                    }
                }
            }
        }
    }
}

void decode_one_macroblock(Macroblock *currMB)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;  
    StorablePicture *dec_picture = currSlice->dec_picture;

    // macroblock decoding **************************************************
    if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
        if (!currMB->is_intra_block) {
            init_cur_imgy(p_Vid, currSlice, PLANE_Y);
            decode_one_component(currMB, PLANE_Y, dec_picture->imgY, dec_picture);
            init_cur_imgy(p_Vid, currSlice, PLANE_U);
            decode_one_component(currMB, PLANE_U, dec_picture->imgUV[0], dec_picture);
            init_cur_imgy(p_Vid, currSlice, PLANE_V);
            decode_one_component(currMB, PLANE_V, dec_picture->imgUV[1], dec_picture);
        } else {
            decode_one_component(currMB, PLANE_Y, dec_picture->imgY, dec_picture);
            decode_one_component(currMB, PLANE_U, dec_picture->imgUV[0], dec_picture);
            decode_one_component(currMB, PLANE_V, dec_picture->imgUV[1], dec_picture);      
        }
        currSlice->is_reset_coeff = FALSE;
        currSlice->is_reset_coeff_cr = FALSE;
    } else
        decode_one_component(currMB, PLANE_Y, dec_picture->imgY, dec_picture);
}
