
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
            int list_offset = currMB->list_offset;

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

    currSlice->max_mb_vmv_r = (currSlice->field_pic_flag || ( currMB->mb_field_decoding_flag )) ? p_Vid->max_vmv_r >> 1 : p_Vid->max_vmv_r;
}


static int mb_pred_skip(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
        copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
        copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
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
    int smb = (currSlice->slice_type == SP_SLICE);
//    int block4x4;
//    int blockoffset = currMB->mb_type == P16x16 ? 16
//                    : currMB->mb_type == P16x8 ? 8
//                    : currMB->mb_type == P8x16 ? 8 : 4;

    set_chroma_vector(currMB);

//    for (block4x4 = 0; block4x4 < 16; block4x4 += blockoffset) {
//        int mv_mode  = currMB->b8mode[block4x4 / 4];
//        int pred_dir = currMB->b8pdir[block4x4 / 4];
//
//        int ioff = ((block4x4 / 4) % 2) * 8 + ((block4x4 % 4) % 2) * 4;
//        int joff = ((block4x4 / 4) / 2) * 8 + ((block4x4 % 4) / 2) * 4;
//
//        int block_size_x = currMB->mb_type == P16x8 ? MB_BLOCK_SIZE : BLOCK_SIZE_8x8;
//        int block_size_y = currMB->mb_type == P16x8 ? MB_BLOCK_SIZE : BLOCK_SIZE_8x8;
//    }

    if (currMB->mb_type == P16x16)
        perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    else if (currMB->mb_type == P16x8) {
        perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, MB_BLOCK_SIZE, BLOCK_SIZE_8x8);
        perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[2], 0, 2, MB_BLOCK_SIZE, BLOCK_SIZE_8x8);
    } else if (currMB->mb_type == P8x16) {
        perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, BLOCK_SIZE_8x8, MB_BLOCK_SIZE);
        perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[1], 2, 0, BLOCK_SIZE_8x8, MB_BLOCK_SIZE);
    } else {
        int block8x8;
        for (block8x8 = 0; block8x8 < 4; block8x8++) {
            int mv_mode  = currMB->b8mode[block8x8];
            int pred_dir = currMB->b8pdir[block8x8];

            int k_start = block8x8 * 4;
            int k_inc = (mv_mode == SMB8x4) ? 2 : 1;
            int k_end = (mv_mode == SMB8x8) ? k_start + 1 :
                        (mv_mode == SMB4x4) ? k_start + 4 :
                                              k_start + k_inc + 1;

            int block_size_x = ( mv_mode == SMB8x4 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
            int block_size_y = ( mv_mode == SMB4x8 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;

            int k;
            for (k = k_start; k < k_end; k += k_inc) {
                int i =  (decode_block_scan[k] & 3);
                int j = ((decode_block_scan[k] >> 2) & 3);
                perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, block_size_x, block_size_y);
            }
        }
    }

    iTransform(currMB, curr_plane, smb);

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

static int mb_pred_b_d8x8temporal(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    int k;
    int block8x8;   // needed for ABT
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int dir8x8 = currSlice->active_sps->direct_8x8_inference_flag;

    set_chroma_vector(currMB);

    int k_inc = dir8x8 ? 1 : BLOCK_MULTIPLE;
    int block_size_x = dir8x8 ? SMB_BLOCK_SIZE : BLOCK_SIZE;
    int block_size_y = dir8x8 ? SMB_BLOCK_SIZE : BLOCK_SIZE;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    for (block8x8 = 0; block8x8 < 4; block8x8++) {
        int pred_dir = currMB->b8pdir[block8x8];

        int k_start = (block8x8 << 2);
        int k_end = k_start + k_inc;

        if (dir8x8)
            pred_dir = get_direct8x8temporal(currMB, dec_picture, block8x8);
        else
            pred_dir = get_direct4x4temporal(currMB, dec_picture, block8x8);

        for (k = k_start; k < k_end; k ++) {
            int i =  (decode_block_scan[k] & 3);
            int j = ((decode_block_scan[k] >> 2) & 3);
            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, block_size_x, block_size_y);
        }
    }

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if (dec_picture->chroma_format_idc != YUV400 && dec_picture->chroma_format_idc != YUV444) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }
    return 1;
}

static int mb_pred_b_d8x8spatial(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    int i4, j4;
    int block8x8;
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    int pred_dir = 0;

    set_chroma_vector(currMB);

    prepare_direct_params(currMB, dec_picture, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    if (l0_rFrame == 0 || l1_rFrame == 0) {
        for (block8x8 = 0; block8x8 < 4; block8x8++) {
            int k_start = (block8x8 << 2);

            int i  =  (decode_block_scan[k_start] & 3);
            int j  = ((decode_block_scan[k_start] >> 2) & 3);
            i4  = currMB->block_x + i;
            j4  = currMB->block_y + j;

            pred_dir = get_direct8x8spatial_eq(currMB, dec_picture, block8x8, &pmvl0, &pmvl1, l0_rFrame, l1_rFrame);

            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, SMB_BLOCK_SIZE, SMB_BLOCK_SIZE);
        }
    } else {

        pred_dir = get_direct8x8spatial_ne(currMB, dec_picture, 0, &pmvl0, &pmvl1, l0_rFrame, l1_rFrame);

        // Now perform Motion Compensation
        perform_mc(currMB, curr_plane, dec_picture, pred_dir, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    }

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if (dec_picture->chroma_format_idc != YUV400 && dec_picture->chroma_format_idc != YUV444) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }

    return 1;
}

static int mb_pred_b_d4x4spatial(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    int k;
    int block8x8;
    Slice *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int pred_dir = 0;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    set_chroma_vector(currMB);

    prepare_direct_params(currMB, dec_picture, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    for (block8x8 = 0; block8x8 < 4; block8x8++) {
        int k_start = (block8x8 << 2);
        int k_end = k_start + BLOCK_MULTIPLE;

        pred_dir = get_direct4x4spatial(currMB, dec_picture, block8x8, &pmvl0, &pmvl1, l0_rFrame, l1_rFrame);

        for (k = k_start; k < k_end; k ++) {
            int i =  (decode_block_scan[k] & 3);
            int j = ((decode_block_scan[k] >> 2) & 3);

            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, BLOCK_SIZE, BLOCK_SIZE);
        }
    }

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, mb_cr_size_x, mb_cr_size_y);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }

    return 1;
}

static int mb_pred_b_inter8x8(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    int block_size_x, block_size_y;
    int k;
    int block8x8;   // needed for ABT
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    set_chroma_vector(currMB);
  
    // prepare direct modes
    if (currSlice->direct_spatial_mv_pred_flag
        && (!(currMB->b8mode[0] && currMB->b8mode[1] && currMB->b8mode[2] && currMB->b8mode[3])))
        prepare_direct_params(currMB, dec_picture, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    for (block8x8=0; block8x8<4; block8x8++) {
        int mv_mode  = currMB->b8mode[block8x8];
        int pred_dir = currMB->b8pdir[block8x8];
        int k_start, k_end, k_inc;

        pred_dir = get_inter8x8(currMB, dec_picture, block8x8);

        if ( mv_mode != 0 ) {
            k_start = (block8x8 << 2);
            k_inc = (mv_mode == SMB8x4) ? 2 : 1;
            k_end = (mv_mode == SMB8x8) ? k_start + 1 : ((mv_mode == SMB4x4) ? k_start + 4 : k_start + k_inc + 1);

            block_size_x = ( mv_mode == SMB8x4 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
            block_size_y = ( mv_mode == SMB4x8 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
        } else {
            k_start = (block8x8 << 2);
            k_end = k_start;
            k_inc = 1;

            if (p_Vid->active_sps->direct_8x8_inference_flag) {
                block_size_x = SMB_BLOCK_SIZE;
                block_size_y = SMB_BLOCK_SIZE;
                k_end ++;
            } else {
                block_size_x = BLOCK_SIZE;
                block_size_y = BLOCK_SIZE;
                k_end += BLOCK_MULTIPLE;
            }
        }

        for (k = k_start; k < k_end; k += k_inc) {
            int i =  (decode_block_scan[k] & 3);
            int j = ((decode_block_scan[k] >> 2) & 3);
            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, block_size_x, block_size_y);
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

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    for (i = 0; i < MB_BLOCK_SIZE; ++i) {
        for (j = 0;j < MB_BLOCK_SIZE ; ++j)
            dec_picture->imgY[currMB->pix_y + i][currMB->pix_x + j] = (imgpel) currSlice->cof[0][i][j];
    }

    if (dec_picture->chroma_format_idc != YUV400 && p_Vid->active_sps->separate_colour_plane_flag == 0) {
        for (k = 0; k < 2; ++k) {
            for (i = 0; i < mb_cr_size_y; ++i) {
                for (j = 0;j < mb_cr_size_x; ++j)
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

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

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
            if (currMB->c_ipred_mode == VERT_PRED_8 || currMB->c_ipred_mode == HOR_PRED_8)
                Inv_Residual_trans_Chroma(currMB, uv) ;
            else {
                for (j = 0; j < mb_cr_size_y; j++)
                    for (i = 0; i < mb_cr_size_x; i++)
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
            if (currSlice->direct_spatial_mv_pred_flag == 0)
                mb_pred_b_d8x8temporal(currMB, curr_plane, currImg, dec_picture);
            else {
                if (currSlice->active_sps->direct_8x8_inference_flag)
                    mb_pred_b_d8x8spatial(currMB, curr_plane, currImg, dec_picture);
                else
                    mb_pred_b_d4x4spatial(currMB, curr_plane, currImg, dec_picture);
            }
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
    VideoParameters *p_Vid = currMB->p_Vid;  
    StorablePicture *dec_picture = currSlice->dec_picture;

    // macroblock decoding **************************************************
    if (currSlice->active_sps->chroma_format_idc == YUV444 && currSlice->active_sps->separate_colour_plane_flag == 0) {
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
