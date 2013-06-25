/*!
 *************************************************************************************
 * \file mb_prediction.c
 *
 * \brief
 *    Macroblock prediction functions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Alexis Michael Tourapis         <alexismt@ieee.org>
 *************************************************************************************
 */

#include "transform.h"
#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "inter_prediction.h"
#include "inter_prediction_mc.h"


static void set_chroma_vector(Macroblock *currMB)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    if (!currSlice->mb_aff_frame_flag) {
        if (currSlice->structure == TOP_FIELD) {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++) {
                    if (p_Vid->structure != currSlice->listX[l][k]->structure)
                        currSlice->chroma_vector_adjustment[l][k] = -2; 
                    else
                        currSlice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        } else if(currSlice->structure == BOTTOM_FIELD) {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++) {
                    if (p_Vid->structure != currSlice->listX[l][k]->structure)
                        currSlice->chroma_vector_adjustment[l][k] = 2; 
                    else
                        currSlice->chroma_vector_adjustment[l][k] = 0; 
                }
            }
        } else {
            int k,l;
            for (l = LIST_0; l <= (LIST_1); l++) {
                for (k = 0; k < currSlice->listXsize[l]; k++)
                    currSlice->chroma_vector_adjustment[l][k] = 0; 
            }
        }
    } else {
        int mb_nr = (currMB->mbAddrX & 0x01);
        int k,l;  

        //////////////////////////
        // find out the correct list offsets
        if (currMB->mb_field) {
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

    currSlice->max_mb_vmv_r = (currSlice->structure != FRAME || ( currMB->mb_field )) ? p_Vid->max_vmv_r >> 1 : p_Vid->max_vmv_r;
}


int mb_pred_skip(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

    copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

    if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
        copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
        copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
    }
    return 1;
}

int mb_pred_sp_skip(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, LIST_0, 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    iTransform(currMB, curr_plane, 1);
    return 1;
}

int mb_pred_p_inter8x8(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    int block8x8;   // needed for ABT
    int i=0, j=0,k;  

    Slice *currSlice = currMB->p_Slice;
    int smb = currSlice->slice_type == SP_SLICE && (currMB->is_intra_block == FALSE);

    set_chroma_vector(currMB);

    for (block8x8=0; block8x8<4; block8x8++) {
        int mv_mode  = currMB->b8mode[block8x8];
        int pred_dir = currMB->b8pdir[block8x8];

        int k_start = (block8x8 << 2);
        int k_inc = (mv_mode == SMB8x4) ? 2 : 1;
        int k_end = (mv_mode == SMB8x8) ? k_start + 1 : ((mv_mode == SMB4x4) ? k_start + 4 : k_start + k_inc + 1);

        int block_size_x = ( mv_mode == SMB8x4 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;
        int block_size_y = ( mv_mode == SMB4x8 || mv_mode == SMB8x8 ) ? SMB_BLOCK_SIZE : BLOCK_SIZE;

        for (k = k_start; k < k_end; k += k_inc) {
            i =  (decode_block_scan[k] & 3);
            j = ((decode_block_scan[k] >> 2) & 3);
            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, block_size_x, block_size_y);
        }
    }

    iTransform(currMB, curr_plane, smb); 

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

int mb_pred_p_inter16x16(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    int smb = (currSlice->slice_type == SP_SLICE);

    set_chroma_vector(currMB);
    perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);
    iTransform(currMB, curr_plane, smb);

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

int mb_pred_p_inter16x8(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    int smb = (currSlice->slice_type == SP_SLICE);

    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, MB_BLOCK_SIZE, BLOCK_SIZE_8x8);
    perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[2], 0, 2, MB_BLOCK_SIZE, BLOCK_SIZE_8x8);
    iTransform(currMB, curr_plane, smb); 
  
    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

int mb_pred_p_inter8x16(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    int smb = (currSlice->slice_type == SP_SLICE);

    set_chroma_vector(currMB);

    perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[0], 0, 0, BLOCK_SIZE_8x8, MB_BLOCK_SIZE);
    perform_mc(currMB, curr_plane, dec_picture, currMB->b8pdir[1], 2, 0, BLOCK_SIZE_8x8, MB_BLOCK_SIZE);
    iTransform(currMB, curr_plane, smb);

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

int mb_pred_b_d8x8temporal(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    int k;
    int block8x8;   // needed for ABT
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    set_chroma_vector(currMB);

    for (block8x8=0; block8x8<4; block8x8++) {
        int pred_dir = currMB->b8pdir[block8x8];

        int k_start = (block8x8 << 2);
        int k_end = k_start + 1;

        pred_dir = get_direct8x8temporal(currMB, dec_picture, block8x8);

        for (k = k_start; k < k_end; k ++) {
            int i =  (decode_block_scan[k] & 3);
            int j = ((decode_block_scan[k] >> 2) & 3);
            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, SMB_BLOCK_SIZE, SMB_BLOCK_SIZE);
        }
    }

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }
    return 1;
}

int mb_pred_b_d4x4temporal(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    int k;
    int block8x8;   // needed for ABT
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    set_chroma_vector(currMB);

    for (block8x8=0; block8x8<4; block8x8++) {
        int pred_dir = currMB->b8pdir[block8x8];

        int k_start = (block8x8 << 2);
        int k_end = k_start + BLOCK_MULTIPLE;

        pred_dir = get_direct4x4temporal(currMB, dec_picture, block8x8);

        for (k = k_start; k < k_end; k ++) {
            int i =  (decode_block_scan[k] & 3);
            int j = ((decode_block_scan[k] >> 2) & 3);
            perform_mc(currMB, curr_plane, dec_picture, pred_dir, i, j, BLOCK_SIZE, BLOCK_SIZE);
        }
    }

    if (currMB->cbp == 0) {
        copy_image_data_16x16(&currImg[currMB->pix_y], currSlice->mb_pred[curr_plane], currMB->pix_x, 0);

        if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }

    return 1;
}

int mb_pred_b_d8x8spatial(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    int i4, j4;
    int block8x8;
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

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

        if ((dec_picture->chroma_format_idc != YUV400) && (dec_picture->chroma_format_idc != YUV444)) {
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }

    return 1;
}

int mb_pred_b_d4x4spatial(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    int k;
    int block8x8;
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    int pred_dir = 0;

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
            copy_image_data(&dec_picture->imgUV[0][currMB->pix_c_y], currSlice->mb_pred[1], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
            copy_image_data(&dec_picture->imgUV[1][currMB->pix_c_y], currSlice->mb_pred[2], currMB->pix_c_x, 0, p_Vid->mb_size[1][0], p_Vid->mb_size[1][1]);
        }
    } else {
        iTransform(currMB, curr_plane, 0); 
        currSlice->is_reset_coeff = FALSE;
    }

    return 1;
}

int mb_pred_b_inter8x8(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
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
