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

#include "global.h"
#include "mbuffer.h"
#include "slice.h"
#include "macroblock.h"
#include "transform.h"
#include "intra_prediction.h"
#include "intra_prediction_common.h"


static void intra_cr_decoding(Macroblock *currMB, int yuv)
{
    VideoParameters *p_Vid = currMB->p_Vid;
    Slice *currSlice = currMB->p_Slice;
    StorablePicture *dec_picture = currSlice->dec_picture;
    imgpel **curUV;
    int uv;
    int b8,b4;
    int ioff, joff;
    int i,j;

    currSlice->intra_pred_chroma(currMB);// last argument is ignored, computes needed data for both uv channels

    void (*itrans_4x4)(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
    itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : itrans4x4_ls;

    for (uv = 0; uv < 2; uv++) {

        curUV = dec_picture->imgUV[uv];

        if (currMB->is_lossless) {
            if (currMB->c_ipred_mode == VERT_PRED_8 || currMB->c_ipred_mode == HOR_PRED_8)
                Inv_Residual_trans_Chroma(currMB, uv) ;
            else {
                for (j = 0; j < p_Vid->mb_cr_size_y; j++)
                    for (i = 0; i < p_Vid->mb_cr_size_x; i++)
                        currSlice->mb_rres[uv + 1][j][i] = currSlice->cof[uv + 1][j][i];
            }
        }

        if (currMB->mb_type != SI4MB && currMB->cbp >> 4) {
            for (b8 = 0; b8 < p_Vid->num_uv_blocks; b8++) {
                for (b4 = 0; b4 < 4; b4++) {
                    joff = subblk_offset_y[yuv][b8][b4];          
                    ioff = subblk_offset_x[yuv][b8][b4];          

                    itrans_4x4(currMB, (ColorPlane) (uv + 1), ioff, joff);

                    copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_rec[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
                }
            }
            currSlice->is_reset_coeff_cr = FALSE;
        } else if (currMB->mb_type == SI4MB) {
            itrans_sp_cr(currMB, uv);

            for (joff  = 0; joff < 8; joff += 4) {
                for (ioff = 0; ioff < 8; ioff += 4) {
                    itrans_4x4(currMB, (ColorPlane) (uv + 1), ioff, joff);

                    copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_rec[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
                }
            }
            currSlice->is_reset_coeff_cr = FALSE;
        } else { 
            for (b8 = 0; b8 < p_Vid->num_uv_blocks; b8++) {
                for (b4 = 0; b4 < 4; b4++) {
                    joff = subblk_offset_y[yuv][b8][b4];
                    ioff = subblk_offset_x[yuv][b8][b4];          

                    copy_image_data_4x4(&curUV[currMB->pix_c_y + joff], &(currSlice->mb_pred[uv + 1][joff]), currMB->pix_c_x + ioff, ioff);
                }
            }
        }
    }
}


void set_intra_prediction_modes(Slice *currSlice)
{ 
    if (currSlice->mb_aff_frame_flag) {
        currSlice->intra_pred_4x4    = intra_pred_4x4_mbaff;
        currSlice->intra_pred_8x8    = intra_pred_8x8_mbaff;
        currSlice->intra_pred_16x16  = intra_pred_16x16_mbaff;    
        currSlice->intra_pred_chroma = intra_pred_chroma_mbaff;
    } else {
        currSlice->intra_pred_4x4    = intra_pred_4x4_normal;  
        currSlice->intra_pred_8x8    = intra_pred_8x8_normal;
        currSlice->intra_pred_16x16  = intra_pred_16x16_normal;
        currSlice->intra_pred_chroma = intra_pred_chroma;   
    }
}

int mb_pred_intra4x4(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    int yuv = dec_picture->chroma_format_idc - 1;
    int i=0, j=0,k, j4=0,i4=0;  
    int j_pos, i_pos;
    int ioff,joff;
    int block8x8;   // needed for ABT

    void (*itrans_4x4)(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
    itrans_4x4 = (currMB->is_lossless == FALSE) ? itrans4x4 : Inv_Residual_trans_4x4;    

    for (block8x8 = 0; block8x8 < 4; block8x8++) {
        for (k = block8x8 * 4; k < block8x8 * 4 + 4; k ++) {
            i =  (decode_block_scan[k] & 3);
            j = ((decode_block_scan[k] >> 2) & 3);

            ioff = (i << 2);
            joff = (j << 2);
            i4   = currMB->block_x + i;
            j4   = currMB->block_y + j;
            j_pos = j4 * BLOCK_SIZE;
            i_pos = i4 * BLOCK_SIZE;

            // PREDICTION
            //===== INTRA PREDICTION =====
            if (currSlice->intra_pred_4x4(currMB, curr_plane, ioff, joff, i4, j4) == SEARCH_SYNC)  /* make 4x4 prediction block mpr from given prediction p_Vid->mb_mode */
                return 1;
            // =============== 4x4 itrans ================
            // -------------------------------------------
            itrans_4x4(currMB, curr_plane, ioff, joff);

            copy_image_data_4x4(&currImg[j_pos], &currSlice->mb_rec[curr_plane][joff], i_pos, ioff);
        }
    }

    // chroma decoding *******************************************************
    if (dec_picture->chroma_format_idc != YUV400 && dec_picture->chroma_format_idc != YUV444)
        intra_cr_decoding(currMB, yuv);

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}


int mb_pred_intra16x16(Macroblock *currMB, ColorPlane curr_plane, StorablePicture *dec_picture)
{
    int yuv = dec_picture->chroma_format_idc - 1;

    currMB->p_Slice->intra_pred_16x16(currMB, curr_plane, currMB->i16mode);
    currMB->ipmode_DPCM = (char) currMB->i16mode; //For residual DPCM
    // =============== 4x4 itrans ================
    // -------------------------------------------
    iMBtrans4x4(currMB, curr_plane, 0);

    // chroma decoding *******************************************************
    if (dec_picture->chroma_format_idc != YUV400 && dec_picture->chroma_format_idc != YUV444) 
        intra_cr_decoding(currMB, yuv);

    currMB->p_Slice->is_reset_coeff = FALSE;
    return 1;
}

int mb_pred_intra8x8(Macroblock *currMB, ColorPlane curr_plane, imgpel **currImg, StorablePicture *dec_picture)
{
    Slice *currSlice = currMB->p_Slice;
    int yuv = dec_picture->chroma_format_idc - 1;
    int block8x8;   // needed for ABT

    void (*itrans_8x8)(struct macroblock_dec *currMB, ColorPlane pl, int ioff, int joff);
    itrans_8x8 = (currMB->is_lossless == FALSE) ? itrans8x8 : Inv_Residual_trans_8x8;

    for (block8x8 = 0; block8x8 < 4; block8x8++) {
        //=========== 8x8 BLOCK TYPE ============
        int ioff = (block8x8 & 0x01) << 3;
        int joff = (block8x8 >> 1  ) << 3;

        //PREDICTION
        currSlice->intra_pred_8x8(currMB, curr_plane, ioff, joff);
        if (currMB->cbp & (1 << block8x8)) 
            itrans_8x8(currMB, curr_plane, ioff, joff);      // use inverse integer transform and make 8x8 block m7 from prediction block mpr
        else
            icopy8x8(currMB, curr_plane, ioff, joff);
    
        copy_image_data_8x8(&currImg[currMB->pix_y + joff], &currSlice->mb_rec[curr_plane][joff], currMB->pix_x + ioff, ioff);
    }
    // chroma decoding *******************************************************
    if (dec_picture->chroma_format_idc != YUV400 && dec_picture->chroma_format_idc != YUV444)
        intra_cr_decoding(currMB, yuv);

    if (currMB->cbp != 0)
        currSlice->is_reset_coeff = FALSE;
    return 1;
}

int mb_pred_ipcm(Macroblock *currMB)
{
    int i, j, k;
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    StorablePicture *dec_picture = currSlice->dec_picture;

    for (i = 0; i < MB_BLOCK_SIZE; ++i) {
        for (j = 0;j < MB_BLOCK_SIZE ; ++j)
            dec_picture->imgY[currMB->pix_y + i][currMB->pix_x + j] = (imgpel) currSlice->cof[0][i][j];
    }

    if (dec_picture->chroma_format_idc != YUV400 && p_Vid->separate_colour_plane_flag == 0) {
        for (k = 0; k < 2; ++k) {
            for (i = 0; i < p_Vid->mb_cr_size_y; ++i) {
                for (j = 0;j < p_Vid->mb_cr_size_x; ++j)
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
    currMB->skip_flag = 0;

    //for deblocking filter CABAC
    currMB->s_cbp[0].blk = 0xFFFF;

    //For CABAC decoding of Dquant
    currSlice->last_dquant = 0;
    currSlice->is_reset_coeff = FALSE;
    currSlice->is_reset_coeff_cr = FALSE;
    return 1;
}
