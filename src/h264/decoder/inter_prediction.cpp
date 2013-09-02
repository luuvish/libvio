#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "transform.h"
#include "inter_prediction.h"

// These variables relate to the subpel accuracy supported by the software (1/4)
#define BLOCK_SIZE_8x8_SP  32  // BLOCK_SIZE8x8 << 2

static inline int iClip1(int high, int x)
{
  x = imax(x, 0);
  x = imin(x, high);

  return x;
}

static inline int RSHIFT_RND(int x, int a)
{
    return (a > 0) ? ((x + (1 << (a-1) )) >> a) : (x << (-a));
}
static inline int RSHIFT_RND_SF(int x, int a)
{
    return ((x + (1 << (a-1) )) >> a);
}

static void mc_prediction(imgpel *mb_pred,
                          imgpel *block, int block_size_y, int block_size_x,
                          mb_t *currMB, ColorPlane pl, short l0_refframe, int pred_dir)
{
    int weight, offset, denom, color_clip;
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    pps_t *pps = currSlice->active_pps;
    if (currSlice->weighted_pred_flag) {
        short ref_idx_wp = l0_refframe;
        int type = currSlice->slice_type;
        if (currMB->mb_field_decoding_flag &&
            ((pps->weighted_pred_flag && (type == P_SLICE || type == SP_SLICE))||
             (pps->weighted_bipred_idc == 1 && (type == B_SLICE))))
            ref_idx_wp >>= 1;
        weight = currSlice->wp_weight[pred_dir][ref_idx_wp][pl];
        offset = currSlice->wp_offset[pred_dir][ref_idx_wp][pl];
        denom  = pl > 0 ? currSlice->chroma_log2_weight_denom : currSlice->luma_log2_weight_denom;
        color_clip = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;
    }

    for (int j = 0; j < block_size_y; j++) {
        for (int i = 0; i < block_size_x; i++) {
            if (currSlice->weighted_pred_flag) {
                int result = RSHIFT_RND((weight * block[i]), denom) + offset;
                mb_pred[i] = (imgpel)iClip3(0, color_clip, result);
            } else
                mb_pred[i] = block[i];
        }
        mb_pred += MB_BLOCK_SIZE;
        block   += MB_BLOCK_SIZE;
    }
}

static void bi_prediction(imgpel *mb_pred, 
                          imgpel *block_l0, imgpel *block_l1, int block_size_y, int block_size_x,
                          mb_t *currMB, ColorPlane pl, short l0_refframe, short l1_refframe)
{
    int weight0, weight1, offset, denom, color_clip;
    VideoParameters *p_Vid = currMB->p_Vid;
    sps_t *sps = p_Vid->active_sps;
    int weighted_bipred_idc = p_Vid->active_pps->weighted_bipred_idc;
    if (weighted_bipred_idc) {
        slice_t *currSlice = currMB->p_Slice;

        int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                          currMB->mbAddrX % 2 ? 4 : 2 : 0;
        int l0_ref_idx  = (currMB->mb_field_decoding_flag && weighted_bipred_idc == 1) ? l0_refframe >> 1: l0_refframe;
        int l1_ref_idx  = (currMB->mb_field_decoding_flag && weighted_bipred_idc == 1) ? l1_refframe >> 1: l1_refframe;
        int wt_list_offset = (weighted_bipred_idc == 2) ? list_offset : 0;
        int *wp_weight0 = currSlice->wbp_weight[LIST_0 + wt_list_offset][l0_ref_idx][l1_ref_idx];
        int *wp_weight1 = currSlice->wbp_weight[LIST_1 + wt_list_offset][l0_ref_idx][l1_ref_idx];
        int *wp_offset0 = currSlice->wp_offset[LIST_0 + wt_list_offset][l0_ref_idx];
        int *wp_offset1 = currSlice->wp_offset[LIST_1 + wt_list_offset][l1_ref_idx];

        weight0 = wp_weight0[pl];
        weight1 = wp_weight1[pl];
        offset  = (wp_offset0[pl] + wp_offset1[pl] + 1) >> 1;
        denom   = pl > 0 ? currSlice->chroma_log2_weight_denom + 1 : currSlice->luma_log2_weight_denom + 1;
        color_clip = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;
    }

    int row_inc = MB_BLOCK_SIZE - block_size_x;
    for (int j = 0; j < block_size_y; j++) {
        for (int i = 0; i < block_size_x; i++) {
            if (weighted_bipred_idc) {
                int result = RSHIFT_RND((weight0 * *(block_l0++) + weight1 * *(block_l1++)), denom);
                *(mb_pred++) = (imgpel) iClip1(color_clip, result + offset);
            } else
                *(mb_pred++) = (imgpel)(((*(block_l0++) + *(block_l1++)) + 1) >> 1);
        }
        mb_pred += row_inc;
        block_l0 += row_inc;
        block_l1 += row_inc;
    }
}

/*!
 ************************************************************************
 * \brief
 *    Integer positions
 ************************************************************************
 */ 
static void get_block_00(imgpel *block, imgpel *cur_img, int span, int block_size_y)
{
  // fastest to just move an entire block, since block is a temp block is a 256 byte block (16x16)
  // writes 2 lines of 16 imgpel 1 to 8 times depending in block_size_y
  int j;
  
  for (j = 0; j < block_size_y; j += 2)
  { 
    memcpy(block, cur_img, MB_BLOCK_SIZE * sizeof(imgpel));
    block += MB_BLOCK_SIZE;
    cur_img += span;
    memcpy(block, cur_img, MB_BLOCK_SIZE * sizeof(imgpel));
    block += MB_BLOCK_SIZE;
    cur_img += span;
  }
}

/*!
 ************************************************************************
 * \brief
 *    Half vertical
 ************************************************************************
 */ 
static void get_luma_02(imgpel **block, imgpel **cur_imgY, int block_size_y, int block_size_x, int x_pos, int shift_x, int max_imgpel_value, int dx, int dy)
{
    imgpel *p0, *p1, *p2, *p3, *p4, *p5;
    imgpel *orig_line, *cur_line;
    int i, j;
    int result;
    int jj = dy == 3 ? 1 : 0;

    if (dy == 0)
        shift_x = 1;
    if (dx == 0 && dy != 0)
        p0 = &(cur_imgY[ - 2][x_pos]);
    for (j = 0; j < block_size_y; j++) {
        if (dx == 0 && dy != 0)
            cur_line = &(cur_imgY[jj++][x_pos]);
        if (dy == 0)
            cur_line = &cur_imgY[j][x_pos + (dx == 3 ? 1 : 0)];
        if (dy == 0)
            p0 = &cur_imgY[j][x_pos - 2];
        p1 = p0 + shift_x;
        p2 = p1 + shift_x;
        p3 = p2 + shift_x;
        p4 = p3 + shift_x;
        p5 = p4 + shift_x;
        orig_line = block[j];

        for (i = 0; i < block_size_x; i++) {
            result  = (*(p0++) + *(p5++)) - 5 * (*(p1++) + *(p4++)) + 20 * (*(p2++) + *(p3++));

            *orig_line = (imgpel) iClip1(max_imgpel_value, ((result + 16)>>5));
            if (dx == 1 || dx == 3 || dy == 1 || dy == 3)
                *orig_line = (imgpel) ((*orig_line + *(cur_line++) + 1 ) >> 1);
            orig_line++;
        }

        if (dx == 0 && dy != 0)
            p0 = p1 - block_size_x;
    }
}

/*!
 ************************************************************************
 * \brief
 *    Qpel horizontal, Hpel vertical (3, 2)
 ************************************************************************
 */ 
static void get_luma_32(imgpel **block, imgpel **cur_imgY, int **tmp_res, int block_size_y, int block_size_x, int x_pos, int shift_x, int max_imgpel_value, int dx, int dy)
{
    int i, j;
    int *tmp_line;
    imgpel *p0, *p1, *p2, *p3, *p4, *p5;        
    int    *x0, *x1, *x2, *x3, *x4, *x5;  
    imgpel *orig_line;  
    int result;      

    int jj = -2;

    if (dx == 2)
        shift_x = 1;
    if (dx != 2 && dy == 2)
        p0 = &(cur_imgY[ -2][x_pos - 2]);
    for (j = 0; j < block_size_y + (dx == 2 ? 5 : 0); j++) {
        if (dx == 2)
            p0 = &(cur_imgY[jj++][x_pos - 2]);
        p1 = p0 + shift_x;
        p2 = p1 + shift_x;
        p3 = p2 + shift_x;
        p4 = p3 + shift_x;
        p5 = p4 + shift_x;
        tmp_line  = tmp_res[j];

        for (i = 0; i < block_size_x + (dy == 2 ? 5 : 0); i++)
            *(tmp_line++) = (*(p0++) + *(p5++)) - 5 * (*(p1++) + *(p4++)) + 20 * (*(p2++) + *(p3++));
        if (dx != 2 && dy == 2)
            p0 = p1 - (block_size_x + 5);
    }

    if (dx == 2)
        jj = dy == 3 ? 3 : 2;
    for (j = 0; j < block_size_y; j++) {
        if (dx == 2)
            tmp_line = tmp_res[jj++];
        if (dx != 2 && dy == 2)
            tmp_line = &tmp_res[j][dx == 3 ? 3 : 2];
        orig_line = block[j];
        x0 = tmp_res[j];
        x1 = dx != 2 && dy == 2 ? x0 + 1 : tmp_res[j + 1];
        x2 = dx != 2 && dy == 2 ? x1 + 1 : tmp_res[j + 2];
        x3 = dx != 2 && dy == 2 ? x2 + 1 : tmp_res[j + 3];
        x4 = dx != 2 && dy == 2 ? x3 + 1 : tmp_res[j + 4];
        x5 = dx != 2 && dy == 2 ? x4 + 1 : tmp_res[j + 5];

        for (i = 0; i < block_size_x; i++) {
            result  = (*(x0++) + *(x5++)) - 5 * (*(x1++) + *(x4++)) + 20 * (*(x2++) + *(x3++));

            *orig_line = (imgpel) iClip1(max_imgpel_value, ((result + 512)>>10));
            if (dx == 1 || dx == 3 || dy == 1 || dy == 3)
                *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((*(tmp_line++) + 16)>>5))+1)>>1);
            orig_line++;
        }
    }
}

/*!
 ************************************************************************
 * \brief
 *    Qpel horizontal, Qpel vertical (3, 1)
 ************************************************************************
 */ 
static void get_luma_31(imgpel **block, imgpel **cur_imgY, int block_size_y, int block_size_x, int x_pos, int shift_x, int max_imgpel_value, int dx, int dy)
{
    /* Diagonal interpolation */
    int i, j;
    imgpel *p0, *p1, *p2, *p3, *p4, *p5;
    imgpel *orig_line;  
    int result;

    int jj = dy == 3 ? 1 : 0;

    for (j = 0; j < block_size_y; j++) {
        p0 = &cur_imgY[jj++][x_pos - 2];
        p1 = p0 + 1;
        p2 = p1 + 1;
        p3 = p2 + 1;
        p4 = p3 + 1;
        p5 = p4 + 1;

        orig_line = block[j];

        for (i = 0; i < block_size_x; i++) {
            result  = (*(p0++) + *(p5++)) - 5 * (*(p1++) + *(p4++)) + 20 * (*(p2++) + *(p3++));

            *(orig_line++) = (imgpel) iClip1(max_imgpel_value, ((result + 16)>>5));
        }
    }

    p0 = &(cur_imgY[-2][x_pos + (dx == 3 ? 1 : 0)]);
    for (j = 0; j < block_size_y; j++) {
        p1 = p0 + shift_x;
        p2 = p1 + shift_x;
        p3 = p2 + shift_x;
        p4 = p3 + shift_x;
        p5 = p4 + shift_x;
        orig_line = block[j];

        for (i = 0; i < block_size_x; i++) {
            result  = (*(p0++) + *(p5++)) - 5 * (*(p1++) + *(p4++)) + 20 * (*(p2++) + *(p3++));

            *orig_line = (imgpel) ((*orig_line + iClip1(max_imgpel_value, ((result + 16) >> 5)) + 1) >> 1);
            orig_line++;
        }
        p0 = p1 - block_size_x;
    }
}

/*!
 ************************************************************************
 * \brief
 *    Interpolation of 1/4 subpixel
 ************************************************************************
 */ 
void get_block_luma(StorablePicture *curr_ref, int x_pos, int y_pos, int block_size_x, int block_size_y,
                    imgpel **block, int shift_x, int maxold_x, int maxold_y,
                    ColorPlane pl, mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int **tmp_res = currSlice->tmp_res;
    int max_imgpel_value = (1 << (pl > 0 ? sps->BitDepthC : sps->BitDepthY)) - 1;
    imgpel no_ref_value = (imgpel) (pl ? (1 << (sps->BitDepthC - 1)) : (1 << (sps->BitDepthY - 1)));

    if (curr_ref->no_ref) {
        memset(block[0],no_ref_value,block_size_y * block_size_x * sizeof(imgpel));
        return;
    }

    imgpel **cur_imgY = (currMB->p_Vid->active_sps->separate_colour_plane_flag
                      && currMB->p_Slice->colour_plane_id > PLANE_Y)
                      ? curr_ref->imgUV[currMB->p_Slice->colour_plane_id-1]
                      : curr_ref->cur_imgY;
    int dx = (x_pos & 3);
    int dy = (y_pos & 3);
    x_pos >>= 2;
    y_pos >>= 2;
    x_pos = iClip3(-18, maxold_x+2, x_pos);
    y_pos = iClip3(-10, maxold_y+2, y_pos);

    if (dx == 0 && dy == 0)
        get_block_00(&block[0][0], &cur_imgY[y_pos][x_pos], curr_ref->iLumaStride, block_size_y);
    else { /* other positions */
        if (dy == 0 || dx == 0) { /* No vertical interpolation */
            get_luma_02(block, &cur_imgY[ y_pos], block_size_y, block_size_x, x_pos, shift_x, max_imgpel_value, dx, dy);
        } else if (dx == 2 || dy == 2) { /* Vertical & horizontal interpolation */
            get_luma_32(block, &cur_imgY[ y_pos], tmp_res, block_size_y, block_size_x, x_pos, shift_x, max_imgpel_value, dx, dy);
        } else {
            get_luma_31(block, &cur_imgY[ y_pos], block_size_y, block_size_x, x_pos, shift_x, max_imgpel_value, dx, dy);
        }
    }
}

/*!
 ************************************************************************
 * \brief
 *    Chroma (X,X)
 ************************************************************************
 */ 
static void get_chroma_XY(imgpel *block, imgpel *cur_img, int span, int block_size_y, int block_size_x, int w00, int w01, int w10, int w11, int total_scale)
{
    imgpel *cur_row = cur_img;
    imgpel *nxt_row = cur_img + span;

    imgpel *cur_line, *cur_line_p1;
    imgpel *blk_line;
    int result;
    int i, j;
    for (j = 0; j < block_size_y; j++) {
        cur_line    = cur_row;
        cur_line_p1 = nxt_row;
        blk_line = block;
        block += 16;
        cur_row = nxt_row;
        nxt_row += span;
        for (i = 0; i < block_size_x; i++) {
            result  = (w00 * *(cur_line++) + w01 * *(cur_line_p1++));
            result += (w10 * *(cur_line  ) + w11 * *(cur_line_p1  ));
            *(blk_line++) = (imgpel) RSHIFT_RND_SF(result, total_scale);
        }
    }
}

static void get_block_chroma(StorablePicture *curr_ref, int x_pos, int y_pos,
                             int maxold_x, int maxold_y, int block_size_x, int vert_block_size,
                             imgpel *block1, imgpel *block2, VideoParameters *p_Vid)
{
    sps_t *sps = p_Vid->active_sps;
    imgpel no_ref_value = (imgpel)(1 << (sps->BitDepthC - 1));

    int shiftpel_x = sps->chroma_format_idc == YUV400 ? 0 :
                     sps->chroma_format_idc == YUV444 ? 2 : 3;
    int shiftpel_y = sps->chroma_format_idc == YUV400 ? 0 :
                     sps->chroma_format_idc == YUV420 ? 3 : 2;
    int total_scale = shiftpel_x + shiftpel_y;

    int subpel_x = sps->chroma_format_idc == YUV400 ? 0 :
                   sps->chroma_format_idc == YUV444 ? 3 : 7;
    int subpel_y = sps->chroma_format_idc == YUV400 ? 0 :
                   sps->chroma_format_idc == YUV420 ? 7 : 3;

    int iChromaPadX = MCBUF_CHROMA_PAD_X;
    int iChromaPadY = MCBUF_CHROMA_PAD_Y;
    if (sps->chroma_format_idc == YUV422)
        iChromaPadY = MCBUF_CHROMA_PAD_Y * 2;
    else if (sps->chroma_format_idc == YUV444) {
        iChromaPadX = MCBUF_LUMA_PAD_X;
        iChromaPadY = MCBUF_LUMA_PAD_Y;
    }

    imgpel *img1,*img2;
    short dx,dy;
    int span = curr_ref->iChromaStride;
    if (curr_ref->no_ref) {
        memset(block1,no_ref_value,vert_block_size * block_size_x * sizeof(imgpel));
        memset(block2,no_ref_value,vert_block_size * block_size_x * sizeof(imgpel));
    } else {
        dx = (short) (x_pos & subpel_x);
        dy = (short) (y_pos & subpel_y);
        x_pos = x_pos >> shiftpel_x;
        y_pos = y_pos >> shiftpel_y;
        //clip MV;
        assert(vert_block_size <= iChromaPadY && block_size_x <= iChromaPadX);
        x_pos = iClip3(-iChromaPadX, maxold_x, x_pos); //16
        y_pos = iClip3(-iChromaPadY, maxold_y, y_pos); //8
        img1 = &curr_ref->imgUV[0][y_pos][x_pos];
        img2 = &curr_ref->imgUV[1][y_pos][x_pos];

        if (dx == 0 && dy == 0) {
            get_block_00(block1, img1, span, vert_block_size);
            get_block_00(block2, img2, span, vert_block_size);
        } else {
            short dxcur = (short) (subpel_x + 1 - dx);
            short dycur = (short) (subpel_y + 1 - dy);
            short w00 = dxcur * dycur;
            short w01 = dxcur * dy;
            short w10 = dx * dycur;
            short w11 = dx * dy;
            get_chroma_XY(block1, img1, span, vert_block_size, block_size_x, w00, w01, w10, w11, total_scale);
            get_chroma_XY(block2, img2, span, vert_block_size, block_size_x, w00, w01, w10, w11, total_scale);
        }
    }
}


static void check_motion_vector_range(const MotionVector *mv, slice_t *pSlice)
{  
    if (mv->mv_x > 8191 || mv->mv_x < -8192)
        fprintf(stderr,"WARNING! Horizontal motion vector %d is out of allowed range {-8192, 8191} in picture %d, macroblock %d\n", mv->mv_x, pSlice->p_Vid->number, pSlice->current_mb_nr);
    if (mv->mv_y > (pSlice->max_mb_vmv_r - 1) || mv->mv_y < (-pSlice->max_mb_vmv_r))
        fprintf(stderr,"WARNING! Vertical motion vector %d is out of allowed range {%d, %d} in picture %d, macroblock %d\n", mv->mv_y, (-pSlice->max_mb_vmv_r), (pSlice->max_mb_vmv_r - 1), pSlice->p_Vid->number, pSlice->current_mb_nr);
}

static int CheckVertMV(mb_t *currMB, int vec_y, int block_size_y)
{
    StorablePicture *dec_picture = currMB->p_Slice->dec_picture;
    int y_pos = vec_y >> 2;
    int maxold_y = (currMB->mb_field_decoding_flag) ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y_m1;

    if (block_size_y <= (MCBUF_LUMA_PAD_Y-4))
        return 0;

    if (y_pos < (2 - MCBUF_LUMA_PAD_Y) ||
        y_pos > (maxold_y + MCBUF_LUMA_PAD_Y - block_size_y - 2))
        return 1;
    else
        return 0;
}

void perform_mc(mb_t *currMB, ColorPlane pl, StorablePicture *dec_picture, int pred_dir, int i, int j, int block_size_x, int block_size_y)
{
    assert (pred_dir <= 2);

    static const int mv_mul = 16;
    int vec1_x, vec1_y, vec2_x, vec2_y;
    VideoParameters *p_Vid = currMB->p_Vid;    
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int i4 = currMB->block_x + i;
    int j4 = currMB->block_y + j;
    int ioff = (i << 2);
    int joff = (j << 2);
    int chroma_format_idc = dec_picture->chroma_format_idc;
    PicMotionParams *mv_info = &dec_picture->mv_info[j4][i4];
    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;

    // vars for get_block_luma
    int shift_x  = dec_picture->iLumaStride;
    int maxold_x = dec_picture->size_x_m1;
    int maxold_y = (currMB->mb_field_decoding_flag) ? (dec_picture->size_y >> 1) - 1 : dec_picture->size_y_m1;   
    imgpel **tmp_block_l0 = currSlice->tmp_block_l0;
    imgpel **tmp_block_l1 = currSlice->tmp_block_l1;
    imgpel **tmp_block_l2 = currSlice->tmp_block_l2;
    imgpel **tmp_block_l3 = currSlice->tmp_block_l3;

    MotionVector *l0_mv_array, *l1_mv_array;
    short l0_refframe, l1_refframe;
    StorablePicture *list0, *list1;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;

    if (pred_dir != 2) {
        l0_mv_array = &mv_info->mv[pred_dir];
        l0_refframe = mv_info->ref_idx[pred_dir];
        list0 = currSlice->listX[pred_dir + list_offset][l0_refframe];
    } else {
        l0_mv_array = &mv_info->mv[LIST_0];
        l1_mv_array = &mv_info->mv[LIST_1];
        l0_refframe = mv_info->ref_idx[LIST_0];
        l1_refframe = mv_info->ref_idx[LIST_1];
        list0 = currSlice->listX[LIST_0 + list_offset][l0_refframe];
        list1 = currSlice->listX[LIST_1 + list_offset][l1_refframe];
    }

    int block_y_aff;
    if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
        block_y_aff = (currMB->block_y - 4 * (currMB->mbAddrX % 2)) / 2;
    else
        block_y_aff = currMB->block_y;

    check_motion_vector_range(l0_mv_array, currSlice);
    vec1_x = i4 * mv_mul + l0_mv_array->mv_x;
    vec1_y = (block_y_aff + j) * mv_mul + l0_mv_array->mv_y;
    if (pred_dir == 2) {
        check_motion_vector_range(l1_mv_array, currSlice);
        vec2_x = i4 * mv_mul + l1_mv_array->mv_x;
        vec2_y = (block_y_aff + j) * mv_mul + l1_mv_array->mv_y;
    }

    if (CheckVertMV(currMB, vec1_y, block_size_y)) {
        get_block_luma(list0, vec1_x, vec1_y,
                       block_size_x, BLOCK_SIZE_8x8, tmp_block_l0, shift_x,
                       maxold_x, maxold_y, pl, currMB);
        get_block_luma(list0, vec1_x, vec1_y+BLOCK_SIZE_8x8_SP,
                       block_size_x, block_size_y-BLOCK_SIZE_8x8, tmp_block_l0+BLOCK_SIZE_8x8, shift_x,
                       maxold_x, maxold_y, pl, currMB);
    } else
        get_block_luma(list0, vec1_x, vec1_y,
                       block_size_x, block_size_y, tmp_block_l0, shift_x,
                       maxold_x, maxold_y, pl, currMB);
    if (pred_dir == 2) {
        if (CheckVertMV(currMB, vec2_y, block_size_y)) {
            get_block_luma(list1, vec2_x, vec2_y,
                           block_size_x, BLOCK_SIZE_8x8, tmp_block_l1, shift_x,
                           maxold_x, maxold_y, pl, currMB);
            get_block_luma(list1, vec2_x, vec2_y+BLOCK_SIZE_8x8_SP,
                           block_size_x, block_size_y-BLOCK_SIZE_8x8, tmp_block_l1 + BLOCK_SIZE_8x8, shift_x,
                           maxold_x, maxold_y, pl, currMB);
        } else
            get_block_luma(list1, vec2_x, vec2_y,
                           block_size_x, block_size_y, tmp_block_l1, shift_x,
                           maxold_x, maxold_y, pl, currMB);
    }

    if (pred_dir != 2)
        mc_prediction(&currSlice->mb_pred[pl][joff][ioff],
                      tmp_block_l0[0], block_size_y, block_size_x,
                      currMB, pl, l0_refframe, pred_dir); 
    else
        bi_prediction(&currSlice->mb_pred[pl][joff][ioff],
                      tmp_block_l0[0], tmp_block_l1[0], block_size_y, block_size_x,
                      currMB, pl, l0_refframe, l1_refframe);

    if (chroma_format_idc != YUV400 && chroma_format_idc != YUV444) {
        int ioff_cr, joff_cr, block_size_y_cr, block_size_x_cr, vec2_y_cr, vec1_y_cr;
        int maxold_x = dec_picture->size_x_cr_m1;
        int maxold_y = currMB->mb_field_decoding_flag ? (dec_picture->size_y_cr >> 1) - 1 : dec_picture->size_y_cr_m1;

        if (mb_cr_size_x == MB_BLOCK_SIZE) {
            ioff_cr = ioff;
            block_size_x_cr = block_size_x;
        } else {
            ioff_cr = ioff >> 1;
            block_size_x_cr = block_size_x >> 1;
        }
        if (mb_cr_size_y == MB_BLOCK_SIZE) {
            joff_cr = joff;
            block_size_y_cr = block_size_y;
        } else {
            joff_cr = joff >> 1;
            block_size_y_cr = block_size_y >> 1;
        }

        if (pred_dir != 2) {
            if (chroma_format_idc == 1)
                vec1_y_cr = vec1_y + currSlice->chroma_vector_adjustment[pred_dir + list_offset][l0_refframe]; 
            else
                vec1_y_cr = vec1_y;
        } else {
            if (chroma_format_idc == 1) {
                vec1_y_cr = vec1_y + currSlice->chroma_vector_adjustment[LIST_0 + list_offset][l0_refframe]; 
                vec2_y_cr = vec2_y + currSlice->chroma_vector_adjustment[LIST_1 + list_offset][l1_refframe]; 
            } else {
                vec1_y_cr = vec1_y;
                vec2_y_cr = vec2_y;
            }
        }

        get_block_chroma(list0, vec1_x, vec1_y_cr,
                         maxold_x, maxold_y, block_size_x_cr, block_size_y_cr,
                         tmp_block_l0[0], tmp_block_l2[0], p_Vid);
        if (pred_dir == 2)
            get_block_chroma(list1, vec2_x, vec2_y_cr,
                             maxold_x, maxold_y, block_size_x_cr, block_size_y_cr,
                             tmp_block_l1[0], tmp_block_l3[0], p_Vid);

        if (pred_dir != 2) {
            mc_prediction(&currSlice->mb_pred[1][joff_cr][ioff_cr],
                          tmp_block_l0[0], block_size_y_cr, block_size_x_cr,
                          currMB, (ColorPlane)1, l0_refframe, pred_dir);
            mc_prediction(&currSlice->mb_pred[2][joff_cr][ioff_cr],
                          tmp_block_l2[0], block_size_y_cr, block_size_x_cr,
                          currMB, (ColorPlane)2, l0_refframe, pred_dir);
        } else {
            bi_prediction(&currSlice->mb_pred[1][joff_cr][ioff_cr],
                          tmp_block_l0[0], tmp_block_l1[0], block_size_y_cr, block_size_x_cr,
                          currMB, (ColorPlane)1, l0_refframe, l1_refframe);
            bi_prediction(&currSlice->mb_pred[2][joff_cr][ioff_cr],
                          tmp_block_l2[0], tmp_block_l3[0], block_size_y_cr, block_size_x_cr,
                          currMB, (ColorPlane)2, l0_refframe, l1_refframe);
        }
    }
}
