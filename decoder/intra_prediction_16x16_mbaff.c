/*!
 *************************************************************************************
 * \file intra16x16_pred_mbaff.c
 *
 * \brief
 *    Functions for intra 16x16 prediction (MBAFF)
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
 *      - Yuri Vatis
 *      - Jan Muenster
 *      - Alexis Michael Tourapis  <alexismt@ieee.org>
 *
 *************************************************************************************
 */
#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "intra_prediction_common.h"
#include "intra_prediction.h"
#include "neighbour.h"
#include "image.h"


static void get_available(int *avail_left, int *avail_up, int *avail_up_right, int *avail_up_left,
                          PixelPos *pix_a, PixelPos *pix_b, PixelPos *pix_c, PixelPos *pix_d,
                          Macroblock *currMB, int ioff, int joff)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    int i;

    for (i = 0; i < 17; ++i)
        getAffNeighbour(currMB, -1, i - 1, p_Vid->mb_size[IS_LUMA], &pix_a[i]);
    getAffNeighbour(currMB, 0, -1, p_Vid->mb_size[IS_LUMA], pix_b);

    if (!p_Vid->active_pps->constrained_intra_pred_flag) {
        *avail_left    = pix_a[1].available;
        *avail_up      = pix_b->available;
        *avail_up_left = pix_a[0].available;
    } else {
        *avail_left    = 1;
        for (i = 1; i < 17; ++i)
            *avail_left &= pix_a[i].available ? currSlice->intra_block[pix_a[i].mb_addr]: 0;
        *avail_up      = pix_b->available ? currSlice->intra_block[pix_b->mb_addr] : 0;
        *avail_up_left = pix_a[0].available ? currSlice->intra_block[pix_a[0].mb_addr]: 0;
    }
}


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 DC prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *
 ***********************************************************************
 */
static void intra16x16_dc_pred_mbaff(Macroblock *currMB, ColorPlane pl)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;

    int s0 = 0, s1 = 0, s2 = 0;
    int i, j;

    imgpel **imgY = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    imgpel **mb_pred = &(currSlice->mb_pred[pl][0]); 

    PixelPos pix_a[17], pix_b, pix_c, pix_d;
    int left_avail, up_avail, right_up_avail, left_up_avail;
    get_available(&left_avail, &up_avail, &right_up_avail, &left_up_avail,
                  &pix_a[0], &pix_b, &pix_c, &pix_d,
                  currMB, 0, 0);

    for (i = 0; i < MB_BLOCK_SIZE; ++i) {
        if (up_avail)
            s1 += imgY[pix_b.pos_y][pix_b.pos_x+i];    // sum hor pix
        if (left_avail)
            s2 += imgY[pix_a[i + 1].pos_y][pix_a[i + 1].pos_x];    // sum vert pix
    }
    if (up_avail && left_avail)
        s0 = (s1 + s2 + 16)>>5;       // no edge
    else if (!up_avail && left_avail)
        s0 = (s2 + 8)>>4;              // upper edge
    else if (up_avail && !left_avail)
        s0 = (s1 + 8)>>4;              // left edge
    else
        s0 = p_Vid->dc_pred_value_comp[pl];                            // top left corner, nothing to predict from

    for (j = 0; j < MB_BLOCK_SIZE; ++j) {
        for (i = 0; i < MB_BLOCK_SIZE; ++i)
            mb_pred[j][i]=(imgpel) s0;
    }
}



/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 vertical prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *
 ***********************************************************************
 */
static void intra16x16_vert_pred_mbaff(Macroblock *currMB, ColorPlane pl)
{
    Slice *currSlice = currMB->p_Slice;
  
    int j;

    imgpel **imgY = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;

    PixelPos pix_a[17], pix_b, pix_c, pix_d;
    int left_avail, up_avail, right_up_avail, left_up_avail;
    get_available(&left_avail, &up_avail, &right_up_avail, &left_up_avail,
                  &pix_a[0], &pix_b, &pix_c, &pix_d,
                  currMB, 0, 0);

    if (!up_avail) {
        error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);
        return;
    }

    imgpel **prd = &currSlice->mb_pred[pl][0];
    imgpel *src = &(imgY[pix_b.pos_y][pix_b.pos_x]);

    for (j = 0; j < MB_BLOCK_SIZE; j+= 4) {
        memcpy(*prd++, src, MB_BLOCK_SIZE * sizeof(imgpel));
        memcpy(*prd++, src, MB_BLOCK_SIZE * sizeof(imgpel));
        memcpy(*prd++, src, MB_BLOCK_SIZE * sizeof(imgpel));
        memcpy(*prd++, src, MB_BLOCK_SIZE * sizeof(imgpel));
    }
}


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 horizontal prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *
 ***********************************************************************
 */
static void intra16x16_hor_pred_mbaff(Macroblock *currMB, ColorPlane pl)
{
    Slice *currSlice = currMB->p_Slice;
    int i, j;

    imgpel **imgY = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    imgpel **mb_pred = &(currSlice->mb_pred[pl][0]); 
    imgpel prediction;

    PixelPos pix_a[17], pix_b, pix_c, pix_d;
    int left_avail, up_avail, right_up_avail, left_up_avail;
    get_available(&left_avail, &up_avail, &right_up_avail, &left_up_avail,
                  &pix_a[0], &pix_b, &pix_c, &pix_d,
                  currMB, 0, 0);

    if (!left_avail) {
        error ("invalid 16x16 intra pred Mode HOR_PRED_16",500);
        return;
    }

    for (j = 0; j < MB_BLOCK_SIZE; ++j) {
        prediction = imgY[pix_a[j+1].pos_y][pix_a[j+1].pos_x];
        for (i = 0; i < MB_BLOCK_SIZE; ++i)
            mb_pred[j][i]= prediction; // store predicted 16x16 block
    }
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 horizontal prediction mode
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *
 ***********************************************************************
 */
static void intra16x16_plane_pred_mbaff(Macroblock *currMB, ColorPlane pl)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
  
    int i, j;

    int ih = 0, iv = 0;
    int ib, ic, iaa;

    imgpel **imgY = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    imgpel **mb_pred = &(currSlice->mb_pred[pl][0]); 
    imgpel *mpr_line;
    int max_imgpel_value = p_Vid->max_pel_value_comp[pl];

    PixelPos pix_a[17], pix_b, pix_c, pix_d;
    int left_avail, up_avail, right_up_avail, left_up_avail;
    get_available(&left_avail, &up_avail, &right_up_avail, &left_up_avail,
                  &pix_a[0], &pix_b, &pix_c, &pix_d,
                  currMB, 0, 0);

    if (!up_avail || !left_up_avail  || !left_avail) {
        error ("invalid 16x16 intra pred Mode PLANE_16",500);
        return;
    }

    mpr_line = &imgY[pix_b.pos_y][pix_b.pos_x+7];
    for (i = 1; i < 8; ++i) {
        ih += i*(mpr_line[i] - mpr_line[-i]);
        iv += i*(imgY[pix_a[8+i].pos_y][pix_a[8+i].pos_x] - imgY[pix_a[8-i].pos_y][pix_a[8-i].pos_x]);
    }

    ih += 8*(mpr_line[8] - imgY[pix_a[0].pos_y][pix_a[0].pos_x]);
    iv += 8*(imgY[pix_a[16].pos_y][pix_a[16].pos_x] - imgY[pix_a[0].pos_y][pix_a[0].pos_x]);

    ib=(5 * ih + 32)>>6;
    ic=(5 * iv + 32)>>6;

    iaa=16 * (mpr_line[8] + imgY[pix_a[16].pos_y][pix_a[16].pos_x]);
    for (j = 0; j < MB_BLOCK_SIZE; ++j) {
        for (i = 0; i < MB_BLOCK_SIZE; ++i)
      mb_pred[j][i] = (imgpel) iClip1(max_imgpel_value, ((iaa + (i - 7) * ib + (j - 7) * ic + 16) >> 5));
    }
}

/*!
 ***********************************************************************
 * \brief
 *    makes and returns 16x16 intra prediction blocks 
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */
int intra_pred_16x16_mbaff(Macroblock *currMB,  //!< Current Macroblock
                          ColorPlane pl,       //!< Current colorplane (for 4:4:4)                         
                          int predmode)        //!< prediction mode
{
    switch (predmode) {
    case VERT_PRED_16:                       // vertical prediction from block above
        intra16x16_vert_pred_mbaff(currMB, pl);
        break;
    case HOR_PRED_16:                        // horizontal prediction from left block
        intra16x16_hor_pred_mbaff(currMB, pl);
        break;
    case DC_PRED_16:                         // DC prediction
        intra16x16_dc_pred_mbaff(currMB, pl);
        break;
    case PLANE_16:// 16 bit integer plan pred
        intra16x16_plane_pred_mbaff(currMB, pl);
        break;
    default:
        printf("illegal 16x16 intra prediction mode input: %d\n",predmode);
        return SEARCH_SYNC;
    }

    return DECODING_OK;
}
