/*!
 *************************************************************************************
 * \file intra4x4_pred.c
 *
 * \brief
 *    Functions for intra 4x4 prediction
 *
 * \author
 *      Main contributors (see contributors.h for copyright, 
 *                         address and affiliation details)
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


static void neighbouring_samples(imgpel *pred, bool *available, Macroblock *currMB, ColorPlane pl, int xO, int yO)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    int *mb_size = (pl) ? p_Vid->mb_size[IS_CHROMA] : p_Vid->mb_size[IS_LUMA];
    PixelPos pix_a, pix_b, pix_c, pix_d, pix_x[4];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->mb_aff_frame_flag;

    if (MbaffFrameFlag == 0) {
        getNonAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < 4; i++)
            getNonAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getNonAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getNonAffNeighbour(currMB, xO + 4, yO - 1, mb_size, &pix_c);
        getNonAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    } else {
        getAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < 4; i++)
            getAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getAffNeighbour(currMB, xO + 4, yO - 1, mb_size, &pix_c);
        getAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    }
    pix_c.available = pix_c.available && !((xO==4) && ((yO==4)||(yO==12)));

    if (constrained_intra_pred_flag) {
        available[0] = pix_a.available ? currSlice->intra_block[pix_a.mb_addr] : 0;
        for (int i = 0; i < 4; i++)
            available[0] &= pix_x[i].available ? currSlice->intra_block[pix_x[i].mb_addr] : 0;
        available[1] = pix_b.available ? currSlice->intra_block[pix_b.mb_addr] : 0;
        available[2] = pix_c.available ? currSlice->intra_block[pix_c.mb_addr] : 0;
        available[3] = pix_d.available ? currSlice->intra_block[pix_d.mb_addr] : 0;
    } else {
        available[0] = pix_a.available;
        available[1] = pix_b.available;
        available[2] = pix_c.available;
        available[3] = pix_d.available;
    }

#define p(x,y) (pred[((y) + 1) * 9 + ((x) + 1)])
    if (available[3]) {
        imgpel *pix = &img[pix_d.pos_y][pix_d.pos_x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        //int width = (pl) ? currSlice->dec_picture->iChromaStride : currSlice->dec_picture->iLumaStride;
        //int dy = MbaffFrameFlag == 1 && p_Vid->mb_data[pix_a.mb_addr].mb_field ? 2 : 1;
//        imgpel *pix = &img[pix_a.pos_y][pix_a.pos_x];
//        if (MbaffFrameFlag == 1 && !currMB->mb_field && (currMB->mbAddrX & 1))
//            pix -= 16 * width;
//        for (int y = 0; y < 4; y++)
//            p(-1, y) = pix[y * dy * width];
        for (int y = 0; y < 4; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[1]) {
        imgpel *pix = &img[pix_b.pos_y][pix_b.pos_x];
        for (int x = 0; x < 4; x++)
            p(x, -1) = pix[x];
        pix = &img[pix_c.pos_y][pix_c.pos_x-4];
        for (int x = 4; x < 8; x++)
            p(x, -1) = available[2] ? pix[x] : p(3, -1);
        available[2] = available[1];
    }
#undef p
}

#define pred4x4L(x,y) (pred[(y) * 16 + (x)])
#define p(x,y) (pix[((y) + 1) * 9 + ((x) + 1)])

static void intra4x4_vert_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            pred4x4L(x, y) = p(x, -1);
        }
    }
}

static void intra4x4_hor_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            pred4x4L(x, y) = p(-1, y);
        }
    }
}

static void intra4x4_dc_pred_normal(imgpel *pred, imgpel *pix, bool *available, int BitDepth)
{
    bool block_available_left = available[0];
    bool block_available_up   = available[1];

    int shift = 1;
    int sum = 0;
    if (block_available_left) {
        for (int y = 0; y < 4; y++) {
            sum += p(-1, y);
        }
        sum += 2;
        shift++;
    }
    if (block_available_up) {
        for (int x = 0; x < 4; x++) {
            sum += p(x, -1);
        }
        sum += 2;
        shift++;
    }
    sum >>= shift;
    if (!block_available_left && !block_available_up)
        sum = 1 << (BitDepth - 1);

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            pred4x4L(x, y) = sum;
        }
    }
}

static void intra4x4_diag_down_left_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if (x == 3 && y == 3)
                pred4x4L(x, y) = (p(x + y, -1) + 3 * p(x + y + 1, -1) + 2) >> 2;
            else
                pred4x4L(x, y) = (p(x + y, -1) + 2 * p(x + y + 1, -1) + p(x + y + 2, -1) + 2) >> 2;
        }
    }
}

static void intra4x4_diag_down_right_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if (x > y)
                pred4x4L(x, y) = (p(x - y - 2, -1) + 2 * p(x - y - 1, -1) + p(x - y, -1) + 2) >> 2;
            else if (x < y)
                pred4x4L(x, y) = (p(-1, y - x - 2) + 2 * p(-1, y - x - 1) + p(-1, y - x) + 2) >> 2;
            else
                pred4x4L(x, y) = (p(0, -1) + 2 * p(-1, -1) + p(-1, 0) + 2) >> 2;
        }
    }
}

static void intra4x4_vert_right_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            int zVR = 2 * x - y;
            if (zVR == 0 || zVR == 2 || zVR == 4 || zVR == 6)
                pred4x4L(x, y) = (p(x - (y >> 1) - 1, -1) + p(x - (y >> 1), -1) + 1) >> 1;
            else if (zVR == 1 || zVR == 3 || zVR == 5)
                pred4x4L(x, y) = (p(x - (y >> 1) - 2, -1) + 2 * p(x - (y >> 1) - 1, -1) + p(x - (y >> 1), -1) + 2) >> 2;
            else if (zVR == -1)
                pred4x4L(x, y) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            else // zVR == -2 || zVR == -3
                pred4x4L(x, y) = (p(-1, y - 1) + 2 * p(-1, y - 2) + p(-1, y - 3) + 2) >> 2;
        }
    }
}

static void intra4x4_hor_down_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            int zHD = 2 * y - x;
            if (zHD == 0 || zHD == 2 || zHD == 4 || zHD == 6)
                pred4x4L(x, y) = (p(-1, y - (x >> 1) - 1) + p(-1, y - (x >> 1)) + 1) >> 1;
            else if (zHD == 1 || zHD == 3 || zHD == 5)
                pred4x4L(x, y) = (p(-1, y - (x >> 1) - 2) + 2 * p(-1, y - (x >> 1) - 1) + p(-1, y - (x >> 1)) + 2) >> 2;
            else if (zHD == -1)
                pred4x4L(x, y) = (p(-1, 0) + 2 * p(-1, -1) + p(0, -1) + 2) >> 2;
            else // zHD == -2 || zHD == -3
                pred4x4L(x, y) = (p(x - 1, -1) + 2 * p(x - 2, -1) + p(x - 3, -1) + 2) >> 2;
        }
    }
}

static void intra4x4_vert_left_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            if (y == 0 || y == 2)
                pred4x4L(x, y) = (p(x + (y >> 1), -1) + p(x + (y >> 1) + 1, -1) + 1) >> 1;
            else // y == 1 || y == 3
                pred4x4L(x, y) = (p(x + (y >> 1), -1) + 2 * p(x + (y >> 1) + 1, -1) + p(x + (y >> 1) + 2, -1) + 2) >> 2;
        }
    }
}

static void intra4x4_hor_up_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            int zHU = x + 2 * y;
            if (zHU == 0 || zHU == 2 || zHU == 4)
                pred4x4L(x, y) = (p(-1, y + (x >> 1)) + p(-1, y + (x >> 1) + 1) + 1) >> 1;
            else if (zHU == 1 || zHU == 3)
                pred4x4L(x, y) = (p(-1, y + (x >> 1)) + 2 * p(-1, y + (x >> 1) + 1) + p(-1, y + (x >> 1) + 2) + 2) >> 2;
            else if (zHU == 5)
                pred4x4L(x, y) = (p(-1, 2) + 3 * p(-1, 3) + 2) >> 2;
            else // zHU >= 5
                pred4x4L(x, y) = p(-1, 3);
        }
    }
}

#undef p
#undef pred4x4L


/*!
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 intra prediction blocks 
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */
int intra_pred_4x4(Macroblock *currMB,    //!< current macroblock
                  ColorPlane pl,         //!< current image plane
                  int ioff,              //!< pixel offset X within MB
                  int joff,              //!< pixel offset Y within MB
                  int img_block_x,       //!< location of block X, multiples of 4
                  int img_block_y)       //!< location of block Y, multiples of 4
{
    VideoParameters *p_Vid = currMB->p_Vid;
    byte predmode = p_Vid->ipredmode[img_block_y][img_block_x];
    currMB->ipmode_DPCM = predmode; //For residual DPCM
    int BitDepth = pl ? p_Vid->bitdepth_chroma : p_Vid->bitdepth_luma;
    Slice *currSlice = currMB->p_Slice;
    imgpel *pred = &currSlice->mb_pred[pl][joff][ioff];
    imgpel pix[9 * 9];
    bool available[4];
    neighbouring_samples(pix, available, currMB, pl, ioff, joff);

    switch (predmode) {
    case Intra_4x4_Vertical:
        if (!available[1]) {
            printf ("warning: Intra_4x4_Vertical prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_Horizontal:
        if (!available[0]) {
            printf ("warning: Intra_4x4_Horizontal prediction mode not allowed at mb %d\n",(int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_DC:
        break;
    case Intra_4x4_Diagonal_Down_Left:
        if (!available[1]) {
            printf ("warning: Intra_4x4_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_Diagonal_Down_Right:
        if (!available[3]) {
            printf ("warning: Intra_4x4_Diagonal_Down_Right prediction mode not allowed at mb %d\n",(int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_Vertical_Right:
        if (!available[3]) {
            printf ("warning: Intra_4x4_Vertical_Right prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_Horizontal_Down:  
        if (!available[3]) {
            printf ("warning: Intra_4x4_Horizontal_Down prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_Vertical_Left:
        if (!available[1]) {
            printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_4x4_Horizontal_Up:
        if (!available[0]) {
            printf ("warning: Intra_4x4_Horizontal_Up prediction mode not allowed at mb %d\n",(int) currSlice->current_mb_nr);
            return DECODING_OK;
        }    
        break;
    default:
        printf("Error: illegal intra_4x4 prediction mode: %d\n", (int) predmode);
        return SEARCH_SYNC;
        break;
    }

    switch (predmode) {
    case Intra_4x4_Vertical:
        intra4x4_vert_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_Horizontal:
        intra4x4_hor_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_DC:
        intra4x4_dc_pred_normal(pred, pix, available, BitDepth);
        break;
    case Intra_4x4_Diagonal_Down_Left:
        intra4x4_diag_down_left_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_Diagonal_Down_Right:
        intra4x4_diag_down_right_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_Vertical_Right:
        intra4x4_vert_right_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_Horizontal_Down:  
        intra4x4_hor_down_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_Vertical_Left:
        intra4x4_vert_left_pred_normal(pred, pix, available);
        break;
    case Intra_4x4_Horizontal_Up:
        intra4x4_hor_up_pred_normal(pred, pix, available);
        break;
    }

    return DECODING_OK;
}
