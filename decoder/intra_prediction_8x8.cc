/*!
 *************************************************************************************
 * \file intra8x8_pred.c
 *
 * \brief
 *    Functions for intra 8x8 prediction
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
#include "intra_prediction.h"
#include "neighbour.h"
#include "image.h"



static void neighbouring_samples(imgpel *pred, bool *available, Macroblock *currMB, ColorPlane pl, int xO, int yO)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    int *mb_size = (pl) ? p_Vid->mb_size[IS_CHROMA] : p_Vid->mb_size[IS_LUMA];
    PixelPos pix_a, pix_b, pix_c, pix_d, pix_x[8];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->mb_aff_frame_flag;

    if (MbaffFrameFlag == 0) {
        getNonAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < 8; i++)
            getNonAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getNonAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getNonAffNeighbour(currMB, xO + 8, yO - 1, mb_size, &pix_c);
        getNonAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    } else {
        getAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < 8; i++)
            getAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getAffNeighbour(currMB, xO + 8, yO - 1, mb_size, &pix_c);
        getAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    }
    pix_c.available = pix_c.available && !(xO == 8 && yO == 8);

    if (constrained_intra_pred_flag) {
        available[0] = pix_a.available ? currSlice->intra_block[pix_a.mb_addr] : 0;
        for (int i = 0; i < 8; i++)
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

    imgpel predLF[17 * 17];

#define plf(x,y) (pred[((y) + 1) * 17 + ((x) + 1)])
#define p(x,y) (predLF[((y) + 1) * 17 + ((x) + 1)])
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
        for (int y = 0; y < 8; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[1]) {
        imgpel *pix = &img[pix_b.pos_y][pix_b.pos_x];
        for (int x = 0; x < 8; x++)
            p(x, -1) = pix[x];
        pix = &img[pix_c.pos_y][pix_c.pos_x-8];
        for (int x = 8; x < 16; x++)
            p(x, -1) = available[2] ? pix[x] : p(7, -1);
        available[2] = available[1];
    }

    if (available[1]) {
        if (available[3])
            plf(0, -1) = (p(-1, -1) + 2 * p(0, -1) + p(1, -1) + 2) >> 2;
        else
            plf(0, -1) = (3 * p(0, -1) + p(1, -1) + 2) >> 2;
        for (int x = 1; x < 15; x++)
            plf(x, -1) = (p(x - 1, -1) + 2 * p(x, -1) + p(x + 1, -1) + 2) >> 2;
        plf(15, -1) = (p(14, -1) + 3 * p(15, -1) + 2) >> 2;
    }
    if (available[3]) {
        if (available[0] && available[1])
            plf(-1, -1) = (p(0, -1) + 2 * p(-1, -1) + p(-1, 0) + 2) >> 2;
        else if (available[1])
            plf(-1, -1) = (3 * p(-1, -1) + p(0, -1) + 2) >> 2;
        else if (available[0])
            plf(-1, -1) = (3 * p(-1, -1) + p(-1, 0) + 2) >> 2;
        else
            plf(-1, -1) = p(-1, -1);
    }
    if (available[0]) {
        if (available[3])
            plf(-1, 0) = (p(-1, -1) + 2 * p(-1, 0) + p(-1, 1) + 2) >> 2;
        else
            plf(-1, 0) = (3 * p(-1, 0) + p(-1, 1) + 2) >> 2;
        for (int y = 1; y < 7; y++)
            plf(-1, y) = (p(-1, y - 1) + 2 * p(-1, y) + p(-1, y + 1) + 2) >> 2;
        plf(-1, 7) = (p(-1, 6) + 3 * p(-1, 7) + 2) >> 2;
    }
#undef p
#undef plf
}


#define pred8x8L(x,y) (pred[(y) * 16 + (x)])
#define plf(x,y) (pix[((y) + 1) * 17 + ((x) + 1)])

static void intra8x8_vert_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            pred8x8L(x, y) = plf(x, -1);
        }
    }
}

static void intra8x8_hor_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            pred8x8L(x, y) = plf(-1, y);
        }
    }
}

static void intra8x8_dc_pred(imgpel *pred, imgpel *pix, bool *available, int BitDepth)
{
    bool block_available_left = available[0];
    bool block_available_up   = available[1];

    int shift = 2;
    int sum = 0;
    if (block_available_left) {
        for (int y = 0; y < 8; y++) {
            sum += plf(-1, y);
        }
        sum += 4;
        shift++;
    }
    if (block_available_up) {
        for (int x = 0; x < 8; x++) {
            sum += plf(x, -1);
        }
        sum += 4;
        shift++;
    }
    sum >>= shift;
    if (!block_available_left && !block_available_up)
        sum = 1 << (BitDepth - 1);

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            pred8x8L(x, y) = sum;
        }
    }
}

static void intra8x8_diag_down_left_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            if (x == 7 && y == 7)
                pred8x8L(x, y) = (plf(x + y, -1) + 3 * plf(x + y + 1, -1) + 2) >> 2;
            else
                pred8x8L(x, y) = (plf(x + y, -1) + 2 * plf(x + y + 1, -1) + plf(x + y + 2, -1) + 2) >> 2;
        }
    }
}

static void intra8x8_diag_down_right_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            if (x > y)
                pred8x8L(x, y) = (plf(x - y - 2, -1) + 2 * plf(x - y - 1, -1) + plf(x - y, -1) + 2) >> 2;
            else if (x < y)
                pred8x8L(x, y) = (plf(-1, y - x - 2) + 2 * plf(-1, y - x - 1) + plf(-1, y - x) + 2) >> 2;
            else
                pred8x8L(x, y) = (plf(0, -1) + 2 * plf(-1, -1) + plf(-1, 0) + 2) >> 2;
        }
    }
}

static void intra8x8_vert_right_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            int zVR = 2 * x - y;
            if (zVR >= 0 && (zVR % 2) == 0)
                pred8x8L(x, y) = (plf(x - (y >> 1) - 1, -1) + plf(x - (y >> 1), -1) + 1) >> 1;
            else if (zVR >= 0 && (zVR % 2) == 1)
                pred8x8L(x, y) = (plf(x - (y >> 1) - 2, -1) + 2 * plf(x - (y >> 1) - 1, -1) + plf(x - (y >> 1), -1) + 2) >> 2;
            else if (zVR == -1)
                pred8x8L(x, y) = (plf(-1, 0) + 2 * plf(-1, -1) + plf(0, -1) + 2) >> 2;
            else // zVR < -1
                pred8x8L(x, y) = (plf(-1, y - 2 * x - 1) + 2 * plf(-1, y - 2 * x - 2) + plf(-1, y - 2 * x - 3) + 2) >> 2;
        }
    }
}

static void intra8x8_hor_down_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            int zHD = 2 * y - x;
            if (zHD >= 0 && (zHD % 2) == 0)
                pred8x8L(x, y) = (plf(-1, y - (x >> 1) - 1) + plf(-1, y - (x >> 1)) + 1) >> 1;
            else if (zHD >= 0 && (zHD % 2) == 1)
                pred8x8L(x, y) = (plf(-1, y - (x >> 1) - 2) + 2 * plf(-1, y - (x >> 1) - 1) + plf(-1, y - (x >> 1)) + 2) >> 2;
            else if (zHD == -1)
                pred8x8L(x, y) = (plf(-1, 0) + 2 * plf(-1, -1) + plf(0, -1) + 2) >> 2;
            else // zHD < -1
                pred8x8L(x, y) = (plf(x - 2 * y - 1, -1) + 2 * plf(x - 2 * y - 2, -1) + plf(x - 2 * y - 3, -1) + 2) >> 2;
        }
    }
}

static void intra8x8_vert_left_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            if ((y % 2) == 0)
                pred8x8L(x, y) = (plf(x + (y >> 1), -1) + plf(x + (y >> 1) + 1, -1) + 1) >> 1;
            else // (y % 2) == 1
                pred8x8L(x, y) = (plf(x + (y >> 1), -1) + 2 * plf(x + (y >> 1) + 1, -1) + plf(x + (y >> 1) + 2, -1) + 2) >> 2;
        }
    }
}

static void intra8x8_hor_up_pred(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            int zHU = x + 2 * y;
            if (zHU < 13 && (zHU % 2) == 0)
                pred8x8L(x, y) = (plf(-1, y + (x >> 1)) + plf(-1, y + (x >> 1) + 1) + 1) >> 1;
            else if (zHU < 13 && (zHU % 2) == 1)
                pred8x8L(x, y) = (plf(-1, y + (x >> 1)) + 2 * plf(-1, y + (x >> 1) + 1) + plf(-1, y + (x >> 1) + 2) + 2) >> 2;
            else if (zHU == 13)
                pred8x8L(x, y) = (plf(-1, 6) + 3 * plf(-1, 7) + 2) >> 2;
            else // zHU > 13
                pred8x8L(x, y) = plf(-1, 7);
        }
    }
}

#undef plf
#undef pred8x8L

/*!
 ************************************************************************
 * \brief
 *    Make intra 8x8 prediction according to all 9 prediction modes.
 *    The routine uses left and upper neighbouring points from
 *    previous coded blocks to do this (if available). Notice that
 *    inaccessible neighbouring points are signalled with a negative
 *    value in the predmode array .
 *
 *  \par Input:
 *     Starting point of current 8x8 block image position
 *
 ************************************************************************
 */
int intra_pred_8x8(Macroblock *currMB,    //!< Current Macroblock
                        ColorPlane pl,         //!< Current color plane
                        int ioff,              //!< ioff
                        int joff)              //!< joff

{  
    int block_x = (currMB->block_x) + (ioff >> 2);
    int block_y = (currMB->block_y) + (joff >> 2);
    byte predmode = currMB->p_Slice->ipredmode[block_y][block_x];
    currMB->ipmode_DPCM = predmode;  //For residual DPCM

    VideoParameters *p_Vid = currMB->p_Vid;
    int BitDepth = pl ? p_Vid->bitdepth_chroma : p_Vid->bitdepth_luma;
    Slice *currSlice = currMB->p_Slice;
    imgpel *pred = &currSlice->mb_pred[pl][joff][ioff];
    imgpel pix[17 * 17];
    bool available[4];
    neighbouring_samples(pix, available, currMB, pl, ioff, joff);

    switch (predmode) {
    case Intra_8x8_Vertical:
        if (!available[1]) {
            printf ("warning: Intra_8x8_Vertical prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_Horizontal:
        if (!available[0]) {
            printf ("warning: Intra_8x8_Horizontal prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_DC:
        break;
    case Intra_8x8_Diagonal_Down_Left:
        if (!available[1]) {
            printf ("warning: Intra_8x8_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_Diagonal_Down_Right:
        if (!available[3]) {
            printf ("warning: Intra_8x8_Diagonal_Down_Right prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_Vertical_Right:
        if (!available[3]) {
            printf ("warning: Intra_8x8_Vertical_Right prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_Horizontal_Down:  
        if (!available[3]) {
            printf ("warning: Intra_8x8_Horizontal_Down prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_Vertical_Left:
        if (!available[1]) {
            printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    case Intra_8x8_Horizontal_Up:
        if (!available[0]) {
            printf ("warning: Intra_8x8_Horizontal_Up prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
            return DECODING_OK;
        }
        break;
    default:
        printf("Error: illegal intra_8x8 prediction mode: %d\n", (int) predmode);
        return SEARCH_SYNC;
        break;
    }

    switch (predmode) {
    case Intra_8x8_Vertical:
        intra8x8_vert_pred(pred, pix, available);
        break;
    case Intra_8x8_Horizontal:
        intra8x8_hor_pred(pred, pix, available);
        break;
    case Intra_8x8_DC:
        intra8x8_dc_pred(pred, pix, available, BitDepth);
        break;
    case Intra_8x8_Diagonal_Down_Left:
        intra8x8_diag_down_left_pred(pred, pix, available);
        break;
    case Intra_8x8_Diagonal_Down_Right:
        intra8x8_diag_down_right_pred(pred, pix, available);
        break;
    case Intra_8x8_Vertical_Right:
        intra8x8_vert_right_pred(pred, pix, available);
        break;
    case Intra_8x8_Horizontal_Down:  
        intra8x8_hor_down_pred(pred, pix, available);
        break;
    case Intra_8x8_Vertical_Left:
        intra8x8_vert_left_pred(pred, pix, available);
        break;
    case Intra_8x8_Horizontal_Up:
        intra8x8_hor_up_pred(pred, pix, available);
        break;
    }

    return DECODING_OK;
}
