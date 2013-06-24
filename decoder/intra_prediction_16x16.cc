/*!
 *************************************************************************************
 * \file intra16x16_pred.c
 *
 * \brief
 *    Functions for intra 16x16 prediction
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
#include "intra_prediction.h"
#include "neighbour.h"
#include "image.h"


static void neighbouring_samples(imgpel *pred, bool *available, Macroblock *currMB, ColorPlane pl, int xO, int yO)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    int *mb_size = (pl) ? p_Vid->mb_size[IS_CHROMA] : p_Vid->mb_size[IS_LUMA];
    PixelPos pix_a, pix_b, pix_d, pix_x[16];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->mb_aff_frame_flag;

    if (MbaffFrameFlag == 0) {
        getNonAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < 16; i++)
            getNonAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getNonAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getNonAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    } else {
        getAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < 16; i++)
            getAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    }

    if (constrained_intra_pred_flag) {
        available[0] = pix_a.available ? currSlice->intra_block[pix_a.mb_addr] : 0;
        for (int i = 0; i < 16; i++)
            available[0] &= pix_x[i].available ? currSlice->intra_block[pix_x[i].mb_addr] : 0;
        available[1] = pix_b.available ? currSlice->intra_block[pix_b.mb_addr] : 0;
        available[2] = 0;
        available[3] = pix_d.available ? currSlice->intra_block[pix_d.mb_addr] : 0;
    } else {
        available[0] = pix_a.available;
        available[1] = pix_b.available;
        available[2] = 0;
        available[3] = pix_d.available;
    }

#define p(x,y) (pred[((y) + 1) * 17 + ((x) + 1)])
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
        for (int y = 0; y < 16; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[1]) {
        imgpel *pix = &img[pix_b.pos_y][pix_b.pos_x];
        for (int x = 0; x < 16; x++)
            p(x, -1) = pix[x];
    }
#undef p
}

static inline int Clip3(int low, int high, int x)
{
  x = imax(x, low);
  x = imin(x, high);

  return x;
}


#define predL(x,y) (pred[(y) * 16 + (x)])
#define p(x,y) (pix[((y) + 1) * 17 + ((x) + 1)])

static void intra16x16_vert_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++) {
            predL(x, y) = p(x, -1);
        }
    }
}

static void intra16x16_hor_pred_normal(imgpel *pred, imgpel *pix, bool *available)
{
    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++) {
            predL(x, y) = p(-1, y);
        }
    }
}

static void intra16x16_dc_pred_normal(imgpel *pred, imgpel *pix, bool *available, int BitDepth)
{
    bool block_available_left = available[0];
    bool block_available_up   = available[1];

    int shift = 3;
    int sum = 0;
    if (block_available_left) {
        for (int y = 0; y < 16; y++) {
            sum += p(-1, y);
        }
        sum += 8;
        shift++;
    }
    if (block_available_up) {
        for (int x = 0; x < 16; x++) {
            sum += p(x, -1);
        }
        sum += 8;
        shift++;
    }
    sum >>= shift;
    if (!block_available_left && !block_available_up)
        sum = 1 << (BitDepth - 1);

    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++) {
            predL(x, y) = sum;
        }
    }
}

static void intra16x16_plane_pred_normal(imgpel *pred, imgpel *pix, bool *available, int BitDepth)
{
    int H = 0;
    for (int x = 0; x < 8; x++) {
        H += (x + 1) * (p(8 + x, -1) - p(6 - x, -1));
    }
    int V = 0;
    for (int y = 0; y < 8; y++) {
        V += (y + 1) * (p(-1, 8 + y) - p(-1, 6 - y));
    }

    int a = 16 * (p(-1, 15) + p(15, -1));
    int b = (5 * H + 32) >> 6;
    int c = (5 * V + 32) >> 6;

#define Clip1(x) (Clip3(0, (1 << BitDepth) - 1, x))
    for (int y = 0; y < 16; y++) {
        for (int x = 0; x < 16; x++) {
            predL(x, y) = Clip1((a + b * (x - 7) + c * (y - 7) + 16) >> 5);
        }
    }
#undef Clip1
}

#undef p
#undef predL

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
int intra_pred_16x16(Macroblock *currMB,  //!< Current Macroblock
                           ColorPlane pl,       //!< Current colorplane (for 4:4:4)                         
                           int predmode)        //!< prediction mode
{
    VideoParameters *p_Vid = currMB->p_Vid;

    int BitDepth = pl ? p_Vid->bitdepth_chroma : p_Vid->bitdepth_luma;
    Slice *currSlice = currMB->p_Slice;
    imgpel *pred = &currSlice->mb_pred[pl][0][0];
    imgpel pix[17 * 17];
    bool available[4];
    neighbouring_samples(pix, available, currMB, pl, 0, 0);

    switch (predmode) {
    case Intra_16x16_Vertical:
        if (!available[1]) {
            error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);
            return DECODING_OK;
        }
        break;
    case Intra_16x16_Horizontal:
        if (!available[0]) {
            error ("invalid 16x16 intra pred Mode HOR_PRED_16",500);
            return DECODING_OK;
        }
        break;
    case Intra_16x16_DC:
        break;
    case Intra_16x16_Plane:
        if (!available[3]) {
            error ("invalid 16x16 intra pred Mode PLANE_16",500);
            return DECODING_OK;
        }
        break;
    default:
        printf("illegal 16x16 intra prediction mode input: %d\n",predmode);
        return SEARCH_SYNC;
    }

    switch (predmode) {
    case Intra_16x16_Vertical:
        intra16x16_vert_pred_normal(pred, pix, available);
        break;
    case Intra_16x16_Horizontal:
        intra16x16_hor_pred_normal(pred, pix, available);
        break;
    case Intra_16x16_DC:
        intra16x16_dc_pred_normal(pred, pix, available, BitDepth);
        break;
    case Intra_16x16_Plane:
        intra16x16_plane_pred_normal(pred, pix, available, BitDepth);
        break;
    }

    return DECODING_OK;
}
