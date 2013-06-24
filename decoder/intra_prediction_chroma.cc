/*!
 *************************************************************************************
 * \file intra_chroma_pred.c
 *
 * \brief
 *    Functions for intra chroma prediction
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
#include "transform.h"
#include "neighbour.h"
#include "image.h"


static void neighbouring_samples(imgpel *pred, bool *available, Macroblock *currMB, ColorPlane pl, int xO, int yO)
{
    Slice *currSlice = currMB->p_Slice;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    int *mb_size = (pl) ? p_Vid->mb_size[IS_CHROMA] : p_Vid->mb_size[IS_LUMA];
    //PixelPos pix_a;
    PixelPos pix_b, pix_d, pix_x[16];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->mb_aff_frame_flag;

    if (MbaffFrameFlag == 0) {
        //getNonAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < p_Vid->mb_cr_size_y; i++)
            getNonAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getNonAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getNonAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    } else {
        //getAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < p_Vid->mb_cr_size_y; i++)
            getAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    }

    if (constrained_intra_pred_flag) {
        available[0] = 1;
        //available[0] = pix_a.available ? currSlice->intra_block[pix_a.mb_addr] : 0;
        for (int i = 0; i < p_Vid->mb_cr_size_y/2; i++)
            available[0] &= pix_x[i].available ? currSlice->intra_block[pix_x[i].mb_addr] : 0;
        available[1] = pix_b.available ? currSlice->intra_block[pix_b.mb_addr] : 0;
        available[2] = 1;
        for (int i = p_Vid->mb_cr_size_y/2; i < p_Vid->mb_cr_size_y; i++)
            available[2] &= pix_x[i].available ? currSlice->intra_block[pix_x[i].mb_addr] : 0;
        available[3] = pix_d.available ? currSlice->intra_block[pix_d.mb_addr] : 0;
    } else {
        available[0] = pix_x[0].available; //pix_a.available;
        available[1] = pix_b.available;
        available[2] = pix_x[0].available; //pix_a.available;
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
        for (int y = 0; y < p_Vid->mb_cr_size_y/2; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[2]) {
        for (int y = p_Vid->mb_cr_size_y/2; y < p_Vid->mb_cr_size_y; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[1]) {
        imgpel *pix = &img[pix_b.pos_y][pix_b.pos_x];
        for (int x = 0; x < p_Vid->mb_cr_size_x; x++)
            p(x, -1) = pix[x];
    }
#undef p
}



static inline int InverseRasterScan(int index, int dx, int dy, int width, int xy)
{
    if (xy == 0)
        return (index % (width / dx)) * dx;
    else
        return (index / (width / dx)) * dy;
}

static inline int Clip3(int low, int high, int x)
{
  x = imax(x, low);
  x = imin(x, high);

  return x;
}

#define predCb(x,y) (pred[0][(y) * 16 + (x)])
#define predCr(x,y) (pred[1][(y) * 16 + (x)])
#define pcb(x,y) (pix[0][((y) + 1) * 17 + ((x) + 1)])
#define pcr(x,y) (pix[1][((y) + 1) * 17 + ((x) + 1)])


static void intra4x4_dc_pred_chroma(imgpel *pred[2], imgpel pix[2][17*17], bool *available, int xO, int yO, int BitDepth)
{
    bool block_available_left = available[0];
    bool block_available_up   = available[1];

    int shift = 1;
    int sum0 = 0, sum1 = 0;
    if (block_available_left) {
        for (int y = 0; y < 4; y++) {
            sum0 += pcb(-1, y + yO);
            sum1 += pcr(-1, y + yO);
        }
        sum0 += 2;
        sum1 += 2;
        shift++;
    }
    if (block_available_up) {
        for (int x = 0; x < 4; x++) {
            sum0 += pcb(x + xO, -1);
            sum1 += pcr(x + xO, -1);
        }
        sum0 += 2;
        sum1 += 2;
        shift++;
    }
    sum0 >>= shift;
    sum1 >>= shift;
    if (!block_available_left && !block_available_up) {
        sum0 = 1 << (BitDepth - 1);
        sum1 = 1 << (BitDepth - 1);
    }

    for (int y = 0; y < 4; y++) {
        for (int x = 0; x < 4; x++) {
            predCb(x + xO, y + yO) = sum0;
            predCr(x + xO, y + yO) = sum1;
        }
    }
}

static void intrapred_chroma_dc(imgpel *pred[2], imgpel pix[2][17*17], bool *available, int size[2], int BitDepth)
{
    int ChromaArrayType = size[0] == 8 && size[1] == 8 ? 1
                        : size[0] == 16 && size[1] == 16 ? 3
                        : 2;

    for (int chroma4x4BlkIdx = 0; chroma4x4BlkIdx < (1 << (ChromaArrayType + 1)); chroma4x4BlkIdx++) {
        int xO = InverseRasterScan(chroma4x4BlkIdx / 4, 8, 8, size[0], 0)
               + InverseRasterScan(chroma4x4BlkIdx % 4, 4, 4, size[0], 0);
        int yO = InverseRasterScan(chroma4x4BlkIdx / 4, 8, 8, size[0], 1)
               + InverseRasterScan(chroma4x4BlkIdx % 4, 4, 4, size[0], 1);

        bool avail[4] = { 0 };
        if ((xO == 0 && yO == 0) || (xO > 0 && yO > 0)) {
            avail[0] = yO > 0 ? available[2] : available[0];
            avail[1] = available[1];
        } else if (xO > 0 && yO == 0) {
            avail[0] = available[1] ? 0 : available[0];
            avail[1] = available[1];
        } else if (xO == 0 && yO > 0) {
            avail[0] = available[2];
            avail[1] = available[2] ? 0 : available[1];
        }

        intra4x4_dc_pred_chroma(pred, pix, avail, xO, yO, BitDepth);
    }
}

static void intrapred_chroma_hor(imgpel *pred[2], imgpel pix[2][17*17], bool *available, int size[2])
{
    for (int y = 0; y < size[1]; y++) {
        for (int x = 0; x < size[0]; x++) {
            predCb(x, y) = pcb(-1, y);
            predCr(x, y) = pcr(-1, y);
        }
    }
}

static void intrapred_chroma_ver(imgpel *pred[2], imgpel pix[2][17*17], bool *available, int size[2])
{
    for (int y = 0; y < size[1]; y++) {
        for (int x = 0; x < size[0]; x++) {
            predCb(x, y) = pcb(x, -1);
            predCr(x, y) = pcr(x, -1);
        }
    }
}

static void intrapred_chroma_plane(imgpel *pred[2], imgpel pix[2][17*17], bool *available, int size[2], int BitDepth)
{
    int ChromaArrayType = size[0] == 8 && size[1] == 8 ? 1
                        : size[0] == 16 && size[1] == 16 ? 3
                        : 2;

    int xCF = ChromaArrayType == 3 ? 4 : 0;
    int yCF = ChromaArrayType != 1 ? 4 : 0;

    int H0 = 0, H1 = 0;
    for (int x = 0; x < 4 + xCF; x++) {
        H0 += (x + 1) * (pcb(4 + xCF + x, -1) - pcb(2 + xCF - x, -1));
        H1 += (x + 1) * (pcr(4 + xCF + x, -1) - pcr(2 + xCF - x, -1));
    }
    int V0 = 0, V1 = 0;
    for (int y = 0; y < 4 + yCF; y++) {
        V0 += (y + 1) * (pcb(-1, 4 + yCF + y) - pcb(-1, 2 + yCF - y));
        V1 += (y + 1) * (pcr(-1, 4 + yCF + y) - pcr(-1, 2 + yCF - y));
    }

    int a0 = 16 * (pcb(-1, size[1] - 1) + pcb(size[0] - 1, -1));
    int a1 = 16 * (pcr(-1, size[1] - 1) + pcr(size[0] - 1, -1));
    int b0 = ((34 - 29 * (ChromaArrayType == 3 ? 1 : 0)) * H0 + 32) >> 6;
    int b1 = ((34 - 29 * (ChromaArrayType == 3 ? 1 : 0)) * H1 + 32) >> 6;
    int c0 = ((34 - 29 * (ChromaArrayType != 1 ? 1 : 0)) * V0 + 32) >> 6;
    int c1 = ((34 - 29 * (ChromaArrayType != 1 ? 1 : 0)) * V1 + 32) >> 6;

#define Clip1(x) (Clip3(0, (1 << BitDepth) - 1, x))
    for (int y = 0; y < size[1]; y++) {
        for (int x = 0; x < size[0]; x++) {
            predCb(x, y) = Clip1((a0 + b0 * (x - 3 - xCF) + c0 * (y - 3 - yCF) + 16) >> 5);
            predCr(x, y) = Clip1((a1 + b1 * (x - 3 - xCF) + c1 * (y - 3 - yCF) + 16) >> 5);
        }
    }
#undef Clip1
}

#undef pcb
#undef pcr
#undef predCb
#undef predCr

/*!
 ************************************************************************
 * \brief
 *    Chroma Intra prediction. Note that many operations can be moved
 *    outside since they are repeated for both components for no reason.
 ************************************************************************
 */
void intra_pred_chroma(Macroblock *currMB)
{
    VideoParameters *p_Vid = currMB->p_Vid;

    int BitDepth = p_Vid->bitdepth_chroma;
    Slice *currSlice = currMB->p_Slice;
    imgpel *pred[2];
    int size[2];
    imgpel pix[2][17 * 17];
    bool available[2][4];
    pred[0] = &currSlice->mb_pred[1][0][0];
    pred[1] = &currSlice->mb_pred[2][0][0];
    size[0] = p_Vid->mb_cr_size_x;
    size[1] = p_Vid->mb_cr_size_y;
    neighbouring_samples(pix[0], available[0], currMB, PLANE_U, 0, 0);
    neighbouring_samples(pix[1], available[1], currMB, PLANE_V, 0, 0);

    switch (currMB->c_ipred_mode) {
    case Intra_Chroma_DC:  
        break;
    case Intra_Chroma_Horizontal: 
        if (!available[0][0]) {
            error("unexpected HOR_PRED_8 chroma intra prediction mode",-1);
            return;
        }
        break;
    case Intra_Chroma_Vertical: 
        if (!available[0][1]) {
            error("unexpected VERT_PRED_8 chroma intra prediction mode",-1);
            return;
        }
        break;
    case Intra_Chroma_Plane: 
        if (!available[0][3]) {
            error("unexpected PLANE_8 chroma intra prediction mode",-1);
            return;
        }
        break;
    default:
        error("illegal chroma intra prediction mode", 600);
        return;
    }

    switch (currMB->c_ipred_mode) {
    case Intra_Chroma_DC:  
        intrapred_chroma_dc(pred, pix, available[0], size, BitDepth);
        break;
    case Intra_Chroma_Horizontal: 
        intrapred_chroma_hor(pred, pix, available[0], size);
        break;
    case Intra_Chroma_Vertical: 
        intrapred_chroma_ver(pred, pix, available[0], size);
        break;
    case Intra_Chroma_Plane: 
        intrapred_chroma_plane(pred, pix, available[0], size, BitDepth);
        break;
    }
}
