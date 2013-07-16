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


static void neighbouring_samples_4x4(imgpel *pred, bool *available, mb_t *currMB, ColorPlane pl, int xO, int yO)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    PixelPos pix_a, pix_b, pix_c, pix_d, pix_x[4];

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { mb_cr_size_x, mb_cr_size_y }
    };
    int *mb_size = mb_size_xy[pl == 0 ? IS_LUMA : IS_CHROMA];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->MbaffFrameFlag;

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
        //int dy = MbaffFrameFlag == 1 && p_Vid->mb_data[pix_a.mb_addr].mb_field_decoding_flag ? 2 : 1;
//        imgpel *pix = &img[pix_a.pos_y][pix_a.pos_x];
//        if (MbaffFrameFlag == 1 && !currMB->mb_field_decoding_flag && (currMB->mbAddrX & 1))
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

static void neighbouring_samples_8x8(imgpel *pred, bool *available, mb_t *currMB, ColorPlane pl, int xO, int yO)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    PixelPos pix_a, pix_b, pix_c, pix_d, pix_x[8];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->MbaffFrameFlag;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { mb_cr_size_x, mb_cr_size_y }
    };
    int *mb_size = mb_size_xy[pl == 0 ? IS_LUMA : IS_CHROMA];

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
        //int dy = MbaffFrameFlag == 1 && p_Vid->mb_data[pix_a.mb_addr].mb_field_decoding_flag ? 2 : 1;
//        imgpel *pix = &img[pix_a.pos_y][pix_a.pos_x];
//        if (MbaffFrameFlag == 1 && !currMB->mb_field_decoding_flag && (currMB->mbAddrX & 1))
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

static void neighbouring_samples_16x16(imgpel *pred, bool *available, mb_t *currMB, ColorPlane pl, int xO, int yO)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    PixelPos pix_a, pix_b, pix_d, pix_x[16];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->MbaffFrameFlag;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { mb_cr_size_x, mb_cr_size_y }
    };
    int *mb_size = mb_size_xy[pl == 0 ? IS_LUMA : IS_CHROMA];

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
        //int dy = MbaffFrameFlag == 1 && p_Vid->mb_data[pix_a.mb_addr].mb_field_decoding_flag ? 2 : 1;
//        imgpel *pix = &img[pix_a.pos_y][pix_a.pos_x];
//        if (MbaffFrameFlag == 1 && !currMB->mb_field_decoding_flag && (currMB->mbAddrX & 1))
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

static void neighbouring_samples_chroma(imgpel *pred, bool *available, mb_t *currMB, ColorPlane pl, int xO, int yO)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    VideoParameters *p_Vid = currMB->p_Vid;
    imgpel **img = (pl) ? currSlice->dec_picture->imgUV[pl - 1] : currSlice->dec_picture->imgY;
    //PixelPos pix_a;
    PixelPos pix_b, pix_d, pix_x[16];

    bool constrained_intra_pred_flag = p_Vid->active_pps->constrained_intra_pred_flag;
    bool MbaffFrameFlag = currSlice->MbaffFrameFlag;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { mb_cr_size_x, mb_cr_size_y }
    };
    int *mb_size = mb_size_xy[pl == 0 ? IS_LUMA : IS_CHROMA];

    if (MbaffFrameFlag == 0) {
        //getNonAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < mb_cr_size_y; i++)
            getNonAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getNonAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getNonAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    } else {
        //getAffNeighbour(currMB, xO - 1, yO    , mb_size, &pix_a);
        for (int i = 0; i < mb_cr_size_y; i++)
            getAffNeighbour(currMB, xO - 1, yO + i, mb_size, &pix_x[i]);
        getAffNeighbour(currMB, xO    , yO - 1, mb_size, &pix_b);
        getAffNeighbour(currMB, xO - 1, yO - 1, mb_size, &pix_d);
    }

    if (constrained_intra_pred_flag) {
        available[0] = 1;
        //available[0] = pix_a.available ? currSlice->intra_block[pix_a.mb_addr] : 0;
        for (int i = 0; i < mb_cr_size_y/2; i++)
            available[0] &= pix_x[i].available ? currSlice->intra_block[pix_x[i].mb_addr] : 0;
        available[1] = pix_b.available ? currSlice->intra_block[pix_b.mb_addr] : 0;
        available[2] = 1;
        for (int i = mb_cr_size_y/2; i < mb_cr_size_y; i++)
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
        //int dy = MbaffFrameFlag == 1 && p_Vid->mb_data[pix_a.mb_addr].mb_field_decoding_flag ? 2 : 1;
//        imgpel *pix = &img[pix_a.pos_y][pix_a.pos_x];
//        if (MbaffFrameFlag == 1 && !currMB->mb_field_decoding_flag && (currMB->mbAddrX & 1))
//            pix -= 16 * width;
//        for (int y = 0; y < 4; y++)
//            p(-1, y) = pix[y * dy * width];
        for (int y = 0; y < mb_cr_size_y/2; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[2]) {
        for (int y = mb_cr_size_y/2; y < mb_cr_size_y; y++)
            p(-1, y) = img[pix_x[y].pos_y][pix_x[y].pos_x];
    }
    if (available[1]) {
        imgpel *pix = &img[pix_b.pos_y][pix_b.pos_x];
        for (int x = 0; x < mb_cr_size_x; x++)
            p(x, -1) = pix[x];
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
 ***********************************************************************
 * \brief
 *    makes and returns 4x4 intra prediction blocks 
 *
 * \return
 *    DECODING_OK   decoding of intra prediction mode was successful            \n
 *    SEARCH_SYNC   search next sync element as errors while decoding occured
 ***********************************************************************
 */
void intra_prediction_t::intra_pred_4x4(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    int block_pos[4][4] = {
        {  0,  1,  4,  5 },
        {  2,  3,  6,  7 },
        {  8,  9, 12, 13 },
        { 10, 11, 14, 15 }
    };
    uint8_t predmode = currMB->Intra4x4PredMode[block_pos[joff/4][ioff/4]];
    currMB->ipmode_DPCM = predmode; //For residual DPCM

    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int BitDepth = pl ? sps->BitDepthC : sps->BitDepthY;
    imgpel *pred = &currSlice->mb_pred[pl][joff][ioff];
    imgpel pix[9 * 9];
    bool available[4];
    neighbouring_samples_4x4(pix, available, currMB, pl, ioff, joff);

    switch (predmode) {
    case Intra_4x4_Vertical:
        if (!available[1])
            printf ("warning: Intra_4x4_Vertical prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_Horizontal:
        if (!available[0])
            printf ("warning: Intra_4x4_Horizontal prediction mode not allowed at mb %d\n",(int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_DC:
        break;
    case Intra_4x4_Diagonal_Down_Left:
        if (!available[1])
            printf ("warning: Intra_4x4_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_Diagonal_Down_Right:
        if (!available[3])
            printf ("warning: Intra_4x4_Diagonal_Down_Right prediction mode not allowed at mb %d\n",(int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_Vertical_Right:
        if (!available[3])
            printf ("warning: Intra_4x4_Vertical_Right prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_Horizontal_Down:  
        if (!available[3])
            printf ("warning: Intra_4x4_Horizontal_Down prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_Vertical_Left:
        if (!available[1])
            printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_4x4_Horizontal_Up:
        if (!available[0])
            printf ("warning: Intra_4x4_Horizontal_Up prediction mode not allowed at mb %d\n",(int) currSlice->current_mb_nr);
        break;
    default:
        printf("Error: illegal intra_4x4 prediction mode: %d\n", (int) predmode);
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
}

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
void intra_prediction_t::intra_pred_8x8(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    uint8_t predmode = currMB->Intra8x8PredMode[joff/8 * 2 + ioff/8];
    currMB->ipmode_DPCM = predmode;  //For residual DPCM

    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int BitDepth = pl ? sps->BitDepthC : sps->BitDepthY;
    imgpel *pred = &currSlice->mb_pred[pl][joff][ioff];
    imgpel pix[17 * 17];
    bool available[4];
    neighbouring_samples_8x8(pix, available, currMB, pl, ioff, joff);

    switch (predmode) {
    case Intra_8x8_Vertical:
        if (!available[1])
            printf ("warning: Intra_8x8_Vertical prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_Horizontal:
        if (!available[0])
            printf ("warning: Intra_8x8_Horizontal prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_DC:
        break;
    case Intra_8x8_Diagonal_Down_Left:
        if (!available[1])
            printf ("warning: Intra_8x8_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_Diagonal_Down_Right:
        if (!available[3])
            printf ("warning: Intra_8x8_Diagonal_Down_Right prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_Vertical_Right:
        if (!available[3])
            printf ("warning: Intra_8x8_Vertical_Right prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_Horizontal_Down:  
        if (!available[3])
            printf ("warning: Intra_8x8_Horizontal_Down prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_Vertical_Left:
        if (!available[1])
            printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    case Intra_8x8_Horizontal_Up:
        if (!available[0])
            printf ("warning: Intra_8x8_Horizontal_Up prediction mode not allowed at mb %d\n", (int) currSlice->current_mb_nr);
        break;
    default:
        printf("Error: illegal intra_8x8 prediction mode: %d\n", (int) predmode);
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
void intra_prediction_t::intra_pred_16x16(mb_t *currMB, ColorPlane pl, int ioff, int joff)
{
    currMB->ipmode_DPCM = currMB->Intra16x16PredMode;

    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int BitDepth = pl ? sps->BitDepthC : sps->BitDepthY;
    imgpel *pred = &currSlice->mb_pred[pl][0][0];
    imgpel pix[17 * 17];
    bool available[4];
    neighbouring_samples_16x16(pix, available, currMB, pl, 0, 0);

    switch (currMB->Intra16x16PredMode) {
    case Intra_16x16_Vertical:
        if (!available[1])
            error ("invalid 16x16 intra pred Mode VERT_PRED_16",500);
        break;
    case Intra_16x16_Horizontal:
        if (!available[0])
            error ("invalid 16x16 intra pred Mode HOR_PRED_16",500);
        break;
    case Intra_16x16_DC:
        break;
    case Intra_16x16_Plane:
        if (!available[3])
            error ("invalid 16x16 intra pred Mode PLANE_16",500);
        break;
    default:
        printf("illegal 16x16 intra prediction mode input: %d\n",currMB->Intra16x16PredMode);
        return;
    }

    switch (currMB->Intra16x16PredMode) {
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
}

/*!
 ************************************************************************
 * \brief
 *    Chroma Intra prediction. Note that many operations can be moved
 *    outside since they are repeated for both components for no reason.
 ************************************************************************
 */
void intra_prediction_t::intra_pred_chroma(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int BitDepth = sps->BitDepthC;
    imgpel *pred[2];
    int size[2];
    imgpel pix[2][17 * 17];
    bool available[2][4];
    pred[0] = &currSlice->mb_pred[1][0][0];
    pred[1] = &currSlice->mb_pred[2][0][0];
    size[0] = sps->MbWidthC;
    size[1] = sps->MbHeightC;
    neighbouring_samples_chroma(pix[0], available[0], currMB, PLANE_U, 0, 0);
    neighbouring_samples_chroma(pix[1], available[1], currMB, PLANE_V, 0, 0);

    switch (currMB->intra_chroma_pred_mode) {
    case Intra_Chroma_DC:  
        break;
    case Intra_Chroma_Horizontal: 
        if (!available[0][0])
            error("unexpected HOR_PRED_8 chroma intra prediction mode",-1);
        break;
    case Intra_Chroma_Vertical: 
        if (!available[0][1])
            error("unexpected VERT_PRED_8 chroma intra prediction mode",-1);
        break;
    case Intra_Chroma_Plane: 
        if (!available[0][3])
            error("unexpected PLANE_8 chroma intra prediction mode",-1);
        break;
    default:
        error("illegal chroma intra prediction mode", 600);
    }

    switch (currMB->intra_chroma_pred_mode) {
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
