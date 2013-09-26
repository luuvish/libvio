/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : intra_prediction.cpp
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"

#include "intra_prediction.h"


namespace vio  {
namespace h264 {


static void neighbouring_samples_4x4(imgpel* pred, bool* available, mb_t* mb, ColorPlane pl, int xO, int yO)
{
    slice_t* slice = mb->p_Slice;
    pps_t* pps = slice->active_pps;

    imgpel** img = pl ? slice->dec_picture->imgUV[pl - 1] : slice->dec_picture->imgY;

    loc_t locA[4];
    mb_t* mbA [4];
    pos_t posA[4];
    for (int i = 0; i < 4; ++i) {
        locA[i] = slice->neighbour.get_location(slice, mb->mbAddrX, {xO - 1, yO + i});
        mbA [i] = slice->neighbour.get_mb      (slice, locA[i]);
        posA[i] = slice->neighbour.get_blkpos  (slice, locA[i]);
        mbA [i] = mbA[i] && mbA[i]->slice_nr == mb->slice_nr ? mbA[i] : nullptr;
    }
    loc_t locB = slice->neighbour.get_location(slice, mb->mbAddrX, {xO    , yO - 1});
    loc_t locC = slice->neighbour.get_location(slice, mb->mbAddrX, {xO + 4, yO - 1});
    loc_t locD = slice->neighbour.get_location(slice, mb->mbAddrX, {xO - 1, yO - 1});
    mb_t* mbB  = slice->neighbour.get_mb      (slice, locB);
    mb_t* mbC  = slice->neighbour.get_mb      (slice, locC);
    mb_t* mbD  = slice->neighbour.get_mb      (slice, locD);
    pos_t posB = slice->neighbour.get_blkpos  (slice, locB);
    pos_t posC = slice->neighbour.get_blkpos  (slice, locC);
    pos_t posD = slice->neighbour.get_blkpos  (slice, locD);
    mbB = mbB && mbB->slice_nr == mb->slice_nr ? mbB : nullptr;
    mbC = mbC && mbC->slice_nr == mb->slice_nr ? mbC : nullptr;
    mbD = mbD && mbD->slice_nr == mb->slice_nr ? mbD : nullptr;

    mbC = mbC && !(xO == 4 && (yO == 4 || yO == 12)) ? mbC : nullptr;

    if (pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < 4; i++)
            available[0] &= mbA[i] && mbA[i]->is_intra_block ? 1 : 0;
        available[1] = mbB && mbB->is_intra_block ? 1 : 0;
        available[2] = mbC && mbC->is_intra_block ? 1 : 0;
        available[3] = mbD && mbD->is_intra_block ? 1 : 0;
    } else {
        available[0] = mbA[0] ? 1 : 0;
        available[1] = mbB    ? 1 : 0;
        available[2] = mbC    ? 1 : 0;
        available[3] = mbD    ? 1 : 0;
    }

#define p(x,y) (pred[((y) + 1) * 9 + ((x) + 1)])
    if (available[3]) {
        imgpel *pix = &img[posD.y][posD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < 4; ++y)
            p(-1, y) = img[posA[y].y][posA[y].x];
    }
    if (available[1]) {
        imgpel *pix = &img[posB.y][posB.x];
        for (int x = 0; x < 4; ++x)
            p(x, -1) = pix[x];
        pix = &img[posC.y][posC.x - 4];
        for (int x = 4; x < 8; x++)
            p(x, -1) = available[2] ? pix[x] : p(3, -1);
        available[2] = available[1];
    }
#undef p
}

static void neighbouring_samples_8x8(imgpel* pred, bool* available, mb_t* mb, ColorPlane pl, int xO, int yO)
{
    slice_t* slice = mb->p_Slice;
    pps_t* pps = slice->active_pps;

    imgpel** img = pl ? slice->dec_picture->imgUV[pl - 1] : slice->dec_picture->imgY;

    loc_t locA[8];
    mb_t* mbA [8];
    pos_t posA[8];
    for (int i = 0; i < 8; ++i) {
        locA[i] = slice->neighbour.get_location(slice, mb->mbAddrX, {xO - 1, yO + i});
        mbA [i] = slice->neighbour.get_mb      (slice, locA[i]);
        posA[i] = slice->neighbour.get_blkpos  (slice, locA[i]);
        mbA [i] = mbA[i] && mbA[i]->slice_nr == mb->slice_nr ? mbA[i] : nullptr;
    }
    loc_t locB = slice->neighbour.get_location(slice, mb->mbAddrX, {xO    , yO - 1});
    loc_t locC = slice->neighbour.get_location(slice, mb->mbAddrX, {xO + 8, yO - 1});
    loc_t locD = slice->neighbour.get_location(slice, mb->mbAddrX, {xO - 1, yO - 1});
    mb_t* mbB  = slice->neighbour.get_mb      (slice, locB);
    mb_t* mbC  = slice->neighbour.get_mb      (slice, locC);
    mb_t* mbD  = slice->neighbour.get_mb      (slice, locD);
    pos_t posB = slice->neighbour.get_blkpos  (slice, locB);
    pos_t posC = slice->neighbour.get_blkpos  (slice, locC);
    pos_t posD = slice->neighbour.get_blkpos  (slice, locD);
    mbB = mbB && mbB->slice_nr == mb->slice_nr ? mbB : nullptr;
    mbC = mbC && mbC->slice_nr == mb->slice_nr ? mbC : nullptr;
    mbD = mbD && mbD->slice_nr == mb->slice_nr ? mbD : nullptr;

    mbC = mbC && !(xO == 8 && yO == 8) ? mbC : nullptr;

    if (pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < 8; ++i)
            available[0] &= mbA[i] && mbA[i]->is_intra_block ? 1 : 0;
        available[1] = mbB && mbB->is_intra_block ? 1 : 0;
        available[2] = mbC && mbC->is_intra_block ? 1 : 0;
        available[3] = mbD && mbD->is_intra_block ? 1 : 0;
    } else {
        available[0] = mbA[0] ? 1 : 0;
        available[1] = mbB    ? 1 : 0;
        available[2] = mbC    ? 1 : 0;
        available[3] = mbD    ? 1 : 0;
    }

    imgpel predLF[17 * 17];

#define plf(x,y) (pred[((y) + 1) * 17 + ((x) + 1)])
#define p(x,y) (predLF[((y) + 1) * 17 + ((x) + 1)])
    if (available[3]) {
        imgpel *pix = &img[posD.y][posD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < 8; ++y)
            p(-1, y) = img[posA[y].y][posA[y].x];
    }
    if (available[1]) {
        imgpel *pix = &img[posB.y][posB.x];
        for (int x = 0; x < 8; ++x)
            p(x, -1) = pix[x];
        pix = &img[posC.y][posC.x - 8];
        for (int x = 8; x < 16; ++x)
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

static void neighbouring_samples_16x16(imgpel* pred, bool* available, mb_t* mb, ColorPlane pl, int xO, int yO)
{
    slice_t* slice = mb->p_Slice;
    pps_t* pps = slice->active_pps;

    imgpel** img = pl ? slice->dec_picture->imgUV[pl - 1] : slice->dec_picture->imgY;

    loc_t locA[16];
    mb_t* mbA [16];
    pos_t posA[16];
    for (int i = 0; i < 16; ++i) {
        locA[i] = slice->neighbour.get_location(slice, mb->mbAddrX, {xO - 1, yO + i});
        mbA [i] = slice->neighbour.get_mb      (slice, locA[i]);
        posA[i] = slice->neighbour.get_blkpos  (slice, locA[i]);
        mbA [i] = mbA[i] && mbA[i]->slice_nr == mb->slice_nr ? mbA[i] : nullptr;
    }
    loc_t locB = slice->neighbour.get_location(slice, mb->mbAddrX, {xO    , yO - 1});
    loc_t locD = slice->neighbour.get_location(slice, mb->mbAddrX, {xO - 1, yO - 1});
    mb_t* mbB  = slice->neighbour.get_mb      (slice, locB);
    mb_t* mbD  = slice->neighbour.get_mb      (slice, locD);
    pos_t posB = slice->neighbour.get_blkpos  (slice, locB);
    pos_t posD = slice->neighbour.get_blkpos  (slice, locD);
    mbB = mbB && mbB->slice_nr == mb->slice_nr ? mbB : nullptr;
    mbD = mbD && mbD->slice_nr == mb->slice_nr ? mbD : nullptr;

    if (pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < 16; ++i)
            available[0] &= mbA[i] && mbA[i]->is_intra_block ? 1 : 0;
        available[1] = mbB && mbB->is_intra_block ? 1 : 0;
        available[2] = 0;
        available[3] = mbD && mbD->is_intra_block ? 1 : 0;
    } else {
        available[0] = mbA[0] ? 1 : 0;
        available[1] = mbB    ? 1 : 0;
        available[2] = 0;
        available[3] = mbD    ? 1 : 0;
    }

#define p(x,y) (pred[((y) + 1) * 17 + ((x) + 1)])
    if (available[3]) {
        imgpel *pix = &img[posD.y][posD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < 16; ++y)
            p(-1, y) = img[posA[y].y][posA[y].x];
    }
    if (available[1]) {
        imgpel *pix = &img[posB.y][posB.x];
        for (int x = 0; x < 16; ++x)
            p(x, -1) = pix[x];
    }
#undef p
}

static void neighbouring_samples_chroma(imgpel* pred, bool* available, mb_t* mb, ColorPlane pl, int xO, int yO)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    imgpel** img = slice->dec_picture->imgUV[pl - 1];

    loc_t locA[16];
    mb_t* mbA [16];
    pos_t posA[16];
    for (int i = 0; i < sps->MbHeightC; ++i) {
        locA[i] = slice->neighbour.get_location_c(slice, mb->mbAddrX, {xO - 1, yO + i});
        mbA [i] = slice->neighbour.get_mb_c      (slice, locA[i]);
        posA[i] = slice->neighbour.get_blkpos_c  (slice, locA[i]);
        mbA [i] = mbA[i] && mbA[i]->slice_nr == mb->slice_nr ? mbA[i] : nullptr;
    }
    loc_t locB = slice->neighbour.get_location_c(slice, mb->mbAddrX, {xO    , yO - 1});
    loc_t locD = slice->neighbour.get_location_c(slice, mb->mbAddrX, {xO - 1, yO - 1});
    mb_t* mbB  = slice->neighbour.get_mb_c      (slice, locB);
    mb_t* mbD  = slice->neighbour.get_mb_c      (slice, locD);
    pos_t posB = slice->neighbour.get_blkpos_c  (slice, locB);
    pos_t posD = slice->neighbour.get_blkpos_c  (slice, locD);
    mbB = mbB && mbB->slice_nr == mb->slice_nr ? mbB : nullptr;
    mbD = mbD && mbD->slice_nr == mb->slice_nr ? mbD : nullptr;

    if (pps->constrained_intra_pred_flag) {
        available[0] = 1;
        for (int i = 0; i < sps->MbHeightC / 2; ++i)
            available[0] &= mbA[i] && mbA[i]->is_intra_block ? 1 : 0;
        available[1] = mbB && mbB->is_intra_block ? 1 : 0;
        available[2] = 1;
        for (int i = sps->MbHeightC / 2; i < sps->MbHeightC; ++i)
            available[2] &= mbA[i] && mbA[i]->is_intra_block ? 1 : 0;
        available[3] = mbD && mbD->is_intra_block ? 1 : 0;
    } else {
        available[0] = mbA[0] ? 1 : 0;
        available[1] = mbB ? 1 : 0;
        available[2] = mbA[0] ? 1 : 0;
        available[3] = mbD ? 1 : 0;
    }

#define p(x,y) (pred[((y) + 1) * 17 + ((x) + 1)])
    if (available[3]) {
        imgpel *pix = &img[posD.y][posD.x];
        p(-1, -1) = pix[0];
    }
    if (available[0]) {
        for (int y = 0; y < sps->MbHeightC / 2; ++y)
            p(-1, y) = img[posA[y].y][posA[y].x];
    }
    if (available[2]) {
        for (int y = sps->MbHeightC / 2; y < sps->MbHeightC; ++y)
            p(-1, y) = img[posA[y].y][posA[y].x];
    }
    if (available[1]) {
        imgpel *pix = &img[posB.y][posB.x];
        for (int x = 0; x < sps->MbWidthC; ++x)
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

#define Clip1(x) (clip3(0, (1 << BitDepth) - 1, x))
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
        int xO = ((chroma4x4BlkIdx / 4) % (size[0] / 8)) * 8 + ((chroma4x4BlkIdx % 4) % (size[0] / 4)) * 4;
        int yO = ((chroma4x4BlkIdx / 4) / (size[0] / 8)) * 8 + ((chroma4x4BlkIdx % 4) / (size[0] / 4)) * 4;

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

#define Clip1(x) (clip3(0, (1 << BitDepth) - 1, x))
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


void IntraPrediction::intra_pred_4x4(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    int BitDepth = sps->BitDepthY;
    imgpel* pred = &slice->mb_pred[pl][joff][ioff];
    imgpel pix[9 * 9];
    bool available[4];
    neighbouring_samples_4x4(pix, available, mb, pl, ioff, joff);

    int i4x4 = ((joff/4) / 2) * 8  + ((joff/4) % 2) * 2 + ((ioff/4) / 2) * 4 + ((ioff/4) % 2);
    uint8_t pred_mode = mb->Intra4x4PredMode[i4x4];

    switch (pred_mode) {
    case Intra_4x4_Vertical:
        if (!available[1])
            printf ("warning: Intra_4x4_Vertical prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_4x4_Horizontal:
        if (!available[0])
            printf ("warning: Intra_4x4_Horizontal prediction mode not allowed at mb %d\n",(int) slice->current_mb_nr);
        break;
    case Intra_4x4_DC:
        break;
    case Intra_4x4_Diagonal_Down_Left:
        if (!available[1])
            printf ("warning: Intra_4x4_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_4x4_Diagonal_Down_Right:
        if (!available[3])
            printf ("warning: Intra_4x4_Diagonal_Down_Right prediction mode not allowed at mb %d\n",(int) slice->current_mb_nr);
        break;
    case Intra_4x4_Vertical_Right:
        if (!available[3])
            printf ("warning: Intra_4x4_Vertical_Right prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_4x4_Horizontal_Down:  
        if (!available[3])
            printf ("warning: Intra_4x4_Horizontal_Down prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_4x4_Vertical_Left:
        if (!available[1])
            printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_4x4_Horizontal_Up:
        if (!available[0])
            printf ("warning: Intra_4x4_Horizontal_Up prediction mode not allowed at mb %d\n",(int) slice->current_mb_nr);
        break;
    default:
        printf("Error: illegal intra_4x4 prediction mode: %d\n", (int) pred_mode);
        break;
    }

    switch (pred_mode) {
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

void IntraPrediction::intra_pred_8x8(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int BitDepth = sps->BitDepthY;
    imgpel* pred = &slice->mb_pred[pl][joff][ioff];
    imgpel pix[17 * 17];
    bool available[4];
    neighbouring_samples_8x8(pix, available, mb, pl, ioff, joff);

    uint8_t pred_mode = mb->Intra8x8PredMode[joff/8 * 2 + ioff/8];

    switch (pred_mode) {
    case Intra_8x8_Vertical:
        if (!available[1])
            printf ("warning: Intra_8x8_Vertical prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_Horizontal:
        if (!available[0])
            printf ("warning: Intra_8x8_Horizontal prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_DC:
        break;
    case Intra_8x8_Diagonal_Down_Left:
        if (!available[1])
            printf ("warning: Intra_8x8_Diagonal_Down_Left prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_Diagonal_Down_Right:
        if (!available[3])
            printf ("warning: Intra_8x8_Diagonal_Down_Right prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_Vertical_Right:
        if (!available[3])
            printf ("warning: Intra_8x8_Vertical_Right prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_Horizontal_Down:  
        if (!available[3])
            printf ("warning: Intra_8x8_Horizontal_Down prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_Vertical_Left:
        if (!available[1])
            printf ("warning: Intra_4x4_Vertical_Left prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    case Intra_8x8_Horizontal_Up:
        if (!available[0])
            printf ("warning: Intra_8x8_Horizontal_Up prediction mode not allowed at mb %d\n", (int) slice->current_mb_nr);
        break;
    default:
        printf("Error: illegal intra_8x8 prediction mode: %d\n", (int) pred_mode);
        break;
    }

    switch (pred_mode) {
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

void IntraPrediction::intra_pred_16x16(mb_t* mb, ColorPlane pl, int ioff, int joff)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int BitDepth = sps->BitDepthY;
    imgpel* pred = &slice->mb_pred[pl][0][0];
    imgpel pix[17 * 17];
    bool available[4];
    neighbouring_samples_16x16(pix, available, mb, pl, 0, 0);

    switch (mb->Intra16x16PredMode) {
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
        printf("illegal 16x16 intra prediction mode input: %d\n", mb->Intra16x16PredMode);
        return;
    }

    switch (mb->Intra16x16PredMode) {
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

void IntraPrediction::intra_pred_chroma(mb_t* mb)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int BitDepth = sps->BitDepthC;
    imgpel* pred[2] = { &slice->mb_pred[1][0][0], &slice->mb_pred[2][0][0] };
    int mb_size[2] = { sps->MbWidthC, sps->MbHeightC };
    imgpel pix[2][17 * 17];
    bool available[2][4];
    neighbouring_samples_chroma(pix[0], available[0], mb, PLANE_U, 0, 0);
    neighbouring_samples_chroma(pix[1], available[1], mb, PLANE_V, 0, 0);

    switch (mb->intra_chroma_pred_mode) {
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

    switch (mb->intra_chroma_pred_mode) {
    case Intra_Chroma_DC:  
        intrapred_chroma_dc(pred, pix, available[0], mb_size, BitDepth);
        break;
    case Intra_Chroma_Horizontal: 
        intrapred_chroma_hor(pred, pix, available[0], mb_size);
        break;
    case Intra_Chroma_Vertical: 
        intrapred_chroma_ver(pred, pix, available[0], mb_size);
        break;
    case Intra_Chroma_Plane: 
        intrapred_chroma_plane(pred, pix, available[0], mb_size, BitDepth);
        break;
    }
}


}
}
