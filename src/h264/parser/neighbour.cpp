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
 *  File      : neighbour.cpp
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
#include "neighbour.h"


loc_t operator + (const loc_t& l, const loc_t& r)
{
    return { l.x + r.x, l.y + r.y };
}

loc_t operator - (const loc_t& l, const loc_t& r)
{
    return { l.x - r.x, l.y - r.y };
}


loc_t neighbour_t::get_location(slice_t* slice, int mbAddr, const pos_t& offset)
{
    sps_t* sps = slice->active_sps;

    loc_t loc {};

    if (slice->MbaffFrameFlag == 0) {
        loc.x = mbAddr % sps->PicWidthInMbs * 16;
        loc.y = mbAddr / sps->PicWidthInMbs * 16;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps->PicWidthInMbs * 16;
        loc.y = (mbAddr / 2) / sps->PicWidthInMbs * 32;
        mb_t* mb = &slice->mb_data[mbAddr];
        loc.x += offset.x;
        if (mb->mb_field_decoding_flag == 0) {
            loc.y += mbAddr % 2 * 16;
            loc.y += offset.y;
        } else {
            loc.y += mbAddr % 2;
            loc.y += offset.y * 2;
        }
    }

    return loc;
}

loc_t neighbour_t::get_location_c(slice_t* slice, int mbAddr, const pos_t& offset)
{
    sps_t* sps = slice->active_sps;

    loc_t loc {};

    if (slice->MbaffFrameFlag == 0) {
        loc.x = mbAddr % sps->PicWidthInMbs * sps->MbWidthC;
        loc.y = mbAddr / sps->PicWidthInMbs * sps->MbHeightC;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps->PicWidthInMbs * sps->MbWidthC;
        loc.y = (mbAddr / 2) / sps->PicWidthInMbs * sps->MbHeightC * 2;
        mb_t* mb = &slice->mb_data[mbAddr];
        loc.x += offset.x;
        if (mb->mb_field_decoding_flag == 0) {
            loc.y += mbAddr % 2 * sps->MbHeightC;
            loc.y += offset.y;
        } else {
            loc.y += mbAddr % 2;
            loc.y += offset.y * 2;
        }
    }

    return loc;
}

mb_t* neighbour_t::get_mb(slice_t* slice, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * 16)
        return nullptr;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * 16)
        return nullptr;

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / 16) * sps->PicWidthInMbs + (loc.x / 16)) :
        ((loc.y / 32) * sps->PicWidthInMbs + (loc.x / 16)) * 2;

    mb_t* mb = &slice->mb_data[mbAddr];
    if (slice->MbaffFrameFlag)
        mb += ((mb->mb_field_decoding_flag == 0) ? (loc.y & 16) : (loc.y & 1)) ? 1 : 0;
    return mb;
}

mb_t* neighbour_t::get_mb_c(slice_t* slice, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * sps->MbWidthC)
        return nullptr;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * sps->MbHeightC)
        return nullptr;

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / sps->MbHeightC      ) * sps->PicWidthInMbs + (loc.x / sps->MbWidthC)) :
        ((loc.y / (sps->MbHeightC * 2)) * sps->PicWidthInMbs + (loc.x / sps->MbWidthC)) * 2;

    mb_t* mb = &slice->mb_data[mbAddr];
    if (slice->MbaffFrameFlag)
        mb += ((mb->mb_field_decoding_flag == 0) ? (loc.y & sps->MbHeightC) : (loc.y & 1)) ? 1 : 0;
    return mb;
}

int neighbour_t::get_mbaddr(slice_t* slice, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * 16)
        return -1;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * 16)
        return -1;

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / 16) * sps->PicWidthInMbs + (loc.x / 16)) :
        ((loc.y / 32) * sps->PicWidthInMbs + (loc.x / 16)) * 2;

    if (slice->MbaffFrameFlag) {
        mb_t* mb = &slice->mb_data[mbAddr];
        mbAddr += ((mb->mb_field_decoding_flag == 0) ? (loc.y & 16) : (loc.y & 1)) ? 1 : 0;
    }
    return mbAddr;
}

pos_t neighbour_t::get_blkpos(slice_t* slice, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    pos_t pos {};

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * 16)
        return pos;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * 16)
        return pos;

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / 16) * sps->PicWidthInMbs + (loc.x / 16)) :
        ((loc.y / 32) * sps->PicWidthInMbs + (loc.x / 16)) * 2;

    pos.x = loc.x;
    pos.y = loc.y;
    if (slice->MbaffFrameFlag) {
        mb_t* mb = &slice->mb_data[mbAddr];
        if (mb->mb_field_decoding_flag)
            pos.y = loc.y / 32 * 32 + (loc.y % 32) / 2 + 16 * (loc.y & 1);
    }
    return pos;
}

pos_t neighbour_t::get_blkpos_c(slice_t* slice, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    pos_t pos {};

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * sps->MbWidthC)
        return pos;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * sps->MbHeightC)
        return pos;

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / sps->MbHeightC      ) * sps->PicWidthInMbs + (loc.x / sps->MbWidthC)) :
        ((loc.y / (sps->MbHeightC * 2)) * sps->PicWidthInMbs + (loc.x / sps->MbWidthC)) * 2;

    pos.x = loc.x;
    pos.y = loc.y;
    if (slice->MbaffFrameFlag) {
        mb_t* mb = &slice->mb_data[mbAddr];
        if (mb->mb_field_decoding_flag)
            pos.y = loc.y / (sps->MbHeightC * 2) * (sps->MbHeightC * 2) +
                    (loc.y % (sps->MbHeightC * 2)) / 2 + sps->MbHeightC * (loc.y & 1);
    }
    return pos;
}

pos_t neighbour_t::get_position(slice_t* slice, int mbAddr, int blkIdx)
{
    sps_t* sps = slice->active_sps;

    pos_t pos {};

    int blkX = ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int blkY = ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    if (slice->MbaffFrameFlag == 0) {
        pos.x = mbAddr % sps->PicWidthInMbs * 16;
        pos.y = mbAddr / sps->PicWidthInMbs * 16;
        pos.x += blkX * 4;
        pos.y += blkY * 4;
    } else {
        pos.x = (mbAddr / 2) % sps->PicWidthInMbs * 16;
        pos.y = (mbAddr / 2) / sps->PicWidthInMbs * 32;
        mb_t* mb = &slice->mb_data[mbAddr];
        pos.x += blkX * 4;
        if (mb->mb_field_decoding_flag == 0) {
            pos.y += mbAddr % 2 * 16;
            pos.y += blkY * 4;
        } else {
            pos.y += mbAddr % 2;
            pos.y += blkY * 8;
        }
    }

    return pos;
}


int neighbour_t::predict_nnz(mb_t* mb, int pl, int i, int j)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    loc_t locA, locB;
    mb_t* mbA, *mbB;
    pos_t posA, posB;
    if (pl == 0 || sps->separate_colour_plane_flag) {
        locA = slice->neighbour.get_location(slice, mb->mbAddrX, {i - 1, j});
        locB = slice->neighbour.get_location(slice, mb->mbAddrX, {i, j - 1});
        mbA  = slice->neighbour.get_mb      (slice, locA);
        mbB  = slice->neighbour.get_mb      (slice, locB);
        posA = slice->neighbour.get_blkpos  (slice, locA);
        posB = slice->neighbour.get_blkpos  (slice, locB);
    } else {
        locA = slice->neighbour.get_location_c(slice, mb->mbAddrX, {i - 1, j});
        locB = slice->neighbour.get_location_c(slice, mb->mbAddrX, {i, j - 1});
        mbA  = slice->neighbour.get_mb_c      (slice, locA);
        mbB  = slice->neighbour.get_mb_c      (slice, locB);
        posA = slice->neighbour.get_blkpos_c  (slice, locA);
        posB = slice->neighbour.get_blkpos_c  (slice, locB);
    }

    mbA = mbA && mbA->slice_nr == mb->slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb->slice_nr ? mbB : nullptr;

    if (pps->constrained_intra_pred_flag && slice->dp_mode == vio::h264::PAR_DP_3) {
        if (mb->is_intra_block) {
            mbA = mbA && mbA->is_intra_block ? mbA : nullptr;
            mbB = mbB && mbB->is_intra_block ? mbB : nullptr;
        }
    }

    int nW = pl == 0 || sps->separate_colour_plane_flag ? 16 : sps->MbWidthC;
    int nH = pl == 0 || sps->separate_colour_plane_flag ? 16 : sps->MbHeightC;

    uint8_t nA = 0;
    if (mbA) {
        //if (mbA->mb_type == PSKIP || mbA->mb_type == BSKIP_DIRECT)
        //    nA = 0;
        //else if (mbA->mb_type != IPCM && (mbA->cbp & 15) == 0)
        //    nA = 0;
        //else if (mbA->mb_type == IPCM)
        //    nA = 16;
        //else
            nA = mbA->nz_coeff[pl][(posA.y % nH) / 4][(posA.x % nW) / 4];
    }

    uint8_t nB = 0;
    if (mbB) {
        //if (mbB->mb_type == PSKIP || mbB->mb_type == BSKIP_DIRECT)
        //    nB = 0;
        //else if (mbB->mb_type != IPCM && (mbB->cbp & 15) == 0)
        //    nB = 0;
        //else if (mbB->mb_type == IPCM)
        //    nB = 16;
        //else
            nB = mbB->nz_coeff[pl][(posB.y % nH) / 4][(posB.x % nW) / 4];
    }

    uint8_t nC = nA + nB;
    if (mbA && mbB)
        nC = (nC + 1) >> 1;
    return nC;
}


void get_neighbors(mb_t* mb, PixelPos* block, int mb_x, int mb_y, int blockshape_x)
{
    slice_t* slice = mb->p_Slice;

    loc_t locA = slice->neighbour.get_location(slice, mb->mbAddrX, {mb_x - 1           , mb_y    });
    loc_t locB = slice->neighbour.get_location(slice, mb->mbAddrX, {mb_x               , mb_y - 1});
    loc_t locC = slice->neighbour.get_location(slice, mb->mbAddrX, {mb_x + blockshape_x, mb_y - 1});
    loc_t locD = slice->neighbour.get_location(slice, mb->mbAddrX, {mb_x - 1           , mb_y - 1});
    mb_t* mbA  = slice->neighbour.get_mb      (slice, locA);
    mb_t* mbB  = slice->neighbour.get_mb      (slice, locB);
    mb_t* mbC  = slice->neighbour.get_mb      (slice, locC);
    mb_t* mbD  = slice->neighbour.get_mb      (slice, locD);
    pos_t posA = slice->neighbour.get_blkpos  (slice, locA);
    pos_t posB = slice->neighbour.get_blkpos  (slice, locB);
    pos_t posC = slice->neighbour.get_blkpos  (slice, locC);
    pos_t posD = slice->neighbour.get_blkpos  (slice, locD);

    mbA = mbA && mbA->slice_nr == mb->slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb->slice_nr ? mbB : nullptr;
    mbC = mbC && mbC->slice_nr == mb->slice_nr ? mbC : nullptr;
    mbD = mbD && mbD->slice_nr == mb->slice_nr ? mbD : nullptr;

    if (mb_y > 0) {
        if (mb_x < 8) { // first column of 8x8 blocks
            if (mb_y == 8) {
                if (blockshape_x == MB_BLOCK_SIZE)
                    mbC = nullptr;
            } else if (mb_x + blockshape_x == 8)
                mbC = nullptr;
        } else if (mb_x + blockshape_x == MB_BLOCK_SIZE)
            mbC = nullptr;
    }

    if (!mbC) {
        mbC  = mbD;
        posC = posD;
    }

    block[0].available = mbA ? 1 : 0;
    block[0].mb_addr   = mbA ? mbA->mbAddrX : 0;
    block[0].x     = (posA.x & 15) / 4;
    block[0].y     = (posA.y & 15) / 4;
    block[0].pos_x = (posA.x     ) / 4;
    block[0].pos_y = (posA.y     ) / 4;

    block[1].available = mbB ? 1 : 0;
    block[1].mb_addr   = mbB ? mbB->mbAddrX : 0;
    block[1].x     = (posB.x & 15) / 4;
    block[1].y     = (posB.y & 15) / 4;
    block[1].pos_x = (posB.x     ) / 4;
    block[1].pos_y = (posB.y     ) / 4;

    block[2].available = mbC ? 1 : 0;
    block[2].mb_addr   = mbC ? mbC->mbAddrX : 0;
    block[2].x     = (posC.x & 15) / 4;
    block[2].y     = (posC.y & 15) / 4;
    block[2].pos_x = (posC.x     ) / 4;
    block[2].pos_y = (posC.y     ) / 4;
}
