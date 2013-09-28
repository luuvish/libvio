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


namespace vio  {
namespace h264 {


loc_t operator + (const loc_t& l, const loc_t& r)
{
    return { l.x + r.x, l.y + r.y };
}

loc_t operator - (const loc_t& l, const loc_t& r)
{
    return { l.x - r.x, l.y - r.y };
}


loc_t Neighbour::get_location(slice_t* slice, bool chroma, int mbAddr, const pos_t& offset)
{
    sps_t* sps = slice->active_sps;

    int maxW = !chroma ? 16 : sps->MbWidthC;
    int maxH = !chroma ? 16 : sps->MbHeightC;

    loc_t loc {};

    if (slice->MbaffFrameFlag == 0) {
        loc.x = mbAddr % sps->PicWidthInMbs * maxW;
        loc.y = mbAddr / sps->PicWidthInMbs * maxH;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps->PicWidthInMbs * maxW;
        loc.y = (mbAddr / 2) / sps->PicWidthInMbs * maxH * 2;
        mb_t* mb = &this->mb_data[mbAddr];
        loc.x += offset.x;
        if (mb->mb_field_decoding_flag == 0) {
            loc.y += mbAddr % 2 * maxH;
            loc.y += offset.y;
        } else {
            loc.y += mbAddr % 2;
            loc.y += offset.y * 2;
        }
    }

    return loc;
}

mb_t* Neighbour::get_mb(slice_t* slice, bool chroma, int mbAddr, const pos_t& offset)
{
    sps_t* sps = slice->active_sps;

    int maxW = !chroma ? 16 : sps->MbWidthC;
    int maxH = !chroma ? 16 : sps->MbHeightC;

    loc_t loc {};

    if (slice->MbaffFrameFlag == 0) {
        loc.x = mbAddr % sps->PicWidthInMbs * maxW;
        loc.y = mbAddr / sps->PicWidthInMbs * maxH;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps->PicWidthInMbs * maxW;
        loc.y = (mbAddr / 2) / sps->PicWidthInMbs * maxH * 2;
        mb_t* mb = &this->mb_data[mbAddr];
        loc.x += offset.x;
        if (mb->mb_field_decoding_flag == 0) {
            loc.y += mbAddr % 2 * maxH;
            loc.y += offset.y;
        } else {
            loc.y += mbAddr % 2;
            loc.y += offset.y * 2;
        }
    }

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * maxW)
        return nullptr;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * maxH)
        return nullptr;

    mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps->PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps->PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    if (slice->MbaffFrameFlag)
        mb += ((mb->mb_field_decoding_flag == 0) ? (loc.y & maxH) : (loc.y & 1)) ? 1 : 0;
    return mb;
}

nb_t Neighbour::get_neighbour(slice_t* slice, bool chroma, int mbAddr, const pos_t& offset)
{
    sps_t* sps = slice->active_sps;

    int maxW = !chroma ? 16 : sps->MbWidthC;
    int maxH = !chroma ? 16 : sps->MbHeightC;

    loc_t loc {};

    if (slice->MbaffFrameFlag == 0) {
        loc.x = mbAddr % sps->PicWidthInMbs * maxW;
        loc.y = mbAddr / sps->PicWidthInMbs * maxH;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps->PicWidthInMbs * maxW;
        loc.y = (mbAddr / 2) / sps->PicWidthInMbs * maxH * 2;
        mb_t* mb = &this->mb_data[mbAddr];
        loc.x += offset.x;
        if (mb->mb_field_decoding_flag == 0) {
            loc.y += mbAddr % 2 * maxH;
            loc.y += offset.y;
        } else {
            loc.y += mbAddr % 2;
            loc.y += offset.y * 2;
        }
    }

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * maxW)
        return {nullptr, 0, 0};
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * maxH)
        return {nullptr, 0, 0};

    mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps->PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps->PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    pos_t pos {loc.x, loc.y};
    if (slice->MbaffFrameFlag) {
        if (mb->mb_field_decoding_flag == 0)
            mb += (loc.y & maxH) ? 1 : 0;
        else {
            mb += (loc.y & 1) ? 1 : 0;
            pos.y = loc.y / (maxH * 2) * (maxH * 2) + (loc.y % (maxH * 2)) / 2 + (loc.y & 1) * maxH;
        }
    }

    return {mb, pos.x, pos.y};
}

mb_t* Neighbour::get_mb(slice_t* slice, bool chroma, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    int maxW = !chroma ? 16 : sps->MbWidthC;
    int maxH = !chroma ? 16 : sps->MbHeightC;

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * maxW)
        return nullptr;
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * maxH)
        return nullptr;

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps->PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps->PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    if (slice->MbaffFrameFlag)
        mb += ((mb->mb_field_decoding_flag == 0) ? (loc.y & maxH) : (loc.y & 1)) ? 1 : 0;
    return mb;
}

nb_t Neighbour::get_neighbour(slice_t* slice, bool chroma, const loc_t& loc)
{
    sps_t* sps = slice->active_sps;

    int maxW = !chroma ? 16 : sps->MbWidthC;
    int maxH = !chroma ? 16 : sps->MbHeightC;

    if (loc.x < 0 || loc.x >= sps->PicWidthInMbs * maxW)
        return {nullptr, 0, 0};
    if (loc.y < 0 || loc.y >= slice->PicHeightInMbs * maxH)
        return {nullptr, 0, 0};

    int mbAddr = (slice->MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps->PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps->PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    pos_t pos {loc.x, loc.y};
    if (slice->MbaffFrameFlag) {
        if (mb->mb_field_decoding_flag == 0)
            mb += (loc.y & maxH) ? 1 : 0;
        else {
            mb += (loc.y & 1) ? 1 : 0;
            pos.y = loc.y / (maxH * 2) * (maxH * 2) + (loc.y % (maxH * 2)) / 2 + (loc.y & 1) * maxH;
        }
    }

    return {mb, pos.x, pos.y};
}


pos_t Neighbour::get_position(slice_t* slice, int mbAddr, int blkIdx)
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
        mb_t* mb = &this->mb_data[mbAddr];
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


int Neighbour::predict_nnz(mb_t* mb, int pl, int i, int j)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    pps_t* pps = slice->active_pps;

    bool chroma = !(pl == 0 || sps->separate_colour_plane_flag);

    nb_t nbA = slice->neighbour.get_neighbour(slice, chroma, mb->mbAddrX, {i - 1, j});
    nb_t nbB = slice->neighbour.get_neighbour(slice, chroma, mb->mbAddrX, {i, j - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb->slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb->slice_nr ? nbB.mb : nullptr;

    if (pps->constrained_intra_pred_flag && slice->parser.dp_mode == vio::h264::PAR_DP_3) {
        if (mb->is_intra_block) {
            nbA.mb = nbA.mb && nbA.mb->is_intra_block ? nbA.mb : nullptr;
            nbB.mb = nbB.mb && nbB.mb->is_intra_block ? nbB.mb : nullptr;
        }
    }

    int nW = !chroma ? 16 : sps->MbWidthC;
    int nH = !chroma ? 16 : sps->MbHeightC;

    uint8_t nA = 0;
    if (nbA.mb) {
        //if (nbA.mb->mb_type == PSKIP || nbA.mb->mb_type == BSKIP_DIRECT)
        //    nA = 0;
        //else if (nbA.mb->mb_type != IPCM && (nbA.mb->cbp & 15) == 0)
        //    nA = 0;
        //else if (nbA.mb->mb_type == IPCM)
        //    nA = 16;
        //else
            nA = nbA.mb->nz_coeff[pl][(nbA.y % nH) / 4][(nbA.x % nW) / 4];
    }

    uint8_t nB = 0;
    if (nbB.mb) {
        //if (nbB.mb->mb_type == PSKIP || nbB.mb->mb_type == BSKIP_DIRECT)
        //    nB = 0;
        //else if (nbB.mb->mb_type != IPCM && (nbB.mb->cbp & 15) == 0)
        //    nB = 0;
        //else if (nbB.mb->mb_type == IPCM)
        //    nB = 16;
        //else
            nB = nbB.mb->nz_coeff[pl][(nbB.y % nH) / 4][(nbB.x % nW) / 4];
    }

    uint8_t nC = nA + nB;
    if (nbA.mb && nbB.mb)
        nC = (nC + 1) >> 1;
    return nC;
}

    
}
}
