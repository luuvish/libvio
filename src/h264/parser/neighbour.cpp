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

#include "decoder.h"


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
    sps_t& sps = *slice->active_sps;
    shr_t& shr = slice->header;

    int maxW = !chroma ? 16 : sps.MbWidthC;
    int maxH = !chroma ? 16 : sps.MbHeightC;

    loc_t loc {};

    if (!shr.MbaffFrameFlag) {
        loc.x = mbAddr % sps.PicWidthInMbs * maxW;
        loc.y = mbAddr / sps.PicWidthInMbs * maxH;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps.PicWidthInMbs * maxW;
        loc.y = (mbAddr / 2) / sps.PicWidthInMbs * maxH * 2;
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
    sps_t& sps = *slice->active_sps;
    shr_t& shr = slice->header;

    int maxW = !chroma ? 16 : sps.MbWidthC;
    int maxH = !chroma ? 16 : sps.MbHeightC;

    loc_t loc {};

    if (!shr.MbaffFrameFlag) {
        loc.x = mbAddr % sps.PicWidthInMbs * maxW;
        loc.y = mbAddr / sps.PicWidthInMbs * maxH;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps.PicWidthInMbs * maxW;
        loc.y = (mbAddr / 2) / sps.PicWidthInMbs * maxH * 2;
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

    if (loc.x < 0 || loc.x >= sps.PicWidthInMbs * maxW)
        return nullptr;
    if (loc.y < 0 || loc.y >= shr.PicHeightInMbs * maxH)
        return nullptr;

    mbAddr = (shr.MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps.PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps.PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    if (shr.MbaffFrameFlag)
        mb += ((mb->mb_field_decoding_flag == 0) ? (loc.y & maxH) : (loc.y & 1)) ? 1 : 0;
    return mb;
}

nb_t Neighbour::get_neighbour(slice_t* slice, bool chroma, int mbAddr, const pos_t& offset)
{
    sps_t& sps = *slice->active_sps;
    shr_t& shr = slice->header;

    int maxW = !chroma ? 16 : sps.MbWidthC;
    int maxH = !chroma ? 16 : sps.MbHeightC;

    loc_t loc {};

    if (!shr.MbaffFrameFlag) {
        loc.x = mbAddr % sps.PicWidthInMbs * maxW;
        loc.y = mbAddr / sps.PicWidthInMbs * maxH;
        loc.x += offset.x;
        loc.y += offset.y;
    } else {
        loc.x = (mbAddr / 2) % sps.PicWidthInMbs * maxW;
        loc.y = (mbAddr / 2) / sps.PicWidthInMbs * maxH * 2;
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

    if (loc.x < 0 || loc.x >= sps.PicWidthInMbs * maxW)
        return {nullptr, 0, 0};
    if (loc.y < 0 || loc.y >= shr.PicHeightInMbs * maxH)
        return {nullptr, 0, 0};

    mbAddr = (shr.MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps.PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps.PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    pos_t pos {loc.x, loc.y};
    if (shr.MbaffFrameFlag) {
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
    sps_t& sps = *slice->active_sps;
    shr_t& shr = slice->header;

    int maxW = !chroma ? 16 : sps.MbWidthC;
    int maxH = !chroma ? 16 : sps.MbHeightC;

    if (loc.x < 0 || loc.x >= sps.PicWidthInMbs * maxW)
        return nullptr;
    if (loc.y < 0 || loc.y >= shr.PicHeightInMbs * maxH)
        return nullptr;

    int mbAddr = (shr.MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps.PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps.PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    if (shr.MbaffFrameFlag)
        mb += ((mb->mb_field_decoding_flag == 0) ? (loc.y & maxH) : (loc.y & 1)) ? 1 : 0;
    return mb;
}

nb_t Neighbour::get_neighbour(slice_t* slice, bool chroma, const loc_t& loc)
{
    sps_t& sps = *slice->active_sps;
    shr_t& shr = slice->header;

    int maxW = !chroma ? 16 : sps.MbWidthC;
    int maxH = !chroma ? 16 : sps.MbHeightC;

    if (loc.x < 0 || loc.x >= sps.PicWidthInMbs * maxW)
        return {nullptr, 0, 0};
    if (loc.y < 0 || loc.y >= shr.PicHeightInMbs * maxH)
        return {nullptr, 0, 0};

    int mbAddr = (shr.MbaffFrameFlag == 0) ?
        ((loc.y / maxH) * sps.PicWidthInMbs + (loc.x / maxW)) :
        ((loc.y / (maxH * 2)) * sps.PicWidthInMbs + (loc.x / maxW)) * 2;

    mb_t* mb = &this->mb_data[mbAddr];
    pos_t pos {loc.x, loc.y};
    if (shr.MbaffFrameFlag) {
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
    sps_t& sps = *slice->active_sps;
    shr_t& shr = slice->header;

    pos_t pos {};

    int blkX = ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int blkY = ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    if (!shr.MbaffFrameFlag) {
        pos.x = mbAddr % sps.PicWidthInMbs * 16;
        pos.y = mbAddr / sps.PicWidthInMbs * 16;
        pos.x += blkX * 4;
        pos.y += blkY * 4;
    } else {
        pos.x = (mbAddr / 2) % sps.PicWidthInMbs * 16;
        pos.y = (mbAddr / 2) / sps.PicWidthInMbs * 32;
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
    sps_t& sps = *slice->active_sps;
    pps_t& pps = *slice->active_pps;

    bool chroma = !(pl == 0 || sps.separate_colour_plane_flag);

    nb_t nbA = slice->neighbour.get_neighbour(slice, chroma, mb->mbAddrX, {i - 1, j});
    nb_t nbB = slice->neighbour.get_neighbour(slice, chroma, mb->mbAddrX, {i, j - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb->slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb->slice_nr ? nbB.mb : nullptr;

    if (pps.constrained_intra_pred_flag && slice->parser.dp_mode == vio::h264::PAR_DP_3) {
        if (mb->is_intra_block) {
            nbA.mb = nbA.mb && nbA.mb->is_intra_block ? nbA.mb : nullptr;
            nbB.mb = nbB.mb && nbB.mb->is_intra_block ? nbB.mb : nullptr;
        }
    }

    int nW = !chroma ? 16 : sps.MbWidthC;
    int nH = !chroma ? 16 : sps.MbHeightC;

    uint8_t nA = 0;
    if (nbA.mb) {
        //if (nbA.mb->mb_type == P_Skip || nbA.mb->mb_type == B_Skip)
        //    nA = 0;
        //else if (nbA.mb->mb_type != I_PCM && (nbA.mb->cbp & 15) == 0)
        //    nA = 0;
        //else if (nbA.mb->mb_type == I_PCM)
        //    nA = 16;
        //else
            nA = nbA.mb->nz_coeff[pl][(nbA.y % nH) / 4][(nbA.x % nW) / 4];
    }

    uint8_t nB = 0;
    if (nbB.mb) {
        //if (nbB.mb->mb_type == P_Skip || nbB.mb->mb_type == B_Skip)
        //    nB = 0;
        //else if (nbB.mb->mb_type != I_PCM && (nbB.mb->cbp & 15) == 0)
        //    nB = 0;
        //else if (nbB.mb->mb_type == I_PCM)
        //    nB = 16;
        //else
            nB = nbB.mb->nz_coeff[pl][(nbB.y % nH) / 4][(nbB.x % nW) / 4];
    }

    uint8_t nC = nA + nB;
    if (nbA.mb && nbB.mb)
        nC = (nC + 1) >> 1;
    return nC;
}



uint8_t intra_4x4_pred_mode(mb_t& mb, int bx, int by)
{
    slice_t& slice = *mb.p_Slice;
    pps_t& pps = *slice.active_pps;
    shr_t& shr = slice.header;

    nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx - 1, by});
    nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx, by - 1});

    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

    //get from array and decode
    if (pps.constrained_intra_pred_flag) {
        nbA.mb = nbA.mb && nbA.mb->is_intra_block ? nbA.mb : nullptr;
        nbB.mb = nbB.mb && nbB.mb->is_intra_block ? nbB.mb : nullptr;
    }
    // !! KS: not sure if the following is still correct...
    if (shr.slice_type == SI_slice) { // need support for MBINTLC1
        nbA.mb = nbA.mb && nbA.mb->mb_type == SI ? nbA.mb : nullptr;
        nbB.mb = nbB.mb && nbB.mb->mb_type == SI ? nbB.mb : nullptr;
    }

    bool dcPredModePredictedFlag = !nbA.mb || !nbB.mb;

    int scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };
    uint8_t intraMxMPredModeA = IntraPrediction::Intra_4x4_DC;
    uint8_t intraMxMPredModeB = IntraPrediction::Intra_4x4_DC;
    if (!dcPredModePredictedFlag) {
        if (nbA.mb->mb_type == I_8x8)
            intraMxMPredModeA = nbA.mb->Intra8x8PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4] / 4];
        else if (nbA.mb->mb_type == I_4x4)
            intraMxMPredModeA = nbA.mb->Intra4x4PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4]];
        if (nbB.mb->mb_type == I_8x8)
            intraMxMPredModeB = nbB.mb->Intra8x8PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4] / 4];
        else if (nbB.mb->mb_type == I_4x4)
            intraMxMPredModeB = nbB.mb->Intra4x4PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4]];
    }

    uint8_t predIntra4x4PredMode = min(intraMxMPredModeA, intraMxMPredModeB);

    return predIntra4x4PredMode;
}

uint8_t intra_8x8_pred_mode(mb_t& mb, int bx, int by)
{
    slice_t& slice = *mb.p_Slice;
    pps_t& pps = *slice.active_pps;

    nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx - 1, by});
    nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {bx, by - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

    //get from array and decode
    if (pps.constrained_intra_pred_flag) {
        nbA.mb = nbA.mb && nbA.mb->is_intra_block ? nbA.mb : nullptr;
        nbB.mb = nbB.mb && nbB.mb->is_intra_block ? nbB.mb : nullptr;
    }

    bool dcPredModePredictedFlag = !nbA.mb || !nbB.mb;

    int scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };
    uint8_t intraMxMPredModeA = IntraPrediction::Intra_8x8_DC;
    uint8_t intraMxMPredModeB = IntraPrediction::Intra_8x8_DC;
    if (!dcPredModePredictedFlag) {
        if (nbA.mb->mb_type == I_8x8)
            intraMxMPredModeA = nbA.mb->Intra8x8PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4] / 4];
        else if (nbA.mb->mb_type == I_4x4)
            intraMxMPredModeA = nbA.mb->Intra4x4PredMode[scan[(nbA.y & 12) + (nbA.x & 15) / 4]];
        if (nbB.mb->mb_type == I_8x8)
            intraMxMPredModeB = nbB.mb->Intra8x8PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4] / 4];
        else if (nbB.mb->mb_type == I_4x4)
            intraMxMPredModeB = nbB.mb->Intra4x4PredMode[scan[(nbB.y & 12) + (nbB.x & 15) / 4]];
    }

    uint8_t predIntra8x8PredMode = min(intraMxMPredModeA, intraMxMPredModeB);

    return predIntra8x8PredMode;
}



CtxIdxInc::CtxIdxInc(mb_t& _mb) :
    sps { *_mb.p_Slice->active_sps },
    pps { *_mb.p_Slice->active_pps },
    slice { *_mb.p_Slice },
    mb { _mb }
{
    this->mb_data = _mb.p_Slice->neighbour.mb_data;
}

CtxIdxInc::~CtxIdxInc()
{
}


int CtxIdxInc::mb_skip_flag()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && !mbA->mb_skip_flag ? 1 : 0;
    int condTermFlagB = mbB && !mbB->mb_skip_flag ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::mb_field_decoding_flag()
{
    shr_t& shr = slice.header;

    int topMbAddr = shr.MbaffFrameFlag ? mb.mbAddrX & ~1 : mb.mbAddrX;

    mb_t* mbA = this->get_mb(&slice, false, topMbAddr, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, topMbAddr, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->mb_field_decoding_flag ? 1 : 0;
    int condTermFlagB = mbB && mbB->mb_field_decoding_flag ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::mb_type_si_slice()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->mb_type != SI ? 1 : 0;
    int condTermFlagB = mbB && mbB->mb_type != SI ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::mb_type_i_slice()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->mb_type != I_4x4 && mbA->mb_type != I_8x8 ? 1 : 0;
    int condTermFlagB = mbB && mbB->mb_type != I_4x4 && mbB->mb_type != I_8x8 ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::mb_type_b_slice()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->mb_type != 0 ? 1 : 0;
    int condTermFlagB = mbB && mbB->mb_type != 0 ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::transform_size_8x8_flag()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->transform_size_8x8_flag ? 1 : 0;
    int condTermFlagB = mbB && mbB->transform_size_8x8_flag ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::intra_chroma_pred_mode()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && mbA->intra_chroma_pred_mode != 0 && mbA->mb_type != I_PCM ? 1 : 0;
    int condTermFlagB = mbB && mbB->intra_chroma_pred_mode != 0 && mbB->mb_type != I_PCM ? 1 : 0;
    int ctxIdxInc = condTermFlagA + condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::ref_idx_l(uint8_t list, uint8_t x0, uint8_t y0)
{
    shr_t& shr = slice.header;

    nb_t nbA = this->get_neighbour(&slice, false, mb.mbAddrX, {x0 * 4 - 1, y0 * 4});
    nb_t nbB = this->get_neighbour(&slice, false, mb.mbAddrX, {x0 * 4, y0 * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

    int condTermFlagA = 0;
    int condTermFlagB = 0;
    int ctxIdxInc;

#define IS_DIRECT(MB) ((MB)->mb_type == 0 && (shr.slice_type == B_slice))
    if (nbA.mb) {
        int ref_idx_lX = slice.dec_picture->mv_info[nbA.y / 4][nbA.x / 4].ref_idx[list];
        int refIdxZeroFlagA = 0;
        if (shr.MbaffFrameFlag && !mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag)
            refIdxZeroFlagA = (ref_idx_lX > 1 ? 0 : 1);
        else
            refIdxZeroFlagA = (ref_idx_lX > 0 ? 0 : 1);

        int mbPartIdxA = ((nbA.y / 4) & 2) + ((nbA.x / 8) & 1);
        int predModeEqualFlagA = 1;
        if (IS_DIRECT(nbA.mb) || nbA.mb->SubMbType[mbPartIdxA] == 0)
            predModeEqualFlagA = 0;
        //else if (nbA.mb->SubMbPredMode[mbPartIdxA] != Pred_LX | BiPred)
        //    predModeEqualFlagA = 0;
        condTermFlagA = (nbA.mb->mb_type == P_Skip || nbA.mb->mb_type == B_Skip ||
                         nbA.mb->is_intra_block ||
                         predModeEqualFlagA == 0 || refIdxZeroFlagA == 1) ? 0 : 1;
    }
    if (nbB.mb) {
        int ref_idx_lX = slice.dec_picture->mv_info[nbB.y / 4][nbB.x / 4].ref_idx[list];
        int refIdxZeroFlagB = 0;
        if (shr.MbaffFrameFlag && !mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag)
            refIdxZeroFlagB = (ref_idx_lX > 1 ? 0 : 1);
        else
            refIdxZeroFlagB = (ref_idx_lX > 0 ? 0 : 1);

        int mbPartIdxB = ((nbB.y / 4) & 2) + ((nbB.x / 8) & 1);
        int predModeEqualFlagB = 1;
        if (IS_DIRECT(nbB.mb) || nbB.mb->SubMbType[mbPartIdxB] == 0)
            predModeEqualFlagB = 0;
        //else if (nbB.mb->SubMbPredMode[mbPartIdxB] != Pred_LX | BiPred)
        //    predModeEqualFlagB = 0;
        condTermFlagB = (nbB.mb->mb_type == P_Skip || nbB.mb->mb_type == B_Skip ||
                         nbB.mb->is_intra_block ||
                         predModeEqualFlagB == 0 || refIdxZeroFlagB == 1) ? 0 : 1;
    }
#undef IS_DIRECT

    ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::mvd_l(uint8_t list, uint8_t x0, uint8_t y0, bool compIdx)
{
    shr_t& shr = slice.header;

    nb_t nbA = this->get_neighbour(&slice, false, mb.mbAddrX, {x0 * 4 - 1, y0 * 4});
    nb_t nbB = this->get_neighbour(&slice, false, mb.mbAddrX, {x0 * 4, y0 * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

    int absMvdCompA = 0;
    int absMvdCompB = 0;

#define IS_DIRECT(MB) ((MB)->mb_type == 0 && (shr.slice_type == B_slice))
    if (nbA.mb) {
        int mbPartIdxA = ((nbA.y / 4) & 2) + ((nbA.x / 8) & 1);
        int predModeEqualFlagA = 1;
        if (IS_DIRECT(nbA.mb) || nbA.mb->SubMbType[mbPartIdxA] == 0)
            predModeEqualFlagA = 0;
        //else if (nbA.mb->SubMbPredMode[mbPartIdxA] != Pred_LX | BiPred)
        //    predModeEqualFlagA = 0;

        if (!(nbA.mb->mb_type == P_Skip || nbA.mb->mb_type == B_Skip ||
              nbA.mb->is_intra_block || predModeEqualFlagA == 0)) {
            auto mvd_lX = (list == 0) ? nbA.mb->mvd_l0 : nbA.mb->mvd_l1;
            absMvdCompA = abs(mvd_lX[(nbA.y & 15) / 4][(nbA.x & 15) / 4][compIdx]);
            if (shr.MbaffFrameFlag && compIdx) {
                if (!mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag)
                    absMvdCompA *= 2;
                else if (mb.mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag)
                    absMvdCompA /= 2;
            }
        }
    }
    if (nbB.mb) {
        int mbPartIdxB = ((nbB.y / 4) & 2) + ((nbB.x / 8) & 1);
        int predModeEqualFlagB = 1;
        if (IS_DIRECT(nbB.mb) || nbB.mb->SubMbType[mbPartIdxB] == 0)
            predModeEqualFlagB = 0;
        //else if (nbB.mb->SubMbPredMode[mbPartIdxB] != Pred_LX | BiPred)
        //    predModeEqualFlagB = 0;

        if (!(nbB.mb->mb_type == P_Skip || nbB.mb->mb_type == B_Skip ||
              nbB.mb->is_intra_block || predModeEqualFlagB == 0)) {
            auto mvd_lX = (list == 0) ? nbB.mb->mvd_l0 : nbB.mb->mvd_l1;
            absMvdCompB = abs(mvd_lX[(nbB.y & 15) / 4][(nbB.x & 15) / 4][compIdx]);
            if (shr.MbaffFrameFlag && compIdx) {
                if (!mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag)
                    absMvdCompB *= 2;
                else if (mb.mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag)
                    absMvdCompB /= 2;
            }
        }
    }
#undef IS_DIRECT

    int absMvdSum = absMvdCompA + absMvdCompB;
    int ctxIdxInc = (absMvdSum < 3 ? 0 : absMvdSum <= 32 ? 1 : 2);

    return ctxIdxInc;
}

int CtxIdxInc::coded_block_pattern_luma(uint8_t x0, uint8_t y0, uint8_t coded_block_pattern)
{
    nb_t nbA = this->get_neighbour(&slice, false, mb.mbAddrX, {x0 * 4 - 1, y0 * 4});
    nb_t nbB = this->get_neighbour(&slice, false, mb.mbAddrX, {x0 * 4, y0 * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

    int cbp_a = 0x3F, cbp_b = 0x3F;
    int cbp_a_idx = 0, cbp_b_idx = 0;
    if (x0 == 0) {
        if (nbA.mb && nbA.mb->mb_type != I_PCM) {
            cbp_a = nbA.mb->CodedBlockPatternLuma;
            cbp_a_idx = (((nbA.y & 15) / 4) & ~1) + 1;
        }
    } else {
        cbp_a = coded_block_pattern;
        cbp_a_idx = y0;
    }
    if (y0 == 0) {
        if (nbB.mb && nbB.mb->mb_type != I_PCM) {
            cbp_b = nbB.mb->CodedBlockPatternLuma;
            cbp_b_idx = (x0 / 2) + 2;
        }
    } else {
        cbp_b = coded_block_pattern;
        cbp_b_idx = (x0 / 2);
    }

    int condTermFlagA = (cbp_a & (1 << cbp_a_idx)) == 0 ? 1 : 0;
    int condTermFlagB = (cbp_b & (1 << cbp_b_idx)) == 0 ? 1 : 0;
    int ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    return ctxIdxInc;
}

int CtxIdxInc::coded_block_pattern_chroma()
{
    mb_t* mbA = this->get_mb(&slice, false, mb.mbAddrX, {-1, 0});
    mb_t* mbB = this->get_mb(&slice, false, mb.mbAddrX, {0, -1});
    mbA = mbA && mbA->slice_nr == mb.slice_nr ? mbA : nullptr;
    mbB = mbB && mbB->slice_nr == mb.slice_nr ? mbB : nullptr;

    int condTermFlagA = mbA && (mbA->mb_type == I_PCM || mbA->CodedBlockPatternChroma) ? 1 : 0;
    int condTermFlagB = mbB && (mbB->mb_type == I_PCM || mbB->CodedBlockPatternChroma) ? 1 : 0;
    int ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    condTermFlagA = mbA && (mbA->mb_type == I_PCM || mbA->CodedBlockPatternChroma == 2) ? 1 : 0;
    condTermFlagB = mbB && (mbB->mb_type == I_PCM || mbB->CodedBlockPatternChroma == 2) ? 1 : 0;
    ctxIdxInc |= (condTermFlagA + 2 * condTermFlagB + 4) << 2;

    return ctxIdxInc;
}


//int CtxIdxInc::coded_block_flag(int pl, bool chroma, bool ac, int blkIdx)
int coded_block_flag_ctxIdxInc(mb_t& mb, int pl, bool chroma, bool ac, int blkIdx)
{
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;

    int i = chroma ? blkIdx % 2 : ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int j = chroma ? blkIdx / 2 : ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    int y_dc = (!chroma && !ac); 
    int u_dc = (chroma && !ac && pl == 1);
    int v_dc = (chroma && !ac && pl == 2);
    int y_ac = (!chroma && ac);
    int u_ac = (chroma && ac && pl == 1);
    //int v_ac = (chroma && ac && pl == 2);

    int temp_pl = (sps.ChromaArrayType == 3 ? pl : 0);
    int bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35);
    int bit_pos_a = 0;
    int bit_pos_b = 0;

    nb_t nbA = slice.neighbour.get_neighbour(&slice, chroma, mb.mbAddrX, {i * 4 - 1, j * 4});
    nb_t nbB = slice.neighbour.get_neighbour(&slice, chroma, mb.mbAddrX, {i * 4, j * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;

    int nW = !chroma ? 16 : sps.MbWidthC;
    int nH = !chroma ? 16 : sps.MbHeightC;
    if (ac) {
        if (nbA.mb)
            bit_pos_a = ((nbA.y % nH) & 12) + (nbA.x % nW) / 4;
        if (nbB.mb)
            bit_pos_b = ((nbB.y % nH) & 12) + (nbB.x % nW) / 4;
    }

    int condTermFlagA = (mb.is_intra_block ? 1 : 0);
    int condTermFlagB = (mb.is_intra_block ? 1 : 0);

    if (nbA.mb) {
        if (nbA.mb->mb_type == I_PCM)
            condTermFlagA = 1;
        else
            condTermFlagA = (nbA.mb->cbp_bits[temp_pl] >> (bit + bit_pos_a)) & 1;
    }
    if (nbB.mb) {
        if (nbB.mb->mb_type == I_PCM)
            condTermFlagB = 1;
        else
            condTermFlagB = (nbB.mb->cbp_bits[temp_pl] >> (bit + bit_pos_b)) & 1;
    }
    int ctxIdxInc = condTermFlagA + 2 * condTermFlagB;

    return ctxIdxInc;
}

void update_coded_block_flag(mb_t* mb, int pl, bool chroma, bool ac, int blkIdx)
{
    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    int i = chroma ? blkIdx % 2 : ((blkIdx / 4) % 2) * 2 + (blkIdx % 4) % 2;
    int j = chroma ? blkIdx / 2 : ((blkIdx / 4) / 2) * 2 + (blkIdx % 4) / 2;

    int y_dc = (!chroma && !ac);
    int u_dc = (chroma && !ac && pl == 1);
    int v_dc = (chroma && !ac && pl == 2);
    int y_ac = (!chroma && ac);
    int u_ac = (chroma && ac && pl == 1);
    //int v_ac = (chroma && ac && pl == 2);

    int temp_pl = (sps->ChromaArrayType == 3 ? pl : 0);
    int cbp = (mb->transform_size_8x8_flag && !chroma && ac ? 0x33 : 0x01);
    int bit = (y_dc ? 0 : y_ac ? 1 : u_dc ? 17 : v_dc ? 18 : u_ac ? 19 : 35) + (ac ? j * 4 + i : 0);

    mb->cbp_bits[temp_pl] |= ((uint64_t)cbp << bit);
}

    
}
}
