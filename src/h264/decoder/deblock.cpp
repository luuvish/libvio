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
 *  File      : deblock.cpp
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
#include "image.h"
#include "neighbour.h"
#include "deblock.h"


namespace vio  {
namespace h264 {


deblock_t deblock;


static inline int compare_mvs(const MotionVector *mv0, const MotionVector *mv1, int mvlimit)
{
    return (abs(mv0->mv_x - mv1->mv_x) >= 4) | (abs(mv0->mv_y - mv1->mv_y) >= mvlimit);
}

static inline int bs_compare_mvs(int blk_x1, int blk_y1, int blk_x2, int blk_y2, int mvlimit, storable_picture* p)
{
    int StrValue;

    pic_motion_params* mv_info_p = &p->mv_info[blk_y1][blk_x1];
    pic_motion_params* mv_info_q = &p->mv_info[blk_y2][blk_x2];
    storable_picture*  ref_p0    = mv_info_p->ref_pic[LIST_0];
    storable_picture*  ref_q0    = mv_info_q->ref_pic[LIST_0];
    storable_picture*  ref_p1    = mv_info_p->ref_pic[LIST_1];
    storable_picture*  ref_q1    = mv_info_q->ref_pic[LIST_1];

    if ((ref_p0 == ref_q0 && ref_p1 == ref_q1) || (ref_p0 == ref_q1 && ref_p1 == ref_q0)) {
        // L0 and L1 reference pictures of p0 are different; q0 as well
        if (ref_p0 != ref_p1) {
            // compare MV for the same reference picture
            if (ref_p0 == ref_q0) {
                StrValue = 
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit);
            } else {
                StrValue =
                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit);
            }
        } else { // L0 and L1 reference pictures of p0 are the same; q0 as well
            StrValue = ((
                compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                & (
                compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)
            ));
        }
    } else
        StrValue = 1;

    return StrValue;
}


void deblock_t::get_strength_ver(mb_t* MbQ, int edge)
{
    uint8_t* Strength = MbQ->strength_ver[edge];
    int StrValue;

    slice_t* slice = MbQ->p_Slice;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

    PixelPos pixP, pixQ;
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
    if (slice->MbaffFrameFlag) {
        getAffNeighbour(MbQ, edge * 4 - 1, 0, mb_size, &pixP);
        getAffNeighbour(MbQ, edge * 4    , 0, mb_size, &pixQ);
    } else {
        getNonAffNeighbour(MbQ, edge * 4 - 1, 0, mb_size, &pixP);
        getNonAffNeighbour(MbQ, edge * 4    , 0, mb_size, &pixQ);
    }
    mb_t* MbP = &MbQ->p_Vid->mb_data[pixP.mb_addr];
    MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag);

    bool verticalEdgeFlag  = 1;
    bool mixedModeEdgeFlag = slice->MbaffFrameFlag &&
                             MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag;

    bool special = MbP->p_Slice->slice_type == SI_slice || MbP->p_Slice->slice_type == SP_slice ||
                   MbQ->p_Slice->slice_type == SI_slice || MbQ->p_Slice->slice_type == SP_slice;
    bool field = slice->MbaffFrameFlag ?
        MbP->mb_field_decoding_flag || MbQ->mb_field_decoding_flag : slice->field_pic_flag;

    bool cond_bS4 = !field ||
                    ((slice->MbaffFrameFlag || slice->field_pic_flag) && verticalEdgeFlag);
    bool cond_bS3 = !mixedModeEdgeFlag || !verticalEdgeFlag;

    if (edge == 0 && cond_bS4 && special) {
        StrValue = 4;
        memset(Strength, StrValue, MB_BLOCK_SIZE * sizeof(uint8_t));
        return;
    }
    if (cond_bS3 && special) {
        StrValue = 3;
        memset(Strength, StrValue, MB_BLOCK_SIZE * sizeof(uint8_t));
        return;
    }

    if (edge > 0 && slice->slice_type == P_slice && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, MB_BLOCK_SIZE * sizeof(uint8_t));
        return;
    }

    for (int idx = 0; idx < MB_BLOCK_SIZE; ++idx) {
        if (slice->MbaffFrameFlag)
            getAffNeighbour(MbQ, edge * 4 - 1, idx, mb_size, &pixP);
        else
            getNonAffNeighbour(MbQ, edge * 4 - 1, idx, mb_size, &pixP);
        mb_t* MbP = &MbQ->p_Vid->mb_data[pixP.mb_addr];

        int blkP = (pixP.y & ~3) + (pixP.x / 4);
        int blkQ = (idx    & ~3) + edge;

        bool intra = MbP->is_intra_block || MbQ->is_intra_block;

        if (edge == 0 && cond_bS4 && intra)
            StrValue = 4;
        else if (cond_bS3 && intra)
            StrValue = 3;
        else if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
                 (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else if (mixedModeEdgeFlag)
            StrValue = 1;
        else if (edge > 0 && (MbQ->mb_type == P16x16 || MbQ->mb_type == P16x8))
            StrValue = 0;
        else {
            int blk_x  = (pixQ.pos_x / 4);
            int blk_y  = (pixQ.pos_y + idx) / 4;
            int blk_x2 = (pixP.pos_x / 4);
            int blk_y2 = (pixP.pos_y / 4);

            StrValue = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, MbQ->p_Vid->dec_picture);
        }

        Strength[idx] = StrValue;
    }
}

void deblock_t::get_strength_hor(mb_t* MbQ, int edge)
{
    uint8_t* Strength = MbQ->strength_hor[edge];
    int StrValue;
    int yQ = (edge < BLOCK_SIZE ? edge * 4 : 1);

    slice_t* slice = MbQ->p_Slice;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

    PixelPos pixP, pixQ;
    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
    if (slice->MbaffFrameFlag) {
        getAffNeighbour(MbQ, 0, yQ - 1, mb_size, &pixP);
        getAffNeighbour(MbQ, 0, yQ    , mb_size, &pixQ);
    } else {
        getNonAffNeighbour(MbQ, 0, yQ - 1, mb_size, &pixP);
        getNonAffNeighbour(MbQ, 0, yQ    , mb_size, &pixQ);
    }
    mb_t* MbP = &MbQ->p_Vid->mb_data[pixP.mb_addr];
    MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag);

    bool verticalEdgeFlag  = 0;
    bool mixedModeEdgeFlag = slice->MbaffFrameFlag &&
                             MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag;

    bool special = MbP->p_Slice->slice_type == SI_slice || MbP->p_Slice->slice_type == SP_slice ||
                   MbQ->p_Slice->slice_type == SI_slice || MbQ->p_Slice->slice_type == SP_slice;
    bool field = slice->MbaffFrameFlag ?
        MbP->mb_field_decoding_flag || MbQ->mb_field_decoding_flag : slice->field_pic_flag;
    bool intra = MbP->is_intra_block || MbQ->is_intra_block;

    bool cond_bS4 = !field ||
                    ((slice->MbaffFrameFlag || slice->field_pic_flag) && verticalEdgeFlag);
    bool cond_bS3 = !mixedModeEdgeFlag || !verticalEdgeFlag;

    if (edge == 0 && cond_bS4 && (special || intra)) {
        StrValue = 4;
        memset(Strength, StrValue, MB_BLOCK_SIZE * sizeof(uint8_t));
        return;
    }
    if (cond_bS3 && (special || intra)) {
        StrValue = 3;
        memset(Strength, StrValue, MB_BLOCK_SIZE * sizeof(uint8_t));
        return;
    }

    if (edge > 0 && edge < BLOCK_SIZE && slice->slice_type == P_slice && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, MB_BLOCK_SIZE * sizeof(uint8_t));
        return;
    }

    for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
        int blkP = (pixP.y & ~3) + (pixP.x / 4);
        int blkQ = (yQ     & ~3) + (idx    / 4);

        if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
            (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else if (mixedModeEdgeFlag)
            StrValue = 1;
        else if (edge > 0 && edge < BLOCK_SIZE && (MbQ->mb_type == P16x16 || MbQ->mb_type == P8x16))
            StrValue = 0;
        else {
            int blk_x  = (pixQ.pos_x + idx) / 4;
            int blk_y  = (pixQ.pos_y / 4);
            int blk_x2 = (pixP.pos_x / 4);
            int blk_y2 = (pixP.pos_y / 4);

            StrValue = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, MbQ->p_Vid->dec_picture);
        }

        pixP.x     += BLOCK_SIZE;
        pixP.pos_x += BLOCK_SIZE;

        memset(Strength + idx, StrValue, BLOCK_SIZE * sizeof(uint8_t));
    }
}


// Table 8-16 Derivation of offset dependent threshold variables a' and b' from indexA and indexB

static const uint8_t TABLE_ALPHA[52] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   4,   4,   5,   6,   7,   8,   9,  10,  12,  13,
     15,  17,  20,  22,  25,  28,  32,  36,  40,  45,  50,  56,  63,
     71,  80,  90, 101, 113, 127, 144, 162, 182, 203, 226, 255, 255
};

static const uint8_t TABLE_BETA[52] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   2,   2,   2,   3,   3,   3,   3,   4,   4,   4,
      6,   6,   7,   7,   8,   8,   9,   9,  10,  10,  11,  11,  12,
     12,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18
};

// Table 8-17 Value of variable t'c0 as a function of indexA and bS

static const uint8_t TABLE_TC0[52][3] = {
    {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 },
    {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 },
    {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 },
    {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 }, {  0,  0,  0 },
    {  0,  0,  0 }, {  0,  0,  1 }, {  0,  0,  1 }, {  0,  0,  1 },
    {  0,  0,  1 }, {  0,  1,  1 }, {  0,  1,  1 }, {  1,  1,  1 },
    {  1,  1,  1 }, {  1,  1,  1 }, {  1,  1,  1 }, {  1,  1,  2 },
    {  1,  1,  2 }, {  1,  1,  2 }, {  1,  1,  2 }, {  1,  2,  3 },
    {  1,  2,  3 }, {  2,  2,  3 }, {  2,  2,  4 }, {  2,  3,  4 },
    {  2,  3,  4 }, {  3,  3,  5 }, {  3,  4,  6 }, {  3,  4,  6 },
    {  4,  5,  7 }, {  4,  5,  8 }, {  4,  6,  9 }, {  5,  7, 10 },
    {  6,  8, 11 }, {  6,  8, 13 }, {  7, 10, 14 }, {  8, 11, 16 },
    {  9, 12, 18 }, { 10, 13, 20 }, { 11, 15, 23 }, { 13, 17, 25 }
};


void deblock_t::deblock_strong(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag)
{
#define p(i) (pixP[- (i) * widthP])
#define q(i) (pixQ[  (i) * widthQ])
    bool filterSamplesFlag = bS != 0 && abs(p(0) - q(0)) < alpha
                                     && abs(p(1) - p(0)) < beta
                                     && abs(q(1) - q(0)) < beta;

    if (filterSamplesFlag && bS == 4) {
        int ap = abs(p(2) - p(0));
        int aq = abs(q(2) - q(0));
        int p0, p1, p2;
        int q0, q1, q2;

        if (chromaStyleFilteringFlag == 0 && ap < beta && abs(p(0) - q(0)) < (alpha >> 2) + 2) {
            p0 = (p(2) + 2 * p(1) + 2 * p(0) + 2 * q(0) + q(1) + 4) >> 3;
            p1 = (p(2) + p(1) + p(0) + q(0) + 2) >> 2;
            p2 = (2 * p(3) + 3 * p(2) + p(1) + p(0) + q(0) + 4) >> 3;
        } else {
            p0 = (2 * p(1) + p(0) + q(1) + 2) >> 2;
            p1 = p(1);
            p2 = p(2);
        }

        if (chromaStyleFilteringFlag == 0 && aq < beta && abs(p(0) - q(0)) < (alpha >> 2) + 2) {
            q0 = (p(1) + 2 * p(0) + 2 * q(0) + 2 * q(1) + q(2) + 4) >> 3;
            q1 = (p(0) + q(0) + q(1) + q(2) + 2) >> 2;
            q2 = (2 * q(3) + 3 * q(2) + q(1) + q(0) + p(0) + 4) >> 3;
        } else {
            q0 = (2 * q(1) + q(0) + p(1) + 2) >> 2;
            q1 = q(1);
            q2 = q(2);
        }

        p(0) = p0;
        p(1) = p1;
        p(2) = p2;
        q(0) = q0;
        q(1) = q1;
        q(2) = q2;
    }
#undef p
#undef q
}

void deblock_t::deblock_normal(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag, int tc0, int BitDepth)
{
#define p(i) (pixP[- (i) * widthP])
#define q(i) (pixQ[  (i) * widthQ])
    bool filterSamplesFlag = bS != 0 && abs(p(0) - q(0)) < alpha
                                     && abs(p(1) - p(0)) < beta
                                     && abs(q(1) - q(0)) < beta;

    if (filterSamplesFlag && bS < 4) {
        int tc, delta;
        int ap, aq;
        int p0, p1;
        int q0, q1;

        ap = abs(p(2) - p(0));
        aq = abs(q(2) - q(0));
        if (chromaStyleFilteringFlag == 0)
            tc = tc0 + (ap < beta ? 1 : 0) + (aq < beta ? 1 : 0);
        else
            tc = tc0 + 1;
        delta = clip3(-tc, tc, ((((q(0) - p(0)) << 2) + (p(1) - q(1)) + 4) >> 3));

#define Clip1(x) (clip3(0, (1 << BitDepth) - 1, x))
        p0 = Clip1(p(0) + delta);
        q0 = Clip1(q(0) - delta);
#undef Clip1

        if (chromaStyleFilteringFlag == 0 && ap < beta)
            p1 = p(1) + clip3(-tc0, tc0, (p(2) + ((p(0) + q(0) + 1) >> 1) - (p(1) << 1)) >> 1);
        else
            p1 = p(1);
        if (chromaStyleFilteringFlag == 0 && aq < beta)
            q1 = q(1) + clip3(-tc0, tc0, (q(2) + ((p(0) + q(0) + 1) >> 1) - (q(1) << 1)) >> 1);
        else
            q1 = q(1);

        p(0) = p0;
        p(1) = p1;
        q(0) = q0;
        q(1) = q1;
    }
#undef p
#undef q
}

void deblock_t::edge_loop(mb_t* MbQ, bool chromaEdgeFlag, ColorPlane pl, bool verticalEdgeFlag,
                          int edge, uint8_t* Strength)
{
    VideoParameters* p_Vid = MbQ->p_Vid;
    sps_t* sps = p_Vid->active_sps;
    bool chromaStyleFilteringFlag = chromaEdgeFlag && (sps->ChromaArrayType != 3);

    imgpel** Img = pl ? p_Vid->dec_picture->imgUV[pl-1] : p_Vid->dec_picture->imgY;
    int width    = chromaEdgeFlag == 0 ? p_Vid->dec_picture->iLumaStride : p_Vid->dec_picture->iChromaStride;
    int nE       = chromaEdgeFlag == 0 ? 16 : (verticalEdgeFlag == 1 ? sps->MbHeightC : sps->MbWidthC);
    int BitDepth = chromaEdgeFlag == 0 ? sps->BitDepthY : sps->BitDepthC;

    //auto _Strength = verticalEdgeFlag == 0 ? MbQ->strength_hor : MbQ->strength_ver;
    //auto Strength = _Strength[edge == 1 ? 4 : edge / 4];

    PixelPos pixP, pixQ;

    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };
    int* mb_size = mb_size_xy[chromaEdgeFlag == 0 ? IS_LUMA : IS_CHROMA];

    for (int pel = 0; pel < nE; ++pel) {
        if (edge == 1)
            MbQ->DeblockCall = 2;
        if (MbQ->p_Slice->MbaffFrameFlag) {
            if (verticalEdgeFlag) {
                getAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);
                getAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
            } else {
                getAffNeighbour(MbQ, pel, edge - 1, mb_size, &pixP);
                getAffNeighbour(MbQ, pel, edge    , mb_size, &pixQ);
            }
        } else {
            if (verticalEdgeFlag) {
                getNonAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);
                getNonAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
            } else {
                getNonAffNeighbour(MbQ, pel, edge - 1, mb_size, &pixP);
                getNonAffNeighbour(MbQ, pel, edge    , mb_size, &pixQ);
            }
        }
        if (edge == 1)
            MbQ->DeblockCall = 1;

        mb_t* MbP = &p_Vid->mb_data[pixP.mb_addr];
        int StrengthIdx;
        if (MbQ->p_Slice->MbaffFrameFlag)
            StrengthIdx = (nE == 8) ? ((MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) ? pel << 1 : ((pel >> 1) << 2) + (pel & 0x01)) : pel;
        else
            StrengthIdx = (nE == 8) ? pel << 1 : pel;
        int bS = Strength[StrengthIdx];

        if (bS > 0) {
            imgpel* SrcPtrP = &Img[pixP.pos_y][pixP.pos_x];
            imgpel* SrcPtrQ = &Img[pixQ.pos_y][pixQ.pos_x];

            int incP, incQ;
            if (verticalEdgeFlag) {
                incP = 1;
                incQ = 1;
            } else {
                incP = (MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag ? 2 * width : width);
                incQ = (MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag ? 2 * width : width);
            }

            int qPp    = pl ? MbP->QpC[pl - 1] : MbP->QpY;
            int qPq    = pl ? MbQ->QpC[pl - 1] : MbQ->QpY;
            int qPav   = (qPp + qPq + 1) >> 1;
            int indexA = clip3(0, 51, qPav + MbQ->p_Slice->FilterOffsetA);
            int indexB = clip3(0, 51, qPav + MbQ->p_Slice->FilterOffsetB);
            int alpha  = TABLE_ALPHA[indexA] * (1 << (BitDepth - 8));
            int beta   = TABLE_BETA [indexB] * (1 << (BitDepth - 8));

            if (bS == 4)
                this->deblock_strong(SrcPtrP, SrcPtrQ, incP, incQ, alpha, beta, bS, chromaStyleFilteringFlag);
            else if (bS > 0) {
                int tc0 = TABLE_TC0[indexA][bS - 1] * (1 << (BitDepth - 8));
                this->deblock_normal(SrcPtrP, SrcPtrQ, incP, incQ, alpha, beta, bS, chromaStyleFilteringFlag, tc0, BitDepth);
            }
        }
    }
}


void deblock_t::deblock_mb(mb_t* mb)
{
    mb_t* MbQ      = mb;
    slice_t* slice = MbQ->p_Slice;
    sps_t*   sps   = slice->active_sps;

    // return, if filter is disabled
    if (slice->disable_deblocking_filter_idc != 1) {
        MbQ->DeblockCall = 1;
        short mb_x, mb_y;
        int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
        get_mb_pos(slice->p_Vid, MbQ->mbAddrX, mb_size, &mb_x, &mb_y);

        if (MbQ->mb_type == I8MB)
            assert(MbQ->transform_size_8x8_flag);

        bool fieldMbInFrameFlag      = slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag;
        bool filterInternalEdgesFlag = slice->disable_deblocking_filter_idc != 1;
        bool filterLeftMbEdgeFlag = mb_x != 0;
        bool filterTopMbEdgeFlag  = mb_y != 0;
//        bool filterLeftMbEdgeFlag    = (slice->disable_deblocking_filter_idc != 1) && (mb_x != 0) &&
//                                       (slice->disable_deblocking_filter_idc != 2 || MbQ->mbAvailA);
//        bool filterTopMbEdgeFlag     = (slice->disable_deblocking_filter_idc != 1) && (mb_y != 0) &&
//                                       (slice->disable_deblocking_filter_idc != 2 || MbQ->mbAvailB);

//      filterLeftMbEdgeFlag = (MbaffFrameFlag == 0 && CurrMbAddr % PicWidthInMbs == 0) ||
//                             (MbaffFrameFlag == 1 && (CurrMbAddr >> 1) % PicWidthInMbs == 0) ||
//                             (slice->disable_deblocking_filter_idc == 1) ||
//                             (slice->disable_deblocking_filter_idc == 2 && mbAddrA == NULL) ? 0 : 1

//      filterTopMbEdgeFlag = (MbaffFrameFlag == 0 && CurrMbAddr < PicWidthInMbs) ||
//                            (MbaffFrameFlag == 1 && (CurrMbAddr >> 1) < PicWidthInMbs && mb->field) ||
//                            (MbaffFrameFlag == 1 && (CurrMbAddr >> 1) < PicWidthInMbs && mb->frame && CurrMbAddr % 2 == 0) ||
//                            (slice->disable_deblocking_filter_idc == 1) ||
//                            (slice->disable_deblocking_filter_idc == 2 && mbAddrB == NULL) ? 0 : 1

        if (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag&& mb_y == MB_BLOCK_SIZE)
            filterTopMbEdgeFlag = 0;

        if (slice->disable_deblocking_filter_idc == 2) {
            // don't filter at slice boundaries
            filterLeftMbEdgeFlag = MbQ->mbAvailA;
            filterTopMbEdgeFlag  = MbQ->mbAvailB;
            if (slice->MbaffFrameFlag && !MbQ->mb_field_decoding_flag && (MbQ->mbAddrX & 0x01))
                filterTopMbEdgeFlag = 1;
        }

        bool filterVerEdgeFlag[2][4] = {
            { filterLeftMbEdgeFlag,
              filterInternalEdgesFlag && !MbQ->transform_size_8x8_flag,
              filterInternalEdgesFlag,
              filterInternalEdgesFlag && !MbQ->transform_size_8x8_flag },
            { filterLeftMbEdgeFlag,
              filterInternalEdgesFlag && (sps->ChromaArrayType != 3 || !MbQ->transform_size_8x8_flag),
              filterInternalEdgesFlag && (sps->ChromaArrayType == 3),
              filterInternalEdgesFlag && (sps->ChromaArrayType == 3 && !MbQ->transform_size_8x8_flag) }
        };
        bool filterHorEdgeFlag[2][4] = {
            { filterTopMbEdgeFlag,
              filterInternalEdgesFlag && !MbQ->transform_size_8x8_flag,
              filterInternalEdgesFlag,
              filterInternalEdgesFlag && !MbQ->transform_size_8x8_flag },
            { filterTopMbEdgeFlag,
              filterInternalEdgesFlag && (sps->ChromaArrayType != 3 || !MbQ->transform_size_8x8_flag),
              filterInternalEdgesFlag && (sps->ChromaArrayType != 1),
              filterInternalEdgesFlag && (sps->ChromaArrayType == 2 ||
                                          (sps->ChromaArrayType == 3 && !MbQ->transform_size_8x8_flag)) }
        };

        if (slice->MbaffFrameFlag)
            CheckAvailabilityOfNeighborsMBAFF(MbQ);

        bool mixedModeEdgeFlag = 0;
        for (int edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            //if ((MbQ->CodedBlockPatternLuma == 0 && MbQ->CodedBlockPatternChroma == 0) &&
            //    (slice->slice_type == P_slice || slice->slice_type == B_slice)) {
            //    if (filterNon8x8LumaEdgesFlag[edge] == 0 && sps->chroma_format_idc != YUV444)
            //        continue;
            //}
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            //if ((MbQ->CodedBlockPatternLuma == 0 && MbQ->CodedBlockPatternChroma == 0) && 
            //    (slice->slice_type == P_slice || slice->slice_type == B_slice)) {
            //    if (filterNon8x8LumaEdgesFlag[edge] == 0 && sps->chroma_format_idc == YUV420)
            //        continue;
            //}

            if (filterVerEdgeFlag[0][edge])
                this->get_strength_ver(MbQ, edge);
            if (filterHorEdgeFlag[0][edge])
                this->get_strength_hor(MbQ, edge);

            if (filterHorEdgeFlag[0][edge] && edge == 0 && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag) {
                mixedModeEdgeFlag = MbQ->mixedModeEdgeFlag;
                MbQ->DeblockCall = 2;
                this->get_strength_hor(MbQ, 4);
                MbQ->DeblockCall = 1;
            }
        }
        MbQ->mixedModeEdgeFlag = mixedModeEdgeFlag;

        // Vertical deblocking
        for (int edge = 0; edge < 4; ++edge) {
            if (filterVerEdgeFlag[0][edge]) {
                uint8_t* Strength = MbQ->strength_ver[edge];
                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1))))
                    this->edge_loop(MbQ, false, PLANE_Y, true, edge * 4, Strength);
            }
            if (sps->ChromaArrayType != 0 && filterVerEdgeFlag[1][edge]) {
                uint8_t* Strength = MbQ->strength_ver[sps->ChromaArrayType < 3 ? edge * 2 : edge];
                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) {
                    this->edge_loop(MbQ, true, PLANE_U, true, edge * 4, Strength);
                    this->edge_loop(MbQ, true, PLANE_V, true, edge * 4, Strength);
                }
            }
        }

        // horizontal deblocking
        for (int edge = 0; edge < 4; ++edge) {
            if (filterHorEdgeFlag[0][edge]) {
                uint8_t* Strength = MbQ->strength_hor[edge];
                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1))))
                    this->edge_loop(MbQ, false, PLANE_Y, false, edge * 4, Strength);
                if (edge == 0 && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag) {
                    uint8_t* Strength = MbQ->strength_hor[4];
                    if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1))))
                        this->edge_loop(MbQ, false, PLANE_Y, false, 1, Strength);
                }
            }
            if (sps->ChromaArrayType != 0 && filterHorEdgeFlag[1][edge]) {
                uint8_t* Strength = MbQ->strength_hor[sps->ChromaArrayType < 2 ? edge * 2 : edge];
                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) {
                    this->edge_loop(MbQ, true, PLANE_U, false, edge * 4, Strength);
                    this->edge_loop(MbQ, true, PLANE_V, false, edge * 4, Strength);
                }
                if (edge == 0 && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag) {
                    uint8_t* Strength = MbQ->strength_hor[4];
                    if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) {
                        this->edge_loop(MbQ, true, PLANE_U, false, 1, Strength);
                        this->edge_loop(MbQ, true, PLANE_V, false, 1, Strength);
                    }
                }
            }
        }
    }

    MbQ->DeblockCall = 0;
}

void deblock_t::deblock_pic(VideoParameters *p_Vid, storable_picture *p)
{
    this->init_neighbors(p_Vid);

    for (int mbAddr = 0; mbAddr < p->PicSizeInMbs; ++mbAddr)
        this->deblock_mb(&p_Vid->mb_data[mbAddr]);
}


void deblock_t::init_neighbors(VideoParameters *p_Vid)
{
    slice_t* slice  = p_Vid->ppSliceList[0];
    sps_t*   sps    = p_Vid->active_sps;
    int      width  = sps->PicWidthInMbs;
    int      height = sps->FrameHeightInMbs / (1 + slice->field_pic_flag);
    int      size   = width * height;

    mb_t* currMB = &p_Vid->mb_data[0];
    // do the top left corner
    currMB->mbup   = NULL;
    currMB->mbleft = NULL;
    currMB++;
    // do top row
    for (int i = 1; i < width; i++) {
        currMB->mbup   = NULL;
        currMB->mbleft = currMB - 1;
        currMB++;
    }

    // do left edge
    for (int i = width; i < size; i += width) {
        currMB->mbup   = currMB - width;
        currMB->mbleft = NULL;   
        currMB += width;
    }
    // do all others
    for (int j = width + 1; j < width * height + 1; j += width) {
        currMB = &p_Vid->mb_data[j];
        for (int i = 1; i < width; i++) {
            currMB->mbup   = currMB - width;
            currMB->mbleft = currMB - 1;
            currMB++;
        }
    }
}

void deblock_t::make_frame_picture_JV(VideoParameters *p_Vid)
{
    sps_t* sps = p_Vid->active_sps;
    p_Vid->dec_picture = p_Vid->dec_picture_JV[0];

    if (p_Vid->dec_picture->used_for_reference) {
        int nsize = (p_Vid->dec_picture->size_y / BLOCK_SIZE) *
                    (p_Vid->dec_picture->size_x / BLOCK_SIZE) * sizeof(pic_motion_params);
        memcpy(&p_Vid->dec_picture->JVmv_info[PLANE_Y][0][0], &p_Vid->dec_picture_JV[PLANE_Y]->mv_info[0][0], nsize);
        memcpy(&p_Vid->dec_picture->JVmv_info[PLANE_U][0][0], &p_Vid->dec_picture_JV[PLANE_U]->mv_info[0][0], nsize);
        memcpy(&p_Vid->dec_picture->JVmv_info[PLANE_V][0][0], &p_Vid->dec_picture_JV[PLANE_V]->mv_info[0][0], nsize);
    }

    // This could be done with pointers and seems not necessary
    for (int uv = 0; uv < 2; uv++) {
        for (int line = 0; line < sps->FrameHeightInMbs * 16; line++) {
            int nsize = sizeof(imgpel) * sps->PicWidthInMbs * 16;
            memcpy(p_Vid->dec_picture->imgUV[uv][line], p_Vid->dec_picture_JV[uv+1]->imgY[line], nsize);
        }
        free_storable_picture(p_Vid->dec_picture_JV[uv+1]);
    }
}

void deblock_t::deblock(VideoParameters *p_Vid, storable_picture *p)
{
    int iDeblockMode = 1;
    //init mb_data;
    for (int j = 0; j < p_Vid->iSliceNumOfCurrPic; j++) {
        if (p_Vid->ppSliceList[j]->disable_deblocking_filter_idc != 1)
            iDeblockMode = 0;
#if (MVC_EXTENSION_ENABLE)
        assert(p_Vid->ppSliceList[j]->view_id == p_Vid->ppSliceList[0]->view_id);
#endif
    }

    if (!iDeblockMode && (0x03 & (1 << p->used_for_reference))) {
        //deblocking for frame or field
        if (p_Vid->active_sps->separate_colour_plane_flag) {
            int colour_plane_id = p_Vid->ppSliceList[0]->colour_plane_id;
            for (int nplane = 0; nplane < MAX_PLANE; ++nplane) {
                p_Vid->ppSliceList[0]->colour_plane_id = nplane;
                p_Vid->mb_data     = p_Vid->mb_data_JV    [nplane];
                p_Vid->dec_picture = p_Vid->dec_picture_JV[nplane];
                this->deblock_pic(p_Vid, p);
            }
            p_Vid->ppSliceList[0]->colour_plane_id = colour_plane_id;
        } else
            this->deblock_pic(p_Vid, p);
    }

    if (p_Vid->active_sps->separate_colour_plane_flag)
        this->make_frame_picture_JV(p_Vid);
}


}
}
