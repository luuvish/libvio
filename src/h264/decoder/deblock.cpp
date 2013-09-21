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


inline int deblock_t::compare_mvs(const MotionVector* mv0, const MotionVector* mv1, int mvlimit)
{
    return (abs(mv0->mv_x - mv1->mv_x) >= 4) | (abs(mv0->mv_y - mv1->mv_y) >= mvlimit);
}

inline int deblock_t::bs_compare_mvs(const pic_motion_params* mv_info_p, const pic_motion_params* mv_info_q, int mvlimit)
{
    int StrValue;

    storable_picture* ref_p0 = mv_info_p->ref_pic[LIST_0];
    storable_picture* ref_q0 = mv_info_q->ref_pic[LIST_0];
    storable_picture* ref_p1 = mv_info_p->ref_pic[LIST_1];
    storable_picture* ref_q1 = mv_info_q->ref_pic[LIST_1];

    if ((ref_p0 == ref_q0 && ref_p1 == ref_q1) || (ref_p0 == ref_q1 && ref_p1 == ref_q0)) {
        // L0 and L1 reference pictures of p0 are different; q0 as well
        if (ref_p0 != ref_p1) {
            // compare MV for the same reference picture
            if (ref_p0 == ref_q0) {
                StrValue = 
                    this->compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                    this->compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit);
            } else {
                StrValue =
                    this->compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                    this->compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit);
            }
        } else { // L0 and L1 reference pictures of p0 are the same; q0 as well
            StrValue = ((
                this->compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) |
                this->compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                & (
                this->compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                this->compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)
            ));
        }
    } else
        StrValue = 1;

    return StrValue;
}


void deblock_t::strength_vertical(mb_t* MbQ, int edge)
{
    uint8_t* Strength = MbQ->strength_ver[edge];
    int StrValue;

    slice_t* slice = MbQ->p_Slice;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;
    pic_motion_params** mv_info = MbQ->p_Vid->dec_picture->mv_info;

    int xI, yI;
    neighbour.get_mb2pos(MbQ->p_Slice, MbQ->mbAddrX, xI, yI);
    int xQ = xI + edge * 4;
    int yQ = MbQ->fieldMbInFrameFlag ? yI / 32 * 32 + (yI % 32) / 2 + 16 * (yI & 1) : yI;

    int mbAddrP;
    neighbour.get_pos2mb(MbQ->p_Slice, xI + (edge - 1) * 4, yI, mbAddrP);
    mb_t* MbP = &MbQ->p_Vid->mb_data[mbAddrP];
    int xP = xQ - 4;
    int yP = yQ;

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

    int blkP_x = (xP & 15) / 4;
    int blkQ   = edge;

    for (int y = 0; y < MB_BLOCK_SIZE; ++y) {
        int y0 = yI + (1 + MbQ->fieldMbInFrameFlag) * y;
        neighbour.get_pos2mb(MbQ->p_Slice, xI + (edge - 1) * 4, y0, mbAddrP);
        MbP = &MbQ->p_Vid->mb_data[mbAddrP];
        bool fieldMbInFrameFlag = slice->MbaffFrameFlag && MbP->mb_field_decoding_flag;
        yP = fieldMbInFrameFlag ? y0 / 32 * 32 + (y0 % 32) / 2 + 16 * (y0 & 1) : y0;

        int blkP = (yP & 12) + blkP_x;

        bool intra = MbP->is_intra_block || MbQ->is_intra_block;

        if (edge == 0 && cond_bS4 && intra)
            StrValue = 4;
        else if (cond_bS3 && intra)
            StrValue = 3;
        else if ((MbQ->cbp_blks[0] & ((uint64_t)1 << (blkQ + (y & ~3)))) != 0 ||
                 (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else if (mixedModeEdgeFlag)
            StrValue = 1;
        else if (edge > 0 && (MbQ->mb_type == P16x16 || MbQ->mb_type == P16x8))
            StrValue = 0;
        else {
            auto mv_info_p = &mv_info[yQ / 4 + y / 4][xQ / 4];
            auto mv_info_q = &mv_info[yP / 4][xP / 4];

            StrValue = this->bs_compare_mvs(mv_info_p, mv_info_q, mvlimit);
        }

        Strength[y] = StrValue;
    }
}

void deblock_t::strength_horizontal(mb_t* MbQ, int edge)
{
    uint8_t* Strength = MbQ->strength_hor[edge];
    int StrValue;

    slice_t* slice = MbQ->p_Slice;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;
    pic_motion_params** mv_info = MbQ->p_Vid->dec_picture->mv_info;

    int dy = (1 + MbQ->fieldMbInFrameFlag);

    int xI, yI;
    neighbour.get_mb2pos(MbQ->p_Slice, MbQ->mbAddrX, xI, yI);
    int xQ = xI;
    int y0 = yI + (edge == 4 ? 1 : dy * edge * 4);
    int yQ = MbQ->fieldMbInFrameFlag ? y0 / 32 * 32 + (y0 % 32) / 2 + 16 * (y0 & 1) : y0;

    int mbAddrP;
    neighbour.get_pos2mb(MbQ->p_Slice, xQ, y0 - 4 * dy, mbAddrP);
    mb_t* MbP = &MbQ->p_Vid->mb_data[mbAddrP];
    int xP = xQ;
    int yP = y0 - 4 * dy;
    yP = MbQ->fieldMbInFrameFlag ? yP / 32 * 32 + (yP % 32) / 2 + 16 * (yP & 1) : yP;

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

    int blkP   = (yP & 12);
    int p_blky = (yP / 4);
    int blkQ   = (yQ & 12);

    for (int x = 0; x < MB_BLOCK_SIZE; x += BLOCK_SIZE) {
        if ((MbQ->cbp_blks[0] & ((uint64_t)1 << (blkQ + x / 4))) != 0 ||
            (MbP->cbp_blks[0] & ((uint64_t)1 << (blkP + x / 4))) != 0)
            StrValue = 2;
        else if (mixedModeEdgeFlag)
            StrValue = 1;
        else if (edge > 0 && edge < BLOCK_SIZE && (MbQ->mb_type == P16x16 || MbQ->mb_type == P8x16))
            StrValue = 0;
        else {
            auto mv_info_p = &mv_info[yQ / 4][xQ / 4 + x / 4];
            auto mv_info_q = &mv_info[p_blky][xP / 4 + x / 4];

            StrValue = this->bs_compare_mvs(mv_info_p, mv_info_q, mvlimit);
        }

        memset(Strength + x, StrValue, BLOCK_SIZE * sizeof(uint8_t));
    }
}

void deblock_t::strength(mb_t* MbQ)
{
    slice_t* slice = MbQ->p_Slice;
    sps_t* sps = slice->active_sps;

    if (slice->disable_deblocking_filter_idc != 1) {
        MbQ->fieldMbInFrameFlag = slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag;

        bool filterLeftMbEdgeFlag = 0;
        bool filterTopMbEdgeFlag  = 0;

        int xI, yI;
        neighbour.get_mb2pos(MbQ->p_Slice, MbQ->mbAddrX, xI, yI);

        int mbAddrL;
        mb_t* MbL = neighbour.get_pos2mb(MbQ->p_Slice, xI - 16, yI, mbAddrL);
        if (MbL) {
            if (slice->disable_deblocking_filter_idc == 0)
                filterLeftMbEdgeFlag = 1;
            if (slice->disable_deblocking_filter_idc == 2)
                filterLeftMbEdgeFlag = MbQ->slice_nr == slice->mb_data[mbAddrL].slice_nr;
        }
        int mbAddrU;
        mb_t* MbU = neighbour.get_pos2mb(MbQ->p_Slice, xI, yI - 16, mbAddrU);
        if (MbU) {
            if (slice->disable_deblocking_filter_idc == 0)
                filterTopMbEdgeFlag = 1;
            if (slice->disable_deblocking_filter_idc == 2)
                filterTopMbEdgeFlag = MbQ->slice_nr == slice->mb_data[mbAddrU].slice_nr;
        }

        bool filterInternalEdgesFlag = slice->disable_deblocking_filter_idc != 1;

        for (int dir = 0; dir < 2; dir++) {
            MbQ->filterVerEdgeFlag[dir][0] = filterLeftMbEdgeFlag;
            MbQ->filterHorEdgeFlag[dir][0] = filterTopMbEdgeFlag;
            for (int edge = 1; edge < 4; edge++) {
                MbQ->filterVerEdgeFlag[dir][edge] = filterInternalEdgesFlag;
                MbQ->filterHorEdgeFlag[dir][edge] = filterInternalEdgesFlag;
            }
        }

        if (MbQ->transform_size_8x8_flag) {
            MbQ->filterVerEdgeFlag[0][1] = MbQ->filterVerEdgeFlag[0][3] = 0;
            MbQ->filterHorEdgeFlag[0][1] = MbQ->filterHorEdgeFlag[0][3] = 0;
        }
        if (sps->ChromaArrayType == 1) {
            MbQ->filterVerEdgeFlag[1][2] = MbQ->filterVerEdgeFlag[1][3] = 0;
            MbQ->filterHorEdgeFlag[1][2] = MbQ->filterHorEdgeFlag[1][3] = 0;
        } else if (sps->ChromaArrayType == 2) {
            MbQ->filterVerEdgeFlag[1][2] = MbQ->filterVerEdgeFlag[1][3] = 0;
        } else if (sps->ChromaArrayType == 3 && MbQ->transform_size_8x8_flag) {
            MbQ->filterVerEdgeFlag[1][1] = MbQ->filterVerEdgeFlag[1][3] = 0;
            MbQ->filterHorEdgeFlag[1][1] = MbQ->filterHorEdgeFlag[1][3] = 0;
        }

        bool mixedModeEdgeFlag = 0;
        for (int edge = 0; edge < 4; ++edge) {
            if (MbQ->filterVerEdgeFlag[0][edge])
                this->strength_vertical(MbQ, edge);
            if (MbQ->filterHorEdgeFlag[0][edge])
                this->strength_horizontal(MbQ, edge);

            if (MbQ->filterHorEdgeFlag[0][edge] &&
                edge == 0 && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag) {
                mixedModeEdgeFlag = MbQ->mixedModeEdgeFlag;
                this->strength_horizontal(MbQ, 4);
            }
        }
        MbQ->mixedModeEdgeFlag = mixedModeEdgeFlag;
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


void deblock_t::filter_strong(imgpel *pixQ, int width, int alpha, int beta, int bS, bool chromaStyleFilteringFlag)
{
#define p(i) (pixQ[- (i + 1) * width])
#define q(i) (pixQ[  (i    ) * width])
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

void deblock_t::filter_normal(imgpel *pixQ, int width, int alpha, int beta, int bS, bool chromaStyleFilteringFlag, int tc0, int BitDepth)
{
#define p(i) (pixQ[- (i + 1) * width])
#define q(i) (pixQ[  (i    ) * width])
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


void deblock_t::filter_edge(mb_t* MbQ, bool chromaEdgeFlag, ColorPlane pl, bool verticalEdgeFlag, bool fieldModeInFrameFilteringFlag, int edge)
{
    VideoParameters* p_Vid = MbQ->p_Vid;
    sps_t* sps = p_Vid->active_sps;
    bool chromaStyleFilteringFlag = chromaEdgeFlag && (sps->ChromaArrayType != 3);

    imgpel** Img = pl ? p_Vid->dec_picture->imgUV[pl-1] : p_Vid->dec_picture->imgY;
    int width    = chromaEdgeFlag == 0 ? p_Vid->dec_picture->iLumaStride : p_Vid->dec_picture->iChromaStride;
    int nE       = chromaEdgeFlag == 0 ? 16 : (verticalEdgeFlag == 1 ? sps->MbHeightC : sps->MbWidthC);
    int BitDepth = chromaEdgeFlag == 0 ? sps->BitDepthY : sps->BitDepthC;

    uint8_t* Strength;
    if (verticalEdgeFlag)
        Strength = MbQ->strength_ver[edge == 1 ? 4 : edge * 4 / (chromaEdgeFlag == 0 ? 16 : sps->MbWidthC)];
    else
        Strength = MbQ->strength_hor[edge == 1 ? 4 : edge * 4 / (chromaEdgeFlag == 0 ? 16 : sps->MbHeightC)];
    uint64_t* bS64 = (uint64_t*)Strength;
    if (bS64[0] == 0 && bS64[1] == 0)
        return;

    int dy = (1 + fieldModeInFrameFilteringFlag);

    int xI, yI;
    neighbour.get_mb2pos(MbQ->p_Slice, MbQ->mbAddrX, xI, yI);
    int xP = chromaEdgeFlag == 0 ? xI : (xI / sps->SubWidthC);
    int yP = chromaEdgeFlag == 0 ? yI : (yI + sps->SubHeightC - 1) / sps->SubHeightC;

    int mbAddrJ;
    int xJ = xI, yJ = yI;
    if (verticalEdgeFlag)
        xJ += (edge - 1) * (chromaEdgeFlag == 0 ? 1 : sps->SubWidthC);
    else
        yJ += dy * (edge - 1) * (chromaEdgeFlag == 0 ? 1 : sps->SubHeightC) - (edge % 2);
    neighbour.get_pos2mb(MbQ->p_Slice, xJ, yJ, mbAddrJ);
    mb_t* MbP = &p_Vid->mb_data[mbAddrJ];

    int incQ = verticalEdgeFlag ? 1 : dy * width;
    int nxtQ = verticalEdgeFlag ? dy * width : 1;

    imgpel* SrcPtrQ = verticalEdgeFlag ? &Img[yP][xP + edge] : &Img[yP + dy * edge - (edge % 2)][xP];

    bool mixed = verticalEdgeFlag && (!MbQ->mb_field_decoding_flag && MbP->mb_field_decoding_flag);

    for (int pel = 0; pel < nE; ++pel) {
        int StrengthIdx = (nE == 8) ? (pel << 1) + (mixed && (pel & 1)) : pel;
        int bS = Strength[StrengthIdx];

        if (bS > 0) {
            if (verticalEdgeFlag) {
                yJ = yI + dy * pel * (chromaEdgeFlag == 0 ? 1 : sps->SubHeightC);
                neighbour.get_pos2mb(MbQ->p_Slice, xJ, yJ + (chromaEdgeFlag && mixed && (pel & 1)), mbAddrJ);
                MbP = &p_Vid->mb_data[mbAddrJ];
            }

            int qPp    = pl ? MbP->QpC[pl - 1] : MbP->QpY;
            int qPq    = pl ? MbQ->QpC[pl - 1] : MbQ->QpY;
            int qPav   = (qPp + qPq + 1) >> 1;
            int indexA = clip3(0, 51, qPav + MbQ->p_Slice->FilterOffsetA);
            int indexB = clip3(0, 51, qPav + MbQ->p_Slice->FilterOffsetB);
            int alpha  = TABLE_ALPHA[indexA] * (1 << (BitDepth - 8));
            int beta   = TABLE_BETA [indexB] * (1 << (BitDepth - 8));

            if (bS == 4)
                this->filter_strong(SrcPtrQ, incQ, alpha, beta, bS, chromaStyleFilteringFlag);
            else if (bS > 0) {
                int tc0 = TABLE_TC0[indexA][bS - 1] * (1 << (BitDepth - 8));
                this->filter_normal(SrcPtrQ, incQ, alpha, beta, bS, chromaStyleFilteringFlag, tc0, BitDepth);
            }
        }
        SrcPtrQ += nxtQ;
    }
}

void deblock_t::filter_vertical(mb_t* MbQ)
{
    slice_t* slice = MbQ->p_Slice;
    sps_t* sps = slice->active_sps;

    if (slice->disable_deblocking_filter_idc != 1) {
        for (int edge = 0; edge < 4; ++edge) {
            if (MbQ->filterVerEdgeFlag[0][edge])
                this->filter_edge(MbQ, false, PLANE_Y, true, MbQ->fieldMbInFrameFlag, edge * 4);
            if (sps->ChromaArrayType != 0 && MbQ->filterVerEdgeFlag[1][edge]) {
                this->filter_edge(MbQ, true, PLANE_U, true, MbQ->fieldMbInFrameFlag, edge * 4);
                this->filter_edge(MbQ, true, PLANE_V, true, MbQ->fieldMbInFrameFlag, edge * 4);
            }
        }
    }
}

void deblock_t::filter_horizontal(mb_t* MbQ)
{
    slice_t* slice = MbQ->p_Slice;
    sps_t* sps = slice->active_sps;

    if (slice->disable_deblocking_filter_idc != 1) {
        for (int edge = 0; edge < 4; ++edge) {
            bool mixed = (edge == 0) && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag;

            if (MbQ->filterHorEdgeFlag[0][edge]) {
                if (!mixed)
                    this->filter_edge(MbQ, false, PLANE_Y, false, MbQ->fieldMbInFrameFlag, edge * 4);
                else {
                    this->filter_edge(MbQ, false, PLANE_Y, false, true, 0);
                    this->filter_edge(MbQ, false, PLANE_Y, false, true, 1);
                }
            }
            if (sps->ChromaArrayType != 0 && MbQ->filterHorEdgeFlag[1][edge]) {
                if (!mixed) {
                    this->filter_edge(MbQ, true, PLANE_U, false, MbQ->fieldMbInFrameFlag, edge * 4);
                    this->filter_edge(MbQ, true, PLANE_V, false, MbQ->fieldMbInFrameFlag, edge * 4);
                } else {
                    this->filter_edge(MbQ, true, PLANE_U, false, true, 0);
                    this->filter_edge(MbQ, true, PLANE_U, false, true, 1);
                    this->filter_edge(MbQ, true, PLANE_V, false, true, 0);
                    this->filter_edge(MbQ, true, PLANE_V, false, true, 1);
                }
            }
        }
    }
}

void deblock_t::deblock_pic(VideoParameters* p_Vid)
{
    slice_t* slice = p_Vid->ppSliceList[0];
    mb_t* mb_data = slice->mb_data;

    for (int mbAddr = 0; mbAddr < slice->PicSizeInMbs; ++mbAddr) {
        mb_t* mb = &mb_data[mbAddr];
        this->strength(mb);
    }
    for (int mbAddr = 0; mbAddr < slice->PicSizeInMbs; ++mbAddr) {
        mb_t* mb = &mb_data[mbAddr];
        this->filter_vertical(mb);
        this->filter_horizontal(mb);
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

static void update_mbaff_macroblock_data(imgpel **cur_img, imgpel (*temp)[16], int x0, int width, int height)
{
    imgpel (*temp_evn)[16] = temp;
    imgpel (*temp_odd)[16] = temp + height; 
    imgpel **temp_img = cur_img;
    int y;

    for (y = 0; y < 2 * height; ++y)
        memcpy(*temp++, (*temp_img++ + x0), width * sizeof(imgpel));

    for (y = 0; y < height; ++y) {
        memcpy((*cur_img++ + x0), *temp_evn++, width * sizeof(imgpel));
        memcpy((*cur_img++ + x0), *temp_odd++, width * sizeof(imgpel));
    }
}

static void MbAffPostProc(VideoParameters *p_Vid)
{
    slice_t* slice = p_Vid->ppSliceList[0];
    sps_t* sps = p_Vid->active_sps;

    storable_picture* dec_picture = p_Vid->dec_picture;
    imgpel**  imgY  = dec_picture->imgY;
    imgpel*** imgUV = dec_picture->imgUV;

    imgpel temp_buffer[32][16];
    int x0, y0;

    for (int mbAddr = 0; mbAddr < dec_picture->PicSizeInMbs; mbAddr += 2) {
        if (dec_picture->motion.mb_field_decoding_flag[mbAddr]) {
            neighbour.get_mb2pos(slice, mbAddr, x0, y0);
            update_mbaff_macroblock_data(imgY + y0, temp_buffer, x0, MB_BLOCK_SIZE, MB_BLOCK_SIZE);

            if (dec_picture->chroma_format_idc != YUV400) {
                x0 = (short) ((x0 * sps->MbWidthC ) >> 4);
                y0 = (short) ((y0 * sps->MbHeightC) >> 4);

                update_mbaff_macroblock_data(imgUV[0] + y0, temp_buffer, x0, sps->MbWidthC, sps->MbHeightC);
                update_mbaff_macroblock_data(imgUV[1] + y0, temp_buffer, x0, sps->MbWidthC, sps->MbHeightC);
            }
        }
    }
}

void deblock_t::deblock(VideoParameters *p_Vid)
{
    storable_picture* p = p_Vid->dec_picture;
    if (p->mb_aff_frame_flag)
        MbAffPostProc(p_Vid);

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
                this->deblock_pic(p_Vid);
            }
            p_Vid->ppSliceList[0]->colour_plane_id = colour_plane_id;
        } else
            this->deblock_pic(p_Vid);
    }

    if (p_Vid->active_sps->separate_colour_plane_flag)
        this->make_frame_picture_JV(p_Vid);
}


}
}
