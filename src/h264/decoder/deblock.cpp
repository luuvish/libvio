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


deblock_t deblock;


#define MAX_QP          51


static mb_t* get_non_aff_neighbor_luma(mb_t *mb, int xN, int yN)
{
    if (xN < 0)
        return mb->mbleft;
    else if (yN < 0)
        return mb->mbup;
    else
        return mb;
}

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


static void get_strength_ver(mb_t *MbQ, int edge)
{
    uint8_t*  Strength = MbQ->strength_ver[edge];
    slice_t*  slice = MbQ->p_Slice;
    int       StrValue;
    BlockPos* PicPos = MbQ->p_Vid->PicPos;
    storable_picture* p = MbQ->p_Vid->dec_picture;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

    int xQ = (edge << 2) - 1;
    mb_t* neighbor = get_non_aff_neighbor_luma(MbQ, xQ, 0);
    mb_t* MbP = edge ? MbQ : neighbor;

/*
    mixedModeEdgeFlag = MbaffFrameFlag == 1 && p->field != q->field

    if (!p->field && !q->field && (p->intra || q->intra))
        bS = 4;
    if (!p->field && !q->field && (p->slice_type == SP || q->slice_type == SI))
        bS = 4;
    if ((MbaffFrameFlag == 1 || field_pic_flag == 1) && verticalEdgeFlag == 1 && (p->intra || q->intra))
        bS = 4;
    if ((MbaffFrameFlag == 1 || field_pic_flag == 1) && verticalEdgeFlag == 1 && (p->slice_type == SP || q->slice_type == SI))
        bS = 4;

    if (mixedModeEdgeFlag == 0 && (p->intra || q->intra))
        bS = 3;
    if (mixedModeEdgeFlag == 0 && (p->slice_type == SP || q->slice_type == SI))
        bS = 3;
    if (mixedModeEdgeFlag == 1 && verticalEdgeFlag == 0 && (p->intra || q->intra))
        bS = 3;
    if (mixedModeEdgeFlag == 1 && verticalEdgeFlag == 0 && (p->slice_type == SP || q->slice_type == SI))
        bS = 3;

    if (p->transform_size_8x8_flag == 1 && p->cbp[0..3] != 0)
        bS = 2;
    if (p->transform_size_8x8_flag == 0 && p->cbp[0..15] != 0)
        bS = 2;
    if (q->transform_size_8x8_flag == 1 && q->cbp[0..3] != 0)
        bS = 2;
    if (q->transform_size_8x8_flag == 0 && q->cbp[0..15] != 0)
        bS = 2;

    if (mixedModeEdgeFlag == 1)
        bS = 1;
    if (mixedModeEdgeFlag == 0 && p->ref_pic_list != q->ref_pic_list)
        bS = 1;
    if (mixedModeEdgeFlag == 0 && abs(p->mv[0] - q->mv[0]) >= 4)
        bS = 1;
    if (mixedModeEdgeFlag == 0 && abs(p->mv[0] - q->mv[0]) >= 4 && abs(p->mv[1] - q->mv[1]) >= 4)
        bS = 1;
    if (mixedModeEdgeFlag == 0 && abs(p->mv[0] - q->mv[1]) >= 4 && abs(p->mv[1] - q->mv[0]) >= 4)
        bS = 1;

    bS = 0;
*/
    if (slice->slice_type == SI_slice || slice->slice_type == SP_slice ||
        MbQ->is_intra_block || MbP->is_intra_block) {
        // Set strength to either 3 or 4 regardless of pixel position
        StrValue = (edge == 0) ? 4 : 3;
        memset(Strength, (byte) StrValue, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && slice->slice_type == P_slice && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, BLOCK_SIZE * sizeof(byte));
        return;
    }

    for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
        int blkP = idx + ((xQ & 15) / 4);
        int blkQ = idx + (edge);

        if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
            (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P16x8))
            StrValue = 0;
        else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
            BlockPos mb = PicPos[MbQ->mbAddrX];

            int blk_x  = (mb.x * 4) + (blkQ & 3);
            int blk_y  = (mb.y * 4) + (blkQ / 4);
            int blk_x2 = neighbor->mb.x * 4 + (xQ & 15) / 4;
            int blk_y2 = neighbor->mb.y * 4 + idx       / 4;

            StrValue = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, p);
        }
        Strength[idx >> 2] = StrValue;
    }
}

static void get_strength_hor(mb_t *MbQ, int edge)
{  
    uint8_t*  Strength = MbQ->strength_hor[edge];
    int       StrValue;
    slice_t*  slice = MbQ->p_Slice;
    BlockPos* PicPos = MbQ->p_Vid->PicPos;
    storable_picture* p = MbQ->p_Vid->dec_picture;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

    int yQ = (edge < BLOCK_SIZE ? (edge << 2) - 1: 0);
    mb_t* neighbor = get_non_aff_neighbor_luma(MbQ, 0, yQ);
    mb_t* MbP = edge ? MbQ : neighbor;

    if (slice->slice_type == SP_SLICE || slice->slice_type == SI_SLICE ||
        MbQ->is_intra_block || MbP->is_intra_block) {
        // Set strength to either 3 or 4 regardless of pixel position
        StrValue = (edge == 0 && !slice->field_pic_flag) ? 4 : 3;
        memset(Strength, (byte) StrValue, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && slice->slice_type == P_SLICE && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, BLOCK_SIZE * sizeof(byte));
        return;
    }

    for (int idx = 0; idx < BLOCK_SIZE; idx++) {
        int blkQ = (yQ +  1) + idx;
        int blkP = (yQ & 12) + idx;

        if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
            (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P8x16))
            StrValue = 0;
        else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
            BlockPos mb = PicPos[MbQ->mbAddrX];

            int blk_x  = (mb.x * 4) + (blkQ & 3);
            int blk_y  = (mb.y * 4) + (blkQ / 4);
            int blk_x2 = neighbor->mb.x * 4 + idx;
            int blk_y2 = neighbor->mb.y * 4 + (yQ & 15) / 4;

            StrValue = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, p);
        }
        *(int*)(Strength + idx) = StrValue;
    }
}

static void get_strength_ver_MBAff(mb_t *MbQ, int edge)
{
    uint8_t* Strength = MbQ->strength_ver[edge/4];
    int StrValue;

    PixelPos pixP;
    VideoParameters *p_Vid = MbQ->p_Vid;
    storable_picture* p = p_Vid->dec_picture;
    slice_t* slice = MbQ->p_Slice;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    if (p->slice_type == SP_SLICE || p->slice_type == SI_SLICE) {
        for (int idx = 0; idx < MB_BLOCK_SIZE; ++idx)
            Strength[idx] = (edge == 0) ? 4 : 3;
    } else {
        getAffNeighbour(MbQ, edge - 1, 0, mb_size, &pixP);

        mb_t* MbP = &p_Vid->mb_data[pixP.mb_addr];
        // Neighboring Frame MBs
        if (!MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) {
            if (MbQ->is_intra_block || MbP->is_intra_block) {
                StrValue = (edge == 0) ? 4 : 3;
                memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
            } else {
                short mb_x, mb_y;
                get_mb_block_pos_mbaff(p_Vid->PicPos, MbQ->mbAddrX, &mb_x, &mb_y);

                for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
                    int blkP = (pixP.y & ~3) + (pixP.x / 4);
                    int blkQ = (idx    & ~3) + (edge   / 4);

                    if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
                        (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
                        StrValue = 2;
                    else if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P16x8))
                        StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
                    else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
                        int blk_x  = (mb_x * 4) + (blkQ & 3);
                        int blk_y  = (mb_y * 4) + (blkQ / 4);
                        int blk_x2 = (pixP.pos_x / 4);
                        int blk_y2 = (pixP.pos_y / 4);

                        StrValue = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, p);
                    }

                    *(int*)(Strength+idx) = StrValue * 0x01010101;

                    pixP.y     += 4;
                    pixP.pos_y += 4;
                }
            }
        } else {
            for (int idx = 0; idx < MB_BLOCK_SIZE; ++idx) {
                getAffNeighbour(MbQ, edge - 1, idx, mb_size, &pixP);

                int blkP = (pixP.y & ~3) + (pixP.x / 4);
                int blkQ = (idx    & ~3) + (edge   / 4);

                mb_t* MbP = &p_Vid->mb_data[pixP.mb_addr];
                MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag); 

                if (MbQ->is_intra_block || MbP->is_intra_block)
                    Strength[idx] = (edge == 0) ? 4 : 3;
                else if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
                         (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
                    Strength[idx] = 2;
                else if (MbQ->mixedModeEdgeFlag)
                    Strength[idx] = 1;
                else {
                    short mb_x, mb_y;
                    get_mb_block_pos_mbaff(p_Vid->PicPos, MbQ->mbAddrX, &mb_x, &mb_y);

                    int blk_x  = (mb_x * 4) + (blkQ & 3);
                    int blk_y  = (mb_y * 4) + (blkQ / 4);
                    int blk_x2 = (pixP.pos_x / 4);
                    int blk_y2 = (pixP.pos_y / 4);

                    Strength[idx] = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, p);
                }
            }
        }
    }
}

static void get_strength_hor_MBAff(mb_t *MbQ, int edge)
{
    uint8_t* Strength = MbQ->strength_hor[edge/4];
    int StrValue;
    int yQ = (edge < MB_BLOCK_SIZE ? edge : 1);

    PixelPos pixP;
    VideoParameters* p_Vid = MbQ->p_Vid;
    storable_picture* p = p_Vid->dec_picture;
    slice_t* slice = MbQ->p_Slice;
    int mvlimit = (slice->field_pic_flag || (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    if (p->slice_type == SI_slice || p->slice_type == SP_slice) {
        for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
            int xQ = idx;
            getAffNeighbour(MbQ, xQ, yQ - 1, mb_size, &pixP);

            mb_t* MbP = &p_Vid->mb_data[pixP.mb_addr];
            MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag);

            StrValue = (edge == 0) && (!MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag) ? 4 : 3;
      
            *(int*)(Strength+idx) = StrValue * 0x01010101;
        }
    } else {
        getAffNeighbour(MbQ, 0, yQ - 1, mb_size, &pixP);

        mb_t* MbP = &p_Vid->mb_data[pixP.mb_addr];
        MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag); 

        // Set intra mode deblocking
        if (MbQ->is_intra_block || MbP->is_intra_block) {
            StrValue = (edge == 0) && (!MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag) ? 4 : 3;
            memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
        } else {
            for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
                int xQ = idx;    
                getAffNeighbour(MbQ, xQ, yQ - 1, mb_size, &pixP);

                int blkP = (pixP.y & ~3) + (pixP.x / 4);
                int blkQ = (yQ     & ~3) + (xQ     / 4);

                if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
                    (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
                    StrValue = 2;
                else if (MbQ->mixedModeEdgeFlag)
                    StrValue = 1;
                else {
                    short mb_x, mb_y;
                    get_mb_block_pos_mbaff(p_Vid->PicPos, MbQ->mbAddrX, &mb_x, &mb_y);

                    int blk_x  = (mb_x * 4) + (blkQ & 3);
                    int blk_y  = (mb_y * 4) + (blkQ / 4);
                    int blk_x2 = (pixP.pos_x / 4);
                    int blk_y2 = (pixP.pos_y / 4);

                    StrValue = bs_compare_mvs(blk_x, blk_y, blk_x2, blk_y2, mvlimit, p);
                }

                *(int*)(Strength+idx) = StrValue * 0x01010101;
            }
        }
    }
}


/*********************************************************************************************************/

// NOTE: In principle, the alpha and beta tables are calculated with the formulas below
//       Alpha( qp ) = 0.8 * (2^(qp/6)  -  1)
//       Beta ( qp ) = 0.5 * qp  -  7

// The tables actually used have been "hand optimized" though (by Anthony Joch). So, the
// table values might be a little different to formula-generated values. Also, the first
// few values of both tables is set to zero to force the filter off at low qp’s

// Table 8-16 Derivation of offset dependent threshold variables a' and b' from indexA and indexB
static const byte ALPHA_TABLE[52] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   4,   4,   5,   6,   7,   8,   9,  10,  12,  13,
     15,  17,  20,  22,  25,  28,  32,  36,  40,  45,  50,  56,  63,
     71,  80,  90, 101, 113, 127, 144, 162, 182, 203, 226, 255, 255
};

static const byte BETA_TABLE[52] = {
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
      0,   0,   0,   2,   2,   2,   3,   3,   3,   3,   4,   4,   4,
      6,   6,   7,   7,   8,   8,   9,   9,  10,  10,  11,  11,  12,
     12,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18
};

static const byte TABLE_TCO[52][5] = {
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},
    { 0, 0, 0, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 1, 1, 1, 1},
    { 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 2, 2},
    { 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 2, 3, 3},
    { 0, 1, 2, 3, 3},{ 0, 2, 2, 3, 3},{ 0, 2, 2, 4, 4},{ 0, 2, 3, 4, 4},
    { 0, 2, 3, 4, 4},{ 0, 3, 3, 5, 5},{ 0, 3, 4, 6, 6},{ 0, 3, 4, 6, 6},
    { 0, 4, 5, 7, 7},{ 0, 4, 5, 8, 8},{ 0, 4, 6, 9, 9},{ 0, 5, 7,10,10},
    { 0, 6, 8,11,11},{ 0, 6, 8,13,13},{ 0, 7,10,14,14},{ 0, 8,11,16,16},
    { 0, 9,12,18,18},{ 0,10,13,20,20},{ 0,11,15,23,23},{ 0,13,17,25,25}
};

static const int pelnum_cr[2][4] =  {{0,8,16,16}, {0,8, 8,16}};  //[dir:0=vert, 1=hor.][yuv_format]


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

void deblock_t::deblock_normal(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag, int chromaEdgeFlag, int BitDepthY, int BitDepthC, int indexA)
{
#define p(i) (pixP[- (i) * widthP])
#define q(i) (pixQ[  (i) * widthQ])
    int  BitDepth = chromaEdgeFlag == 0 ? BitDepthY : BitDepthC;
    bool filterSamplesFlag = bS != 0 && abs(p(0) - q(0)) < alpha
                                     && abs(p(1) - p(0)) < beta
                                     && abs(q(1) - q(0)) < beta;

    if (filterSamplesFlag && bS < 4) {
        int tc0, tc, delta;
        int ap, aq;
        int p0, p1;
        int q0, q1;

        if (chromaEdgeFlag == 0)
            tc0 = TABLE_TCO[indexA][bS] * (1 << (BitDepthY - 8));
        else
            tc0 = TABLE_TCO[indexA][bS] * (1 << (BitDepthC - 8));

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

/*
    verticalEdgeFlag;
    chromaEdgeFlag;
    chromaStyleFilteringFlag = chromaEdgeFlag && (ChromaArrayType != 3)

    qPav = (qPp + qPq + 1) >> 1;

    indexA = clip3(0, 51, qPav + filterOffsetA);
    indexB = clip3(0, 51, qPav + filterOffsetB);

    if (chromaEdgeFlag == 0) {
        alpha = TABLE_ALPHA[indexA] * (1 << (BitDepthY - 8));
        beta  = TABLE_BETA [indexB] * (1 << (BitDepthY - 8));
    } else {
        alpha = TABLE_ALPHA[indexA] * (1 << (BitDepthC - 8));
        beta  = TABLE_BETA [indexB] * (1 << (BitDepthC - 8));
    }

    filterSamplesFlag = bS != 0 && abs(p[0] - q[0]) < alpha
                                && abs(p[1] - p[0]) < beta
                                && abs(q[1] - q[0]) < beta;
*/

void deblock_t::edge_loop(bool verticalEdgeFlag, bool chromaEdgeFlag, bool chromaStyleFilteringFlag,
                          mb_t *MbQ, byte *Strength, ColorPlane pl, int edge)
{
    VideoParameters *p_Vid = MbQ->p_Vid;
    sps_t *sps = p_Vid->active_sps;

    imgpel** Img = pl ? p_Vid->dec_picture->imgUV[pl-1] : p_Vid->dec_picture->imgY;
    int width = chromaEdgeFlag == 0 ? p_Vid->dec_picture->iLumaStride : p_Vid->dec_picture->iChromaStride;
    int PelNum = pl ? pelnum_cr[verticalEdgeFlag ? 0 : 1][sps->chroma_format_idc] : MB_BLOCK_SIZE;
    int bitdepth_scale = 1 << (pl ? sps->bit_depth_chroma_minus8 : sps->bit_depth_luma_minus8);
    PixelPos pixP, pixQ;

    int yQ = (edge < MB_BLOCK_SIZE ? edge : 1);

    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };
    int *mb_size = mb_size_xy[chromaEdgeFlag == 0 ? IS_LUMA : IS_CHROMA];

    if (MbQ || MbQ->p_Slice->disable_deblocking_filter_idc == 0) {
        for (int pel = 0; pel < PelNum; ++pel) {
            if (MbQ->p_Slice->MbaffFrameFlag) {
                if (verticalEdgeFlag) {
                    getAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);
                    getAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
                } else {
                    getAffNeighbour(MbQ, pel, yQ - 1, mb_size, &pixP);
                    getAffNeighbour(MbQ, pel, yQ    , mb_size, &pixQ);
                }
            } else {
                if (verticalEdgeFlag) {
                    getNonAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);
                    getNonAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
                } else {
                    getNonAffNeighbour(MbQ, pel, yQ - 1, mb_size, &pixP);
                    getNonAffNeighbour(MbQ, pel, yQ    , mb_size, &pixQ);
                }
            }

            mb_t *MbP = &p_Vid->mb_data[pixP.mb_addr];
            int StrengthIdx;
            if (MbQ->p_Slice->MbaffFrameFlag)
                StrengthIdx = (PelNum == 8) ? ((MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) ? pel << 1 : ((pel >> 1) << 2) + (pel & 0x01)) : pel;
            else
                StrengthIdx = PelNum == 8 ? pel >> 1 : pel >> 2;
            int bS = Strength[StrengthIdx];

            if (bS != 0) {
                imgpel *SrcPtrP = &Img[pixP.pos_y][pixP.pos_x];
                imgpel *SrcPtrQ = &Img[pixQ.pos_y][pixQ.pos_x];
                int incP, incQ;

                if (verticalEdgeFlag) {
                    incQ = 1;
                    incP = 1;
                } else if (MbQ->p_Slice->MbaffFrameFlag) {
                    incQ = ((MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag) ? 2 * width : width);
                    incP = ((MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) ? 2 * width : width);
                } else {
                    incP = width;
                    incQ = width;
                }

                // Average QP of the two blocks
                int QP = pl? ((MbP->QpC[pl-1] + MbQ->QpC[pl-1] + 1) >> 1)
                            : (MbP->QpY + MbQ->QpY + 1) >> 1;
                int indexA = clip3(0, MAX_QP, QP + MbQ->p_Slice->FilterOffsetA);
                int indexB = clip3(0, MAX_QP, QP + MbQ->p_Slice->FilterOffsetB);
                int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
                int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

                if (bS == 4)
                    this->deblock_strong(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag);
                else if (bS > 0)
                    this->deblock_normal(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag, pl, sps->BitDepthY, sps->BitDepthC, indexA);
            }
        }
    }
}


void deblock_t::edge_loop_luma_ver(ColorPlane pl, uint8_t* Strength, mb_t *MbQ, int edge)
{
    this->edge_loop(1, 0, 0, MbQ, Strength, pl, edge);
}

void deblock_t::edge_loop_luma_hor(ColorPlane pl, uint8_t* Strength, mb_t *MbQ, int edge)
{
    this->edge_loop(0, 0, 0, MbQ, Strength, pl, edge);
}

void deblock_t::edge_loop_chroma_ver(ColorPlane pl, uint8_t* Strength, mb_t *MbQ, int edge)
{
    this->edge_loop(1, 1, 1, MbQ, Strength, pl, edge);
}

void deblock_t::edge_loop_chroma_hor(ColorPlane pl, uint8_t* Strength, mb_t *MbQ, int edge)
{
    this->edge_loop(0, 1, 1, MbQ, Strength, pl, edge);
}


static const char chroma_edge[2][4][2] = { //[dir][edge][yuv_format]
    {{  0,  0 }, { -4, -4 }, {  4,  4 }, { -4, -4 }},
    {{  0,  0 }, { -4,  4 }, {  4,  8 }, { -4, 12 }}
};

void deblock_t::deblock_mb(mb_t* mb)
{
    mb_t* MbQ      = mb;
    slice_t* slice = MbQ->p_Slice;
    sps_t*   sps   = slice->active_sps;
    storable_picture* p = slice->p_Vid->dec_picture;

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

        int filterNon8x8LumaEdgesFlag[4] = {1,1,1,1};


//      filterLeftMbEdgeFlag = (MbaffFrameFlag == 0 && CurrMbAddr % PicWidthInMbs == 0) ||
//                             (MbaffFrameFlag == 1 && (CurrMbAddr >> 1) % PicWidthInMbs == 0) ||
//                             (slice->disable_deblocking_filter_idc == 1) ||
//                             (slice->disable_deblocking_filter_idc == 2 && mbAddrA == NULL) ? 0 : 1

//      filterTopMbEdgeFlag = (MbaffFrameFlag == 0 && CurrMbAddr < PicWidthInMbs) ||
//                            (MbaffFrameFlag == 1 && (CurrMbAddr >> 1) < PicWidthInMbs && mb->field) ||
//                            (MbaffFrameFlag == 1 && (CurrMbAddr >> 1) < PicWidthInMbs && mb->frame && CurrMbAddr % 2 == 0) ||
//                            (slice->disable_deblocking_filter_idc == 1) ||
//                            (slice->disable_deblocking_filter_idc == 2 && mbAddrB == NULL) ? 0 : 1

        filterNon8x8LumaEdgesFlag[1] =
        filterNon8x8LumaEdgesFlag[3] = !(MbQ->transform_size_8x8_flag);

        if (slice->MbaffFrameFlag && MbQ->mb_field_decoding_flag&& mb_y == MB_BLOCK_SIZE)
            filterTopMbEdgeFlag = 0;

        if (slice->disable_deblocking_filter_idc == 2) {
            // don't filter at slice boundaries
            filterLeftMbEdgeFlag = MbQ->mbAvailA;
            filterTopMbEdgeFlag  = MbQ->mbAvailB;
            if (slice->MbaffFrameFlag && !MbQ->mb_field_decoding_flag && (MbQ->mbAddrX & 0x01))
                filterTopMbEdgeFlag = 1;
        }

        if (slice->MbaffFrameFlag)
            CheckAvailabilityOfNeighborsMBAFF(MbQ);

        // Vertical deblocking
        for (int edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            if ((MbQ->CodedBlockPatternLuma == 0 && MbQ->CodedBlockPatternChroma == 0) &&
                (slice->slice_type == P_slice || slice->slice_type == B_slice)) {
                if (filterNon8x8LumaEdgesFlag[edge] == 0 && sps->chroma_format_idc != YUV444)
                    continue;
                else if (edge > 0) {
                    if (((MbQ->mb_type == PSKIP && slice->slice_type == P_slice) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P16x8)))
                        continue;
                    else if ((edge & 0x01) && ((MbQ->mb_type == P8x16) || (slice->slice_type == B_slice && MbQ->mb_type == BSKIP_DIRECT && sps->direct_8x8_inference_flag)))
                        continue;
                }
            }

            if (edge || filterLeftMbEdgeFlag) {
                uint8_t* Strength = MbQ->strength_ver[edge];
                if (slice->MbaffFrameFlag) // Strength for 4 blks in 1 stripe
                    get_strength_ver_MBAff(MbQ, edge * 4);
                else
                    get_strength_ver(MbQ, edge);

                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) { // only if one of the 16 Strength bytes is != 0
                //if ((*((int *) Strength))) { // only if one of the 16 Strength bytes is != 0
                //if (p_Strength64[0] || p_Strength64[1]) { // only if one of the 16 Strength bytes is != 0
                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        this->edge_loop_luma_ver(PLANE_Y, Strength, MbQ, edge * 4);
                        if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
                            this->edge_loop_luma_ver(PLANE_U, Strength, MbQ, edge * 4);
                            this->edge_loop_luma_ver(PLANE_V, Strength, MbQ, edge * 4);
                        }
                    }
                    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422) {
                        int edge_cr = chroma_edge[0][edge][sps->chroma_format_idc == YUV422];
                        if (edge_cr >= 0) {
                            this->edge_loop_chroma_ver(PLANE_U, Strength, MbQ, edge_cr);
                            this->edge_loop_chroma_ver(PLANE_V, Strength, MbQ, edge_cr);
                        }
                    }
                }        
            }
        }

        // horizontal deblocking  
        for (int edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            if ((MbQ->CodedBlockPatternLuma == 0 && MbQ->CodedBlockPatternChroma == 0) && 
                (slice->slice_type == P_slice || slice->slice_type == B_slice)) {
                if (filterNon8x8LumaEdgesFlag[edge] == 0 && sps->chroma_format_idc == YUV420)
                    continue;
                else if (edge > 0) {
                    if (((MbQ->mb_type == PSKIP && slice->slice_type == P_slice) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P8x16)))
                        continue;
                    else if ((edge & 0x01) && ((MbQ->mb_type == P16x8) || (slice->slice_type == B_slice && MbQ->mb_type == BSKIP_DIRECT && sps->direct_8x8_inference_flag)))
                        continue;
                }
            }

            if (edge || filterTopMbEdgeFlag) {
                uint8_t* Strength = MbQ->strength_hor[edge];
                if (slice->MbaffFrameFlag) // Strength for 4 blks in 1 stripe
                    get_strength_hor_MBAff(MbQ, edge * 4);
                else
                    get_strength_hor(MbQ, edge);

                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) { // only if one of the 16 Strength bytes is != 0
                //if (p_Strength64[0] || p_Strength64[1]) { // only if one of the 16 Strength bytes is != 0
                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        this->edge_loop_luma_hor(PLANE_Y, Strength, MbQ, edge * 4);
                        if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
                            this->edge_loop_luma_hor(PLANE_U, Strength, MbQ, edge * 4);
                            this->edge_loop_luma_hor(PLANE_V, Strength, MbQ, edge * 4);
                        }
                    }
                    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422) {
                        int edge_cr = chroma_edge[1][edge][p->chroma_format_idc == YUV422];
                        if (edge_cr >= 0) {
                            this->edge_loop_chroma_hor(PLANE_U, Strength, MbQ, edge_cr);
                            this->edge_loop_chroma_hor(PLANE_V, Strength, MbQ, edge_cr);
                        }
                    }
                }

                if (!edge && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag) {
                    // this is the extra horizontal edge between a frame macroblock pair and a field above it
                    MbQ->DeblockCall = 2;

                    uint8_t* Strength = MbQ->strength_hor[4];
                    if (slice->MbaffFrameFlag)
                        get_strength_hor_MBAff(MbQ, MB_BLOCK_SIZE); // Strength for 4 blks in 1 stripe
                    else
                        get_strength_hor(MbQ, 4); // Strength for 4 blks in 1 stripe

                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        this->edge_loop_luma_hor(PLANE_Y, Strength, MbQ, MB_BLOCK_SIZE);
                        if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
                            this->edge_loop_luma_hor(PLANE_U, Strength, MbQ, MB_BLOCK_SIZE);
                            this->edge_loop_luma_hor(PLANE_V, Strength, MbQ, MB_BLOCK_SIZE);
                        }
                    }
                    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422) {
                        int edge_cr = chroma_edge[1][edge][sps->chroma_format_idc == YUV422];
                        if (edge_cr >= 0) {
                            this->edge_loop_chroma_hor(PLANE_U, Strength, MbQ, MB_BLOCK_SIZE);
                            this->edge_loop_chroma_hor(PLANE_V, Strength, MbQ, MB_BLOCK_SIZE);
                        }
                    }

                    MbQ->DeblockCall = 1;
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
