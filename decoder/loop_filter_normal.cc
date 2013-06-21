
/*!
 *************************************************************************************
 * \file loop_filter_normal.c
 *
 * \brief
 *    Loop filter to reduce blocking artifacts on a macroblock level (normal).
 *    The filter strength is QP dependent.
 *
 * \author
 *    Contributors:
 *    - Peter List       Peter.List@t-systems.de:  Original code                                 (13-Aug-2001)
 *    - Jani Lainema     Jani.Lainema@nokia.com:   Some bug fixing, removal of recursiveness     (16-Aug-2001)
 *    - Peter List       Peter.List@t-systems.de:  inplace filtering and various simplifications (10-Jan-2002)
 *    - Anthony Joch     anthony@ubvideo.com:      Simplified switching between filters and
 *                                                 non-recursive default filter.                 (08-Jul-2002)
 *    - Cristina Gomila  cristina.gomila@thomson.net: Simplification of the chroma deblocking
 *                                                    from JVT-E089                              (21-Nov-2002)
 *    - Alexis Michael Tourapis atour@dolby.com:   Speed/Architecture improvements               (08-Feb-2007)
 *************************************************************************************
 */

#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "image.h"
#include "neighbour.h"
#include "loop_filter.h"
#include "loop_filter_common.h"

static void get_strength_ver    (Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p);
static void get_strength_hor    (Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p);
static void edge_loop_luma_ver  (ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p);
static void edge_loop_luma_hor  (ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p);
static void edge_loop_chroma_ver(imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p);
static void edge_loop_chroma_hor(imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p);


void set_loop_filter_functions_normal(VideoParameters *p_Vid)
{
  p_Vid->GetStrengthVer    = get_strength_ver;
  p_Vid->GetStrengthHor    = get_strength_hor;
  p_Vid->EdgeLoopLumaVer   = edge_loop_luma_ver;
  p_Vid->EdgeLoopLumaHor   = edge_loop_luma_hor;
  p_Vid->EdgeLoopChromaVer = edge_loop_chroma_ver;
  p_Vid->EdgeLoopChromaHor = edge_loop_chroma_hor;
}


static Macroblock* get_non_aff_neighbor_luma(Macroblock *mb, int xN, int yN)
{
    if (xN < 0)
        return mb->mbleft;
    else if (yN < 0)
        return mb->mbup;
    else
        return mb;
}

static Macroblock* get_non_aff_neighbor_chroma(Macroblock *mb, int xN, int yN, int block_width, int block_height)
{
    if (xN < 0) {
        if (yN < block_height)
            return mb->mbleft;
        else
            return NULL;
    } else if (xN < block_width) {
        if (yN < 0)
            return mb->mbup;
        else if (yN < block_height)
            return mb;
        else
            return NULL;
    } else
        return NULL;
}

#define get_x_luma(x) (x & 15)
#define get_y_luma(y) (y & 15)
#define get_pos_x_luma(mb,x) (mb->pix_x + (x & 15))
#define get_pos_y_luma(mb,y) (mb->pix_y + (y & 15))
#define get_pos_x_chroma(mb,x,max) (mb->pix_c_x + (x & max))
#define get_pos_y_chroma(mb,y,max) (mb->pix_c_y + (y & max))

  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */
static void get_strength_ver(Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p)
{
    byte *Strength = MbQ->strength_ver[edge];
    Slice *currSlice = MbQ->p_Slice;
    int     StrValue;
    BlockPos *PicPos = MbQ->p_Vid->PicPos;

    int xQ = (edge << 2) - 1;
    Macroblock *neighbor = get_non_aff_neighbor_luma(MbQ, xQ, 0);
    Macroblock *MbP = (edge) ? MbQ : neighbor;

    if (currSlice->slice_type == SP_SLICE || currSlice->slice_type == SI_SLICE ||
        MbQ->is_intra_block || (MbP->is_intra_block && edge == 0)) {
        // Set strength to either 3 or 4 regardless of pixel position
        StrValue = (edge == 0) ? 4 : 3;
        memset(Strength, (byte) StrValue, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && currSlice->slice_type == P_SLICE && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P16x8)) {
        int blkP, blkQ, idx;
        for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
            blkQ = idx + (edge);
            blkP = idx + (get_x_luma(xQ) >> 2);
            if ((MbQ->s_cbp[0].blk & (i64_power2(blkQ) | i64_power2(blkP))) != 0)
                StrValue = 2;
            else
                StrValue = 0; // if internal edge of certain types, then we already know StrValue should be 0
            Strength[idx >> 2] = StrValue;
        }
        return;
    }

    int      blkP, blkQ, idx;
    BlockPos mb = PicPos[ MbQ->mbAddrX ];
    mb.x <<= BLOCK_SHIFT;
    mb.y <<= BLOCK_SHIFT;

    for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
        blkQ = idx  + (edge);
        blkP = idx  + (get_x_luma(xQ) >> 2);
        if ((MbQ->s_cbp[0].blk & i64_power2(blkQ)) != 0 ||
            (MbP->s_cbp[0].blk & i64_power2(blkP)) != 0)
            StrValue = 2;
        else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
            int blk_y  = mb.y + (blkQ >> 2);
            int blk_x  = mb.x + (blkQ  & 3);
            int blk_y2 = (short)(get_pos_y_luma(neighbor,  0) + idx) >> 2;
            int blk_x2 = (short)(get_pos_x_luma(neighbor, xQ)      ) >> 2;
            PicMotionParams *mv_info_p = &p->mv_info[blk_y ][blk_x ];            
            PicMotionParams *mv_info_q = &p->mv_info[blk_y2][blk_x2];            
            StorablePicturePtr ref_p0 = mv_info_p->ref_pic[LIST_0];
            StorablePicturePtr ref_q0 = mv_info_q->ref_pic[LIST_0];            
            StorablePicturePtr ref_p1 = mv_info_p->ref_pic[LIST_1];
            StorablePicturePtr ref_q1 = mv_info_q->ref_pic[LIST_1];

            if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) || ((ref_p0==ref_q1) && (ref_p1==ref_q0))) {
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
                            && (
                            compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                            compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)
                        ));
                }
            } else
                StrValue = 1;
        }
        Strength[idx >> 2] = StrValue;
    }
}

  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */
static void get_strength_hor(Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p)
{  
    byte  *Strength = MbQ->strength_hor[edge];
    int    StrValue;
    Slice *currSlice = MbQ->p_Slice;
    BlockPos *PicPos = MbQ->p_Vid->PicPos;

    int yQ = (edge < BLOCK_SIZE ? (edge << 2) - 1: 0);
    Macroblock *neighbor = get_non_aff_neighbor_luma(MbQ, 0, yQ);
    Macroblock *MbP = (edge) ? MbQ : neighbor;

    if (currSlice->slice_type == SP_SLICE || currSlice->slice_type == SI_SLICE ||
        MbQ->is_intra_block || (MbP->is_intra_block && edge == 0)) {
        // Set strength to either 3 or 4 regardless of pixel position
        StrValue = (edge == 0 && p->structure == FRAME) ? 4 : 3;
        memset(Strength, (byte) StrValue, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && currSlice->slice_type == P_SLICE && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P8x16)) {
        int blkP, blkQ, idx;
        for (idx = 0; idx < BLOCK_SIZE; idx++) {
            blkQ = (yQ + 1) + idx;
            blkP = (get_y_luma(yQ) & 0xFFFC) + idx;
            if ((MbQ->s_cbp[0].blk & (i64_power2(blkQ) | i64_power2(blkP))) != 0)
                StrValue = 2;
            else
                StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
            *(int*)(Strength + idx) = StrValue;
        }
        return;
    }

    int      blkP, blkQ, idx;
    BlockPos mb = PicPos[ MbQ->mbAddrX ];
    mb.x <<= 2;
    mb.y <<= 2;

    for (idx = 0; idx < BLOCK_SIZE; idx++) {
        blkQ = (yQ + 1) + idx;
        blkP = (get_y_luma(yQ) & 0xFFFC) + idx;

        if (((MbQ->s_cbp[0].blk & i64_power2(blkQ)) != 0) || ((MbP->s_cbp[0].blk & i64_power2(blkP)) != 0))
            StrValue = 2;
        else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
            int blk_y  = mb.y + (blkQ >> 2);
            int blk_x  = mb.x + (blkQ  & 3);
            int blk_y2 = get_pos_y_luma(neighbor,yQ) >> 2;
            int blk_x2 = ((short)(get_pos_x_luma(neighbor,0)) >> 2) + idx;
            PicMotionParams *mv_info_p = &p->mv_info[blk_y ][blk_x ];
            PicMotionParams *mv_info_q = &p->mv_info[blk_y2][blk_x2];
            StorablePicturePtr ref_p0 = mv_info_p->ref_pic[LIST_0];
            StorablePicturePtr ref_q0 = mv_info_q->ref_pic[LIST_0];
            StorablePicturePtr ref_p1 = mv_info_p->ref_pic[LIST_1];
            StorablePicturePtr ref_q1 = mv_info_q->ref_pic[LIST_1];            

            if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1)) || ((ref_p0==ref_q1) && (ref_p1==ref_q0))) {
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
                        && (
                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) |
                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)
                    ));
                }
            } else
                StrValue = 1;
        }
        *(int*)(Strength + idx) = StrValue;
    }
}

static inline int Abs(int x)
{
  static const int INT_BITS = (sizeof(int) * CHAR_BIT) - 1;
  int y = x >> INT_BITS;
  return (x ^ y) - y;
}

static inline int Clip3(int low, int high, int x)
{
  x = imax(x, low);
  x = imin(x, high);

  return x;
}

static const byte TABLE_TCO[52][5] = {
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 1, 1, 1, 1},
    { 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 2, 3, 3},
    { 0, 1, 2, 3, 3},{ 0, 2, 2, 3, 3},{ 0, 2, 2, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 3, 3, 5, 5},{ 0, 3, 4, 6, 6},{ 0, 3, 4, 6, 6},
    { 0, 4, 5, 7, 7},{ 0, 4, 5, 8, 8},{ 0, 4, 6, 9, 9},{ 0, 5, 7,10,10},{ 0, 6, 8,11,11},{ 0, 6, 8,13,13},{ 0, 7,10,14,14},{ 0, 8,11,16,16},
    { 0, 9,12,18,18},{ 0,10,13,20,20},{ 0,11,15,23,23},{ 0,13,17,25,25}
};

/*!
 *****************************************************************************************
 * \brief
 *    Vertical Deblocking with Strength = 4
 *****************************************************************************************
 */
static void luma_deblock_strong(imgpel *pixP, imgpel *pixQ, int width, int alpha, int beta)
{
#define p(i) (pixP[- (i) * width])
#define q(i) (pixQ[  (i) * width])
    int  bS = 4;
    bool chromaStyleFilteringFlag = 0;
    bool filterSamplesFlag = bS != 0 && Abs(p(0) - q(0)) < alpha
                                     && Abs(p(1) - p(0)) < beta
                                     && Abs(q(1) - q(0)) < beta;

    if (filterSamplesFlag && bS == 4) {
        int ap = Abs(p(2) - p(0));
        int aq = Abs(q(2) - q(0));
        int p0, p1, p2;
        int q0, q1, q2;

        if (chromaStyleFilteringFlag == 0 && ap < beta && Abs(p(0) - q(0)) < (alpha >> 2) + 2) {
            p0 = (p(2) + 2 * p(1) + 2 * p(0) + 2 * q(0) + q(1) + 4) >> 3;
            p1 = (p(2) + p(1) + p(0) + q(0) + 2) >> 2;
            p2 = (2 * p(3) + 3 * p(2) + p(1) + p(0) + q(0) + 4) >> 3;
        } else {
            p0 = (2 * p(1) + p(0) + q(1) + 2) >> 2;
            p1 = p(1);
            p2 = p(2);
        }

        if (chromaStyleFilteringFlag == 0 && aq < beta && Abs(p(0) - q(0)) < (alpha >> 2) + 2) {
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

static void luma_deblock_normal(imgpel *pixP, imgpel *pixQ, int width, int alpha, int beta, int chromaEdgeFlag, int BitDepthY, int BitDepthC, int indexA, int bS)
{
#define p(i) (pixP[- (i) * width])
#define q(i) (pixQ[  (i) * width])
    int  BitDepth = chromaEdgeFlag == 0 ? BitDepthY : BitDepthC;
    bool chromaStyleFilteringFlag = 0;
    bool filterSamplesFlag = bS != 0 && Abs(p(0) - q(0)) < alpha
                                     && Abs(p(1) - p(0)) < beta
                                     && Abs(q(1) - q(0)) < beta;

    if (filterSamplesFlag && bS < 4) {
        int tc0, tc, delta;
        int ap, aq;
        int p0, p1;
        int q0, q1;

        if (chromaEdgeFlag == 0)
            tc0 = TABLE_TCO[indexA][bS] * (1 << (BitDepthY - 8));
        else
            tc0 = TABLE_TCO[indexA][bS] * (1 << (BitDepthC - 8));

        ap = Abs(p(2) - p(0));
        aq = Abs(q(2) - q(0));
        if (chromaStyleFilteringFlag == 0)
            tc = tc0 + (ap < beta ? 1 : 0) + (aq < beta ? 1 : 0);
        else
            tc = tc0 + 1;
        delta = Clip3(-tc, tc, ((((q(0) - p(0)) << 2) + (p(1) - q(1)) + 4) >> 3));

#define Clip1(x) (Clip3(0, (1 << BitDepth) - 1, x))
        p0 = Clip1(p(0) + delta);
        q0 = Clip1(q(0) - delta);
#undef Clip1

        if (chromaStyleFilteringFlag == 0 && ap < beta)
            p1 = p(1) + Clip3(-tc0, tc0, (p(2) + ((p(0) + q(0) + 1) >> 1) - (p(1) << 1)) >> 1);
        else
            p1 = p(1);
        if (chromaStyleFilteringFlag == 0 && aq < beta)
            q1 = q(1) + Clip3(-tc0, tc0, (q(2) + ((p(0) + q(0) + 1) >> 1) - (q(1) << 1)) >> 1);
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

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Frame or Field coded MBs 
 *****************************************************************************************
 */
static void edge_loop_luma_ver(ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p)
{
    VideoParameters *p_Vid = MbQ->p_Vid;
    Macroblock *MbP = get_non_aff_neighbor_luma(MbQ, edge - 1, 0);

    // chromaStyleFilteringFlag = chromaEdgeFlag && (ChromaArrayType != 3)

    if (MbP || MbQ->DFDisableIdc == 0) {
/*
        qPav = (qPp + qPq + 1) >> 1;

        indexA = Clip3(0, 51, qPav + filterOffsetA);
        indexB = Clip3(0, 51, qPav + filterOffsetB);

        if (chromaEdgeFlag == 0) {
            alpha = TABLE_ALPHA[indexA] * (1 << (BitDepthY - 8));
            beta  = TABLE_BETA [indexB] * (1 << (BitDepthY - 8));
        } else {
            alpha = TABLE_ALPHA[indexA] * (1 << (BitDepthC - 8));
            beta  = TABLE_BETA [indexB] * (1 << (BitDepthC - 8));
        }

        filterSamplesFlag = bS != 0 && Abs(p[0] - q[0]) < alpha
                                    && Abs(p[1] - p[0]) < beta
                                    && Abs(q[1] - q[0]) < beta;
*/
        int bitdepth_scale   = pl ? p_Vid->bitdepth_scale[IS_CHROMA] :
                                    p_Vid->bitdepth_scale[IS_LUMA];

        // Average QP of the two blocks
        int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) :
                      (MbP->qp + MbQ->qp + 1) >> 1;

        int indexA = iClip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
        int indexB = iClip3(0, MAX_QP, QP + MbQ->DFBetaOffset);

        int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
        int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

        if ((Alpha | Beta) != 0) {
            int pos_x1 = get_pos_x_luma(MbP, (edge - 1));
            imgpel **cur_img = &Img[get_pos_y_luma(MbP, 0)];
            int pel;

            for (pel = 0; pel < MB_BLOCK_SIZE; pel += 4) {
                for (int i = 0; i < BLOCK_SIZE; ++i) {
                    if (*Strength == 4)    // INTRA strong filtering
                        luma_deblock_strong(&cur_img[i][pos_x1], &cur_img[i][pos_x1+1], 1, Alpha, Beta);
                    else if (*Strength != 0) // normal filtering
                        luma_deblock_normal(&cur_img[i][pos_x1], &cur_img[i][pos_x1+1], 1, Alpha, Beta, pl, p_Vid->bitdepth_luma, p_Vid->bitdepth_chroma, indexA, *Strength);
                }
                cur_img += 4;
                Strength++;
            }
        }
    }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Frame or Field coded MBs 
 *****************************************************************************************
 */
static void edge_loop_luma_hor(ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p)
{
    VideoParameters *p_Vid = MbQ->p_Vid;
    int ypos = (edge < MB_BLOCK_SIZE ? edge - 1: 0);
    Macroblock *MbP = get_non_aff_neighbor_luma(MbQ, 0, ypos); 

    if (MbP || MbQ->DFDisableIdc == 0) {
        int bitdepth_scale   = pl ? p_Vid->bitdepth_scale[IS_CHROMA] : p_Vid->bitdepth_scale[IS_LUMA];

        // Average QP of the two blocks
        int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) : (MbP->qp + MbQ->qp + 1) >> 1;

        int indexA = iClip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
        int indexB = iClip3(0, MAX_QP, QP + MbQ->DFBetaOffset);

        int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
        int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

        if ((Alpha | Beta )!= 0)
        {
            int width = p->iLumaStride; //p->size_x;
            imgpel *imgP = &Img[get_pos_y_luma(MbP, ypos)][get_pos_x_luma(MbP, 0)];
            imgpel *imgQ = imgP + width;
            int pel;

            for (pel = 0; pel < BLOCK_SIZE; pel++) {
                for (int i = 0; i < BLOCK_SIZE; ++i) {
                    if (*Strength == 4)    // INTRA strong filtering
                        luma_deblock_strong(&imgP[i], &imgQ[i], width, Alpha, Beta);
                    else if (*Strength != 0) // normal filtering
                        luma_deblock_normal(&imgP[i], &imgQ[i], width, Alpha, Beta, pl, p_Vid->bitdepth_luma, p_Vid->bitdepth_chroma, indexA, *Strength);
                }
                imgP += 4;
                imgQ += 4;
                Strength ++;
            }
        }
    }
}


/*!
 *****************************************************************************************
 * \brief
 *    Filters chroma block edge for Frame or Field coded pictures
 *****************************************************************************************
 */
static void edge_loop_chroma_ver(imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p)
{
  VideoParameters *p_Vid = MbQ->p_Vid;  

  int block_width  = p_Vid->mb_cr_size_x;
  int block_height = p_Vid->mb_cr_size_y;
  int xQ = edge - 1;
  int yQ = 0;  

  Macroblock *MbP = get_non_aff_neighbor_chroma(MbQ,xQ,yQ,block_width,block_height); 

  if (MbP || (MbQ->DFDisableIdc == 0))
  {
    int      bitdepth_scale   = p_Vid->bitdepth_scale[IS_CHROMA];
    int      max_imgpel_value = p_Vid->max_pel_value_comp[uv + 1];

    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;

    // Average QP of the two blocks
    int QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    int Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta    = BETA_TABLE [indexB] * bitdepth_scale;

    if ((Alpha | Beta) != 0)
    {
      const int PelNum = pelnum_cr[0][p->chroma_format_idc];
      const     byte *ClipTab = CLIP_TAB[indexA];

      int pel;
      int pos_x1 = get_pos_x_chroma(MbP, xQ, (block_width - 1));
      imgpel **cur_img = &Img[get_pos_y_chroma(MbP,yQ, (block_height - 1))];

      for( pel = 0 ; pel < PelNum ; ++pel )
      {
        int Strng = Strength[(PelNum == 8) ? (pel >> 1) : (pel >> 2)];

        if( Strng != 0)
        {
          imgpel *SrcPtrP = *cur_img + pos_x1;
          imgpel *SrcPtrQ = SrcPtrP + 1;
          int edge_diff = *SrcPtrQ - *SrcPtrP;

          if ( iabs( edge_diff ) < Alpha ) 
          {
            imgpel R1  = *(SrcPtrQ + 1);
            if ( iabs(*SrcPtrQ - R1) < Beta )  
            {
              imgpel L1  = *(SrcPtrP - 1);
              if ( iabs(*SrcPtrP - L1) < Beta )
              {
                if( Strng == 4 )    // INTRA strong filtering
                {
                  *SrcPtrP = (imgpel) ( ((L1 << 1) + *SrcPtrP + R1 + 2) >> 2 );
                  *SrcPtrQ = (imgpel) ( ((R1 << 1) + *SrcPtrQ + L1 + 2) >> 2 );
                }
                else
                {
                  int tc0  = ClipTab[ Strng ] * bitdepth_scale + 1;
                  int dif = iClip3( -tc0, tc0, ( ((edge_diff) << 2) + (L1 - R1) + 4) >> 3 );

                  if (dif != 0)
                  {
                    *SrcPtrP = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrP + dif );
                    *SrcPtrQ = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrQ - dif );
                  }
                }
              }
            }
          }
        }
        cur_img++;
      }     
    }
  }
}


/*!
 *****************************************************************************************
 * \brief
 *    Filters chroma block edge for Frame or Field coded pictures
 *****************************************************************************************
 */
static void edge_loop_chroma_hor(imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p)
{
  VideoParameters *p_Vid = MbQ->p_Vid;  
  int block_width = p_Vid->mb_cr_size_x;
  int block_height = p_Vid->mb_cr_size_y;
  int xQ = 0;
  int yQ = (edge < 16 ? edge - 1: 0);

  Macroblock *MbP = get_non_aff_neighbor_chroma(MbQ,xQ,yQ,block_width,block_height);

  if (MbP || (MbQ->DFDisableIdc == 0))
  {
    int      bitdepth_scale   = p_Vid->bitdepth_scale[IS_CHROMA];
    int      max_imgpel_value = p_Vid->max_pel_value_comp[uv + 1];

    int AlphaC0Offset = MbQ->DFAlphaC0Offset;
    int BetaOffset = MbQ->DFBetaOffset;
    int width = p->iChromaStride; //p->size_x_cr;

    // Average QP of the two blocks
    int QP = (MbP->qpc[uv] + MbQ->qpc[uv] + 1) >> 1;

    int indexA = iClip3(0, MAX_QP, QP + AlphaC0Offset);
    int indexB = iClip3(0, MAX_QP, QP + BetaOffset);

    int Alpha   = ALPHA_TABLE[indexA] * bitdepth_scale;
    int Beta    = BETA_TABLE [indexB] * bitdepth_scale;

    if ((Alpha | Beta) != 0)
    {
      const int PelNum = pelnum_cr[1][p->chroma_format_idc];
      const     byte *ClipTab = CLIP_TAB[indexA];

      int pel;

      imgpel *imgP = &Img[get_pos_y_chroma(MbP,yQ, (block_height-1))][get_pos_x_chroma(MbP,xQ, (block_width - 1))];
      imgpel *imgQ = imgP + width ;

      for( pel = 0 ; pel < PelNum ; ++pel )
      {
        int Strng = Strength[(PelNum == 8) ? (pel >> 1) : (pel >> 2)];

        if( Strng != 0)
        {
          imgpel *SrcPtrP = imgP;
          imgpel *SrcPtrQ = imgQ;
          int edge_diff = *imgQ - *imgP;

          if ( iabs( edge_diff ) < Alpha ) 
          {
            imgpel R1  = *(SrcPtrQ + width);
            if ( iabs(*SrcPtrQ - R1) < Beta )  
            {
              imgpel L1  = *(SrcPtrP - width);
              if ( iabs(*SrcPtrP - L1) < Beta )
              {
                if( Strng == 4 )    // INTRA strong filtering
                {
                  *SrcPtrP = (imgpel) ( ((L1 << 1) + *SrcPtrP + R1 + 2) >> 2 );
                  *SrcPtrQ = (imgpel) ( ((R1 << 1) + *SrcPtrQ + L1 + 2) >> 2 );
                }
                else
                {
                  int tc0  = ClipTab[ Strng ] * bitdepth_scale + 1;
                  int dif = iClip3( -tc0, tc0, ( ((edge_diff) << 2) + (L1 - R1) + 4) >> 3 );

                  if (dif != 0)
                  {
                    *SrcPtrP = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrP + dif );
                    *SrcPtrQ = (imgpel) iClip1 ( max_imgpel_value, *SrcPtrQ - dif );
                  }
                }
              }
            }
          }
        }
        imgP++;
        imgQ++;
      }
    }
  }
}
