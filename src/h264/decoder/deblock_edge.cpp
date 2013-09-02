#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "image.h"
#include "neighbour.h"
#include "deblock.h"
#include "deblock_common.h"

#define MAX_QP          51

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
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},{ 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 0, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 0, 1, 1, 1},{ 0, 1, 1, 1, 1},
    { 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 1, 1},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 1, 2, 2},{ 0, 1, 2, 3, 3},
    { 0, 1, 2, 3, 3},{ 0, 2, 2, 3, 3},{ 0, 2, 2, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 2, 3, 4, 4},{ 0, 3, 3, 5, 5},{ 0, 3, 4, 6, 6},{ 0, 3, 4, 6, 6},
    { 0, 4, 5, 7, 7},{ 0, 4, 5, 8, 8},{ 0, 4, 6, 9, 9},{ 0, 5, 7,10,10},{ 0, 6, 8,11,11},{ 0, 6, 8,13,13},{ 0, 7,10,14,14},{ 0, 8,11,16,16},
    { 0, 9,12,18,18},{ 0,10,13,20,20},{ 0,11,15,23,23},{ 0,13,17,25,25}
};

static const int pelnum_cr[2][4] =  {{0,8,16,16}, {0,8, 8,16}};  //[dir:0=vert, 1=hor.][yuv_format]


static mb_t* get_non_aff_neighbor_luma(mb_t *mb, int xN, int yN)
{
    if (xN < 0)
        return mb->mbleft;
    else if (yN < 0)
        return mb->mbup;
    else
        return mb;
}

static mb_t* get_non_aff_neighbor_chroma(mb_t *mb, int xN, int yN, int block_width, int block_height)
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


/*!
 *****************************************************************************************
 * \brief
 *    Vertical Deblocking with Strength = 4
 *****************************************************************************************
 */
static void deblock_strong(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag)
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

static void deblock_normal(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag, int chromaEdgeFlag, int BitDepthY, int BitDepthC, int indexA)
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

static void deblock_mb(bool verticalEdgeFlag, bool chromaEdgeFlag, bool chromaStyleFilteringFlag,
                       mb_t *MbP, mb_t *MbQ, byte *Strength,
                       ColorPlane pl, imgpel **Img, int edge, StorablePicture *p)
{
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
    VideoParameters *p_Vid = MbQ->p_Vid;
    sps_t *sps = p_Vid->active_sps;
    int width = chromaEdgeFlag == 0 ? p->iLumaStride : p->iChromaStride;
    int PelNum = pl ? pelnum_cr[verticalEdgeFlag ? 0 : 1][p->chroma_format_idc] : MB_BLOCK_SIZE;
    int bitdepth_scale = 1 << (pl ? sps->bit_depth_chroma_minus8 : sps->bit_depth_luma_minus8);
    PixelPos pixP, pixQ;

    int yQ = (edge < 16 ? edge - 1: 0);

    int xQ = edge - 1;

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { mb_cr_size_x, mb_cr_size_y }
    };
    int *mb_size = mb_size_xy[chromaEdgeFlag == 0 ? IS_LUMA : IS_CHROMA];

    //if (MbQ->DFDisableIdc == 0) {
    if (MbP || MbQ->DFDisableIdc == 0) {

        int edge_x = verticalEdgeFlag ? edge : 0;
        int edge_y = verticalEdgeFlag ? 0 : edge;
        int step_x = verticalEdgeFlag ? 1 : width;
        int step_y = verticalEdgeFlag ? width : 1;

        if (verticalEdgeFlag) {
            xQ = edge - 1;
            yQ = 0;
        } else {
            xQ = 0;
            yQ = (edge < MB_BLOCK_SIZE ? edge - 1: 0);
        }

        imgpel *imgQ = chromaEdgeFlag == 0 ? &Img[MbQ->pix_y + edge_y][MbQ->pix_x + edge_x]
                                           : &Img[MbQ->pix_c_y + edge_y][MbQ->pix_c_x + edge_x];
        imgpel *imgP = &imgQ[-step_x];

        for (int pel = 0; pel < PelNum; pel++) {
            if (verticalEdgeFlag) {
                getNonAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);     
                getNonAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
            } else {
                getNonAffNeighbour(MbQ, pel, yQ - 1, mb_size, &pixP);     
                getNonAffNeighbour(MbQ, pel, yQ    , mb_size, &pixQ);
            }

            if (chromaEdgeFlag == 0)
                MbP = get_non_aff_neighbor_luma(MbQ, xQ, yQ);
            else
                MbP = get_non_aff_neighbor_chroma(MbQ, xQ, yQ, mb_cr_size_x, mb_cr_size_y);

            mb_t *MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            int bS = Strength[PelNum == 8 ? pel >> 1 : pel >> 2];

            //if (pixP.available && bS != 0) {
            if (bS != 0) {
                //imgpel *SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);
                //imgpel *SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
                imgpel *SrcPtrP = &imgP[pel*step_y];
                imgpel *SrcPtrQ = &imgQ[pel*step_y];
                int incP, incQ;

                if (verticalEdgeFlag) {
                    incP = 1;
                    incQ = 1;
                } else {
                    incP = width;
                    incQ = width;
                }

                // Average QP of the two blocks
                int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1) :
                              (MbP->qp + MbQ->qp + 1) >> 1;
                int indexA = clip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
                int indexB = clip3(0, MAX_QP, QP + MbQ->DFBetaOffset);
                int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
                int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

                if (bS == 4)
                    deblock_strong(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag);
                else if (bS != 0)
                    deblock_normal(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag, pl, sps->BitDepthY, sps->BitDepthC, indexA);
            }
        }
    }
}

static void deblock_mb_mbaff(bool verticalEdgeFlag, bool chromaEdgeFlag, bool chromaStyleFilteringFlag,
                       mb_t *MbP, mb_t *MbQ, byte *Strength,
                       ColorPlane pl, imgpel **Img, int edge, StorablePicture *p)
{
    VideoParameters *p_Vid = MbQ->p_Vid;
    sps_t *sps = p_Vid->active_sps;

    int width = chromaEdgeFlag == 0 ? p->iLumaStride : p->iChromaStride;
    int PelNum = pl ? pelnum_cr[verticalEdgeFlag ? 0 : 1][p->chroma_format_idc] : MB_BLOCK_SIZE;
    int bitdepth_scale = 1 << (pl ? sps->bit_depth_chroma_minus8 : sps->bit_depth_luma_minus8);
    PixelPos pixP, pixQ;

    int yQ = (edge < MB_BLOCK_SIZE ? edge : 1);

    int mb_cr_size_x = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV444 ? 16 : 8;
    int mb_cr_size_y = sps->chroma_format_idc == YUV400 ? 0 :
                       sps->chroma_format_idc == YUV420 ? 8 : 16;
    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { mb_cr_size_x, mb_cr_size_y }
    };
    int *mb_size = mb_size_xy[chromaEdgeFlag == 0 ? IS_LUMA : IS_CHROMA];

    if (MbQ->DFDisableIdc == 0) {
        for (int pel = 0; pel < PelNum; ++pel) {
            if (verticalEdgeFlag) {
                getAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);     
                getAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
            } else {
                getAffNeighbour(MbQ, pel, yQ - 1, mb_size, &pixP);     
                getAffNeighbour(MbQ, pel, yQ    , mb_size, &pixQ);
            }

            mb_t *MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            int StrengthIdx = (PelNum == 8) ? ((MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) ? pel << 1 : ((pel >> 1) << 2) + (pel & 0x01)) : pel;
            int bS = Strength[StrengthIdx];

            if (pixP.available && bS != 0) {
                imgpel *SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);
                imgpel *SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
                int incP, incQ;

                if (verticalEdgeFlag) {
                    incQ = 1;
                    incP = 1;
                } else {
                    incQ = ((MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag) ? 2 * width : width);
                    incP = ((MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) ? 2 * width : width);
                }

                // Average QP of the two blocks
                int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1)
                            : (MbP->qp + MbQ->qp + 1) >> 1;
                int indexA = clip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
                int indexB = clip3(0, MAX_QP, QP + MbQ->DFBetaOffset);
                int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
                int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

                if (bS == 4)
                    deblock_strong(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag);
                else if (bS > 0)
                    deblock_normal(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag, pl, sps->BitDepthY, sps->BitDepthC, indexA);
            }
        }
    }
}


static void edge_loop_luma_ver_MBAff(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, StorablePicture *p)
{
    deblock_mb_mbaff(1, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

static void edge_loop_luma_hor_MBAff(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, StorablePicture *p)
{
    deblock_mb_mbaff(0, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

static void edge_loop_chroma_ver_MBAff(imgpel** Img, byte *Strength, mb_t *MbQ, int edge, int uv, StorablePicture *p)
{
    deblock_mb_mbaff(1, 1, 1, MbQ, MbQ, Strength, (ColorPlane)(uv+1), Img, edge, p);
}

static void edge_loop_chroma_hor_MBAff(imgpel** Img, byte *Strength, mb_t *MbQ, int edge, int uv, StorablePicture *p)
{
    deblock_mb_mbaff(0, 1, 1, MbQ, MbQ, Strength, (ColorPlane)(uv+1), Img, edge, p);
}

static void edge_loop_luma_ver(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, StorablePicture *p)
{
    deblock_mb(1, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

static void edge_loop_luma_hor(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, StorablePicture *p)
{
    deblock_mb(0, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

static void edge_loop_chroma_ver(imgpel** Img, byte *Strength, mb_t *MbQ, int edge, int uv, StorablePicture *p)
{
    deblock_mb(1, 1, 1, MbQ, MbQ, Strength, (ColorPlane)(uv+1), Img, edge, p);
}

static void edge_loop_chroma_hor(imgpel** Img, byte *Strength, mb_t *MbQ, int edge, int uv, StorablePicture *p)
{
    deblock_mb(0, 1, 1, MbQ, MbQ, Strength, (ColorPlane)(uv+1), Img, edge, p);
}


void set_loop_filter_functions_normal(VideoParameters *p_Vid)
{
  p_Vid->GetStrengthVer    = get_strength_ver;
  p_Vid->GetStrengthHor    = get_strength_hor;
  p_Vid->EdgeLoopLumaVer   = edge_loop_luma_ver;
  p_Vid->EdgeLoopLumaHor   = edge_loop_luma_hor;
  p_Vid->EdgeLoopChromaVer = edge_loop_chroma_ver;
  p_Vid->EdgeLoopChromaHor = edge_loop_chroma_hor;
}

void set_loop_filter_functions_mbaff(VideoParameters *p_Vid)
{
  p_Vid->EdgeLoopLumaVer   = edge_loop_luma_ver_MBAff;
  p_Vid->EdgeLoopLumaHor   = edge_loop_luma_hor_MBAff;
  p_Vid->EdgeLoopChromaVer = edge_loop_chroma_ver_MBAff;
  p_Vid->EdgeLoopChromaHor = edge_loop_chroma_hor_MBAff;
}
