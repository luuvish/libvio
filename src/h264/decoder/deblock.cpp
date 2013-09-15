#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "image.h"
#include "neighbour.h"
#include "deblock.h"


deblock_t deblock;


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

static void deblock_mb(bool verticalEdgeFlag, bool chromaEdgeFlag, bool chromaStyleFilteringFlag,
                       mb_t *MbP, mb_t *MbQ, byte *Strength,
                       ColorPlane pl, imgpel **Img, int edge, storable_picture *p)
{
    VideoParameters *p_Vid = MbQ->p_Vid;
    sps_t *sps = p_Vid->active_sps;

    int width = chromaEdgeFlag == 0 ? p->iLumaStride : p->iChromaStride;
    int PelNum = pl ? pelnum_cr[verticalEdgeFlag ? 0 : 1][p->chroma_format_idc] : MB_BLOCK_SIZE;
    int bitdepth_scale = 1 << (pl ? sps->bit_depth_chroma_minus8 : sps->bit_depth_luma_minus8);
    PixelPos pixP, pixQ;

    int yQ = (edge < MB_BLOCK_SIZE ? edge : 1);

    int mb_size_xy[2][2] = {
        { MB_BLOCK_SIZE, MB_BLOCK_SIZE },
        { sps->MbWidthC, sps->MbHeightC }
    };
    int *mb_size = mb_size_xy[chromaEdgeFlag == 0 ? IS_LUMA : IS_CHROMA];

    if (MbP || MbQ->p_Slice->disable_deblocking_filter_idc == 0) {
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
                    deblock_strong(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag);
                else if (bS > 0)
                    deblock_normal(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag, pl, sps->BitDepthY, sps->BitDepthC, indexA);
            }
        }
    }
}


void deblock_t::edge_loop_luma_ver(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, storable_picture *p)
{
    deblock_mb(1, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

void deblock_t::edge_loop_luma_hor(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, storable_picture *p)
{
    deblock_mb(0, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

void deblock_t::edge_loop_chroma_ver(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, storable_picture *p)
{
    deblock_mb(1, 1, 1, MbQ, MbQ, Strength, pl, Img, edge, p);
}

void deblock_t::edge_loop_chroma_hor(ColorPlane pl, imgpel** Img, byte *Strength, mb_t *MbQ, int edge, storable_picture *p)
{
    deblock_mb(0, 1, 1, MbQ, MbQ, Strength, pl, Img, edge, p);
}


static const char chroma_edge[2][4][2] = { //[dir][edge][yuv_format]
    {{  0,  0 }, { -4, -4 }, {  4,  4 }, { -4, -4 }},
    {{  0,  0 }, { -4,  4 }, {  4,  8 }, { -4, 12 }}
};

void deblock_t::DeblockMb(VideoParameters *p_Vid, storable_picture *p, int MbQAddr)
{
    mb_t*    MbQ   = &p_Vid->mb_data[MbQAddr];
    slice_t* slice = MbQ->p_Slice;
    sps_t*   sps   = p_Vid->active_sps;

    // return, if filter is disabled
    if (slice->disable_deblocking_filter_idc != 1) {
        byte* Strength;
        byte  pStrength[16];

        int mvlimit = (slice->field_pic_flag || (p->mb_aff_frame_flag && MbQ->mb_field_decoding_flag)) ? 2 : 4;

        MbQ->DeblockCall = 1;
        short mb_x, mb_y;
        int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };
        get_mb_pos(p_Vid, MbQAddr, mb_size, &mb_x, &mb_y);

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

        if (p->mb_aff_frame_flag && MbQ->mb_field_decoding_flag&& mb_y == MB_BLOCK_SIZE)
            filterTopMbEdgeFlag = 0;

        if (MbQ->p_Slice->disable_deblocking_filter_idc == 2) {
            // don't filter at slice boundaries
            filterLeftMbEdgeFlag = MbQ->mbAvailA;
            filterTopMbEdgeFlag  = MbQ->mbAvailB;
            if (p->mb_aff_frame_flag && !MbQ->mb_field_decoding_flag && (MbQAddr & 0x01))
                filterTopMbEdgeFlag = 1;
        }

        if (p->mb_aff_frame_flag == 1) 
            CheckAvailabilityOfNeighborsMBAFF(MbQ);

        // Vertical deblocking
        for (int edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            if ((MbQ->CodedBlockPatternLuma == 0 && MbQ->CodedBlockPatternChroma == 0) &&
                (slice->slice_type == P_SLICE || slice->slice_type == B_SLICE)) {
                if (filterNon8x8LumaEdgesFlag[edge] == 0 && sps->chroma_format_idc != YUV444)
                    continue;
                else if (edge > 0) {
                    if (((MbQ->mb_type == PSKIP && slice->slice_type == P_SLICE) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P16x8)))
                        continue;
                    else if ((edge & 0x01) && ((MbQ->mb_type == P8x16) || (slice->slice_type == B_SLICE && MbQ->mb_type == BSKIP_DIRECT && sps->direct_8x8_inference_flag)))
                        continue;
                }
            }

            if (edge || filterLeftMbEdgeFlag) {
                if (p->mb_aff_frame_flag) {
                    // Strength for 4 blks in 1 stripe
                    Strength = pStrength;
                    get_strength_ver_MBAff(Strength, MbQ, edge * 4, mvlimit, p);
                }
                else {
                    Strength = MbQ->strength_ver[edge];
                    get_strength_ver(MbQ, edge, mvlimit, p);
                }

                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) { // only if one of the 16 Strength bytes is != 0
                //if ((*((int *) Strength))) { // only if one of the 16 Strength bytes is != 0
                //if (p_Strength64[0] || p_Strength64[1]) { // only if one of the 16 Strength bytes is != 0
                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        this->edge_loop_luma_ver(PLANE_Y, p->imgY, Strength, MbQ, edge * 4, p);
                        if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
                            this->edge_loop_luma_ver(PLANE_U, p->imgUV[0], Strength, MbQ, edge * 4, p);
                            this->edge_loop_luma_ver(PLANE_V, p->imgUV[1], Strength, MbQ, edge * 4, p);
                        }
                    }
                    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422) {
                        int edge_cr = chroma_edge[0][edge][p->chroma_format_idc == YUV422];
                        if (edge_cr >= 0) {
                            this->edge_loop_chroma_ver(PLANE_U, p->imgUV[0], Strength, MbQ, edge_cr, p);
                            this->edge_loop_chroma_ver(PLANE_V, p->imgUV[1], Strength, MbQ, edge_cr, p);
                        }
                    }
                }        
            }
        }

        // horizontal deblocking  
        for (int edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            if ((MbQ->CodedBlockPatternLuma == 0 && MbQ->CodedBlockPatternChroma == 0) && 
                (slice->slice_type == P_SLICE || slice->slice_type == B_SLICE)) {
                if (filterNon8x8LumaEdgesFlag[edge] == 0 && sps->chroma_format_idc == YUV420)
                    continue;
                else if (edge > 0) {
                    if (((MbQ->mb_type == PSKIP && slice->slice_type == P_SLICE) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P8x16)))
                        continue;
                    else if ((edge & 0x01) && ((MbQ->mb_type == P16x8) || (slice->slice_type == B_SLICE && MbQ->mb_type == BSKIP_DIRECT && sps->direct_8x8_inference_flag)))
                        continue;
                }
            }

            if (edge || filterTopMbEdgeFlag) {
                if (p->mb_aff_frame_flag) {
                    // Strength for 4 blks in 1 stripe
                    Strength = pStrength;
                    get_strength_hor_MBAff(Strength, MbQ, edge * 4, mvlimit, p);
                } else {
                    Strength = MbQ->strength_hor[edge];
                    get_strength_hor(MbQ, edge, mvlimit, p);
                }

                if ((*((int64_t *) Strength)) || ((*(((int64_t *) Strength) + 1)))) { // only if one of the 16 Strength bytes is != 0
                //if (p_Strength64[0] || p_Strength64[1]) { // only if one of the 16 Strength bytes is != 0
                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        this->edge_loop_luma_hor(PLANE_Y, p->imgY, Strength, MbQ, edge * 4, p);
                        if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
                            this->edge_loop_luma_hor(PLANE_U, p->imgUV[0], Strength, MbQ, edge * 4, p);
                            this->edge_loop_luma_hor(PLANE_V, p->imgUV[1], Strength, MbQ, edge * 4, p);
                        }
                    }
                    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422) {
                        int edge_cr = chroma_edge[1][edge][p->chroma_format_idc == YUV422];
                        if (edge_cr >= 0) {
                            this->edge_loop_chroma_hor(PLANE_U, p->imgUV[0], Strength, MbQ, edge_cr, p);
                            this->edge_loop_chroma_hor(PLANE_V, p->imgUV[1], Strength, MbQ, edge_cr, p);
                        }
                    }
                }

                if (!edge && !MbQ->mb_field_decoding_flag && MbQ->mixedModeEdgeFlag) {
                    // this is the extra horizontal edge between a frame macroblock pair and a field above it
                    MbQ->DeblockCall = 2;

                    if (p->mb_aff_frame_flag)
                        get_strength_hor_MBAff(Strength, MbQ, MB_BLOCK_SIZE, mvlimit, p); // Strength for 4 blks in 1 stripe
                    else
                        get_strength_hor(MbQ, 4, mvlimit, p); // Strength for 4 blks in 1 stripe

                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        this->edge_loop_luma_hor(PLANE_Y, p->imgY, Strength, MbQ, MB_BLOCK_SIZE, p);
                        if (sps->chroma_format_idc == YUV444 && !sps->separate_colour_plane_flag) {
                            this->edge_loop_luma_hor(PLANE_U, p->imgUV[0], Strength, MbQ, MB_BLOCK_SIZE, p);
                            this->edge_loop_luma_hor(PLANE_V, p->imgUV[1], Strength, MbQ, MB_BLOCK_SIZE, p);
                        }
                    }
                    if (sps->chroma_format_idc == YUV420 || sps->chroma_format_idc == YUV422) {
                        int edge_cr = chroma_edge[1][edge][p->chroma_format_idc == YUV422];
                        if (edge_cr >= 0) {
                            this->edge_loop_chroma_hor(PLANE_U, p->imgUV[0], Strength, MbQ, MB_BLOCK_SIZE, p);
                            this->edge_loop_chroma_hor(PLANE_V, p->imgUV[1], Strength, MbQ, MB_BLOCK_SIZE, p);
                        }
                    }

                    MbQ->DeblockCall = 1;
                }
            }
        }
    }

    MbQ->DeblockCall = 0;
}

void deblock_t::DeblockPicture(VideoParameters *p_Vid, storable_picture *p)
{
    this->init_neighbors(p_Vid);

    for (int i = 0; i < p->PicSizeInMbs; ++i)
        this->DeblockMb(p_Vid, p, i);
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
        int nsize = (p_Vid->dec_picture->size_y/BLOCK_SIZE)*(p_Vid->dec_picture->size_x/BLOCK_SIZE)*sizeof(pic_motion_params);
        memcpy( &(p_Vid->dec_picture->JVmv_info[PLANE_Y][0][0]), &(p_Vid->dec_picture_JV[PLANE_Y]->mv_info[0][0]), nsize);
        memcpy( &(p_Vid->dec_picture->JVmv_info[PLANE_U][0][0]), &(p_Vid->dec_picture_JV[PLANE_U]->mv_info[0][0]), nsize);
        memcpy( &(p_Vid->dec_picture->JVmv_info[PLANE_V][0][0]), &(p_Vid->dec_picture_JV[PLANE_V]->mv_info[0][0]), nsize);
    }

    // This could be done with pointers and seems not necessary
    for (int uv = 0; uv < 2; uv++) {
        for (int line = 0; line < sps->FrameHeightInMbs * 16; line++) {
            int nsize = sizeof(imgpel) * sps->PicWidthInMbs * 16;
            memcpy( p_Vid->dec_picture->imgUV[uv][line], p_Vid->dec_picture_JV[uv+1]->imgY[line], nsize );
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
                DeblockPicture(p_Vid, p);
            }
            p_Vid->ppSliceList[0]->colour_plane_id = colour_plane_id;
        } else {
            DeblockPicture(p_Vid, p);
        }
    }

    if (p_Vid->active_sps->separate_colour_plane_flag)
        this->make_frame_picture_JV(p_Vid);
}
