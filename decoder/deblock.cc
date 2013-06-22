
/*!
 *************************************************************************************
 * \file loopFilter.c
 *
 * \brief
 *    Filter to reduce blocking artifacts on a macroblock level.
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
#include "deblock.h"
#include "deblock_common.h"


static const char chroma_edge[2][4][4] = //[dir][edge][yuv_format]
{ { {-4, 0, 0, 0},
    {-4,-4,-4, 4},
    {-4, 4, 4, 8},
    {-4,-4,-4, 12}},

  { {-4, 0,  0,  0},
    {-4,-4,  4,  4},
    {-4, 4,  8,  8},
    {-4,-4, 12, 12}}};



// likely already set - see testing via asserts
static void init_neighbors(VideoParameters *p_Vid)
{
    int i, j;
    int width = p_Vid->PicWidthInMbs;
    int height = p_Vid->PicHeightInMbs;
    int size = p_Vid->PicSizeInMbs;
    Macroblock *currMB = &p_Vid->mb_data[0];
    // do the top left corner
    currMB->mbup = NULL;
    currMB->mbleft = NULL;
    currMB++;
    // do top row
    for (i = 1; i < width; i++) {
        currMB->mbup = NULL;
        currMB->mbleft = currMB - 1;
        currMB++;
    }

    // do left edge
    for (i = width; i < size; i += width) {
        currMB->mbup = currMB - width;
        currMB->mbleft = NULL;   
        currMB += width;
    }
    // do all others
    for (j = width + 1; j < width * height + 1; j += width) {
        currMB = &p_Vid->mb_data[j];
        for (i = 1; i < width; i++) {
            currMB->mbup   = currMB - width;
            currMB->mbleft = currMB - 1;
            currMB++;
        }
    }
}


/*!
 *****************************************************************************************
 * \brief
 *    Deblocking filter for one macroblock.
 *****************************************************************************************
 */

static void DeblockMb(VideoParameters *p_Vid, StorablePicture *p, int MbQAddr)
{
    Macroblock *MbQ = &(p_Vid->mb_data[MbQAddr]) ; // current Mb

    // return, if filter is disabled
    if (MbQ->DFDisableIdc == 1) 
        MbQ->DeblockCall = (Boolean)0;
    else {
        int           edge;

        byte         *Strength;
        byte          pStrength[16];
        short         mb_x, mb_y;

        int           filterNon8x8LumaEdgesFlag[4] = {1,1,1,1};
        int           filterLeftMbEdgeFlag;
        int           filterTopMbEdgeFlag;
        int           edge_cr;

        imgpel     **imgY = p->imgY;
        imgpel   ***imgUV = p->imgUV;
        Slice  *currSlice = MbQ->p_Slice;
        int       mvlimit = ((p->structure!=FRAME) || (p->mb_aff_frame_flag && MbQ->mb_field)) ? 2 : 4;

        seq_parameter_set_rbsp_t *active_sps = p_Vid->active_sps;

        MbQ->DeblockCall = (Boolean)1;
        get_mb_pos (p_Vid, MbQAddr, p_Vid->mb_size[IS_LUMA], &mb_x, &mb_y);

        if (MbQ->mb_type == I8MB)
            assert(MbQ->luma_transform_size_8x8_flag);

//      fieldMbInFrameFlag = MbaffFrameFlag == 1 && mb_field_decoding_flag == 1

//      filterInternalEdgesFlag = slice->disable_deblocking_filter_idc == 0

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
        filterNon8x8LumaEdgesFlag[3] = !(MbQ->luma_transform_size_8x8_flag);

        filterLeftMbEdgeFlag = (mb_x != 0);
        filterTopMbEdgeFlag  = (mb_y != 0);

        if (p->mb_aff_frame_flag && mb_y == MB_BLOCK_SIZE && MbQ->mb_field)
            filterTopMbEdgeFlag = 0;

        if (MbQ->DFDisableIdc == 2) {
            // don't filter at slice boundaries
            filterLeftMbEdgeFlag = MbQ->mbAvailA;
            // if this the bottom of a frame macroblock pair then always filter the top edge
            filterTopMbEdgeFlag  = (p->mb_aff_frame_flag && !MbQ->mb_field && (MbQAddr & 0x01)) ? 1 : MbQ->mbAvailB;
        }

        if (p->mb_aff_frame_flag == 1) 
            CheckAvailabilityOfNeighborsMBAFF(MbQ);

        // Vertical deblocking
        for (edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            if (MbQ->cbp == 0 && (currSlice->slice_type == P_SLICE || currSlice->slice_type == B_SLICE)) {
                if (filterNon8x8LumaEdgesFlag[edge] == 0 && active_sps->chroma_format_idc != YUV444)
                    continue;
                else if (edge > 0) {
                    if (((MbQ->mb_type == PSKIP && currSlice->slice_type == P_SLICE) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P16x8)))
                        continue;
                    else if ((edge & 0x01) && ((MbQ->mb_type == P8x16) || (currSlice->slice_type == B_SLICE && MbQ->mb_type == BSKIP_DIRECT && active_sps->direct_8x8_inference_flag)))
                        continue;
                }
            }

            if (edge || filterLeftMbEdgeFlag) {
                if (p->mb_aff_frame_flag) {
                    // Strength for 4 blks in 1 stripe
                    Strength = pStrength;
                    get_strength_ver_MBAff(Strength, MbQ, edge << 2, mvlimit, p);
                }
                else {
                    Strength = MbQ->strength_ver[edge];
                    p_Vid->GetStrengthVer(MbQ, edge, mvlimit, p);
                }

                if ((*((int64 *) Strength)) || ((*(((int64 *) Strength) + 1)))) { // only if one of the 16 Strength bytes is != 0
                //if ((*((int *) Strength))) { // only if one of the 16 Strength bytes is != 0
                //if (p_Strength64[0] || p_Strength64[1]) { // only if one of the 16 Strength bytes is != 0
                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        p_Vid->EdgeLoopLumaVer( PLANE_Y, imgY, Strength, MbQ, edge << 2, p);
                        if (currSlice->chroma444_not_separate) {
                            p_Vid->EdgeLoopLumaVer(PLANE_U, imgUV[0], Strength, MbQ, edge << 2, p);
                            p_Vid->EdgeLoopLumaVer(PLANE_V, imgUV[1], Strength, MbQ, edge << 2, p);
                        }
                    }
                    if (active_sps->chroma_format_idc==YUV420 || active_sps->chroma_format_idc==YUV422) {
                        edge_cr = chroma_edge[0][edge][p->chroma_format_idc];
                        if (imgUV != NULL && edge_cr >= 0) {
                            p_Vid->EdgeLoopChromaVer( imgUV[0], Strength, MbQ, edge_cr, 0, p);
                            p_Vid->EdgeLoopChromaVer( imgUV[1], Strength, MbQ, edge_cr, 1, p);
                        }
                    }
                }        
            }
        }//end edge

        // horizontal deblocking  
        for (edge = 0; edge < 4; ++edge) {
            // If cbp == 0 then deblocking for some macroblock types could be skipped
            if (MbQ->cbp == 0 && (currSlice->slice_type == P_SLICE || currSlice->slice_type == B_SLICE)) {
                if (filterNon8x8LumaEdgesFlag[edge] == 0 && active_sps->chroma_format_idc==YUV420)
                    continue;
                else if (edge > 0) {
                    if (((MbQ->mb_type == PSKIP && currSlice->slice_type == P_SLICE) || (MbQ->mb_type == P16x16) || (MbQ->mb_type == P8x16)))
                        continue;
                    else if ((edge & 0x01) && ((MbQ->mb_type == P16x8) || (currSlice->slice_type == B_SLICE && MbQ->mb_type == BSKIP_DIRECT && active_sps->direct_8x8_inference_flag)))
                        continue;
                }
            }

            if (edge || filterTopMbEdgeFlag) {
                if (p->mb_aff_frame_flag) {
                    // Strength for 4 blks in 1 stripe
                    Strength = pStrength;
                    get_strength_hor_MBAff(Strength, MbQ, edge << 2, mvlimit, p);
                } else {
                    Strength = MbQ->strength_hor[edge];
                    p_Vid->GetStrengthHor(MbQ, edge, mvlimit, p);
                }

                if ((*((int64 *) Strength)) || ((*(((int64 *) Strength) + 1)))) { // only if one of the 16 Strength bytes is != 0
                //if (p_Strength64[0] || p_Strength64[1]) { // only if one of the 16 Strength bytes is != 0
                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        p_Vid->EdgeLoopLumaHor( PLANE_Y, imgY, Strength, MbQ, edge << 2, p) ;
                        if (currSlice->chroma444_not_separate) {
                            p_Vid->EdgeLoopLumaHor(PLANE_U, imgUV[0], Strength, MbQ, edge << 2, p);
                            p_Vid->EdgeLoopLumaHor(PLANE_V, imgUV[1], Strength, MbQ, edge << 2, p);
                        }
                    }
                    if (active_sps->chroma_format_idc==YUV420 || active_sps->chroma_format_idc==YUV422) {
                        edge_cr = chroma_edge[1][edge][p->chroma_format_idc];
                        if (imgUV != NULL && edge_cr >= 0) {
                            p_Vid->EdgeLoopChromaHor( imgUV[0], Strength, MbQ, edge_cr, 0, p);
                            p_Vid->EdgeLoopChromaHor( imgUV[1], Strength, MbQ, edge_cr, 1, p);
                        }
                    }
                }

                if (!edge && !MbQ->mb_field && MbQ->mixedModeEdgeFlag) {
                    // this is the extra horizontal edge between a frame macroblock pair and a field above it
                    MbQ->DeblockCall = (Boolean)2;
                    if (p->mb_aff_frame_flag)
                        get_strength_hor_MBAff(Strength, MbQ, MB_BLOCK_SIZE, mvlimit, p); // Strength for 4 blks in 1 stripe
                    else
                        p_Vid->GetStrengthHor(MbQ, 4, mvlimit, p); // Strength for 4 blks in 1 stripe

                    if (filterNon8x8LumaEdgesFlag[edge]) {
                        p_Vid->EdgeLoopLumaHor(PLANE_Y, imgY, Strength, MbQ, MB_BLOCK_SIZE, p) ;
                        if (currSlice->chroma444_not_separate) {
                            p_Vid->EdgeLoopLumaHor(PLANE_U, imgUV[0], Strength, MbQ, MB_BLOCK_SIZE, p) ;
                            p_Vid->EdgeLoopLumaHor(PLANE_V, imgUV[1], Strength, MbQ, MB_BLOCK_SIZE, p) ;
                        }
                    }
                    if (active_sps->chroma_format_idc==YUV420 || active_sps->chroma_format_idc==YUV422) {
                        edge_cr = chroma_edge[1][edge][p->chroma_format_idc];
                        if (imgUV != NULL && edge_cr >= 0) {
                            p_Vid->EdgeLoopChromaHor( imgUV[0], Strength, MbQ, MB_BLOCK_SIZE, 0, p) ;
                            p_Vid->EdgeLoopChromaHor( imgUV[1], Strength, MbQ, MB_BLOCK_SIZE, 1, p) ;
                        }
                    }

                    MbQ->DeblockCall = (Boolean)1;
                }
            }
        }//end edge  

        MbQ->DeblockCall = (Boolean)0;
    }
}

static void DeblockPicture(VideoParameters *p_Vid, StorablePicture *p)
{
    int i;
    for (i = 0; i < p->PicSizeInMbs; ++i)
        DeblockMb(p_Vid, p, i);
}

static void make_frame_picture_JV(VideoParameters *p_Vid)
{
    int uv, line;
    int nsize;
    p_Vid->dec_picture = p_Vid->dec_picture_JV[0];

    if (p_Vid->dec_picture->used_for_reference) {
        nsize = (p_Vid->dec_picture->size_y/BLOCK_SIZE)*(p_Vid->dec_picture->size_x/BLOCK_SIZE)*sizeof(PicMotionParams);
        memcpy( &(p_Vid->dec_picture->JVmv_info[PLANE_Y][0][0]), &(p_Vid->dec_picture_JV[PLANE_Y]->mv_info[0][0]), nsize);
        memcpy( &(p_Vid->dec_picture->JVmv_info[PLANE_U][0][0]), &(p_Vid->dec_picture_JV[PLANE_U]->mv_info[0][0]), nsize);
        memcpy( &(p_Vid->dec_picture->JVmv_info[PLANE_V][0][0]), &(p_Vid->dec_picture_JV[PLANE_V]->mv_info[0][0]), nsize);
    }

    // This could be done with pointers and seems not necessary
    for (uv = 0; uv < 2; uv++) {
        for (line = 0; line < p_Vid->height; line++) {
            nsize = sizeof(imgpel) * p_Vid->width;
            memcpy( p_Vid->dec_picture->imgUV[uv][line], p_Vid->dec_picture_JV[uv+1]->imgY[line], nsize );
        }
        free_storable_picture(p_Vid->dec_picture_JV[uv+1]);
    }
}

void init_Deblock(VideoParameters *p_Vid, int mb_aff_frame_flag)
{
    if (p_Vid->yuv_format == YUV444 && p_Vid->separate_colour_plane_flag) {
        change_plane_JV(p_Vid, PLANE_Y, NULL);
        init_neighbors(p_Dec->p_Vid);
        change_plane_JV(p_Vid, PLANE_U, NULL);
        init_neighbors(p_Dec->p_Vid);
        change_plane_JV(p_Vid, PLANE_V, NULL);
        init_neighbors(p_Dec->p_Vid);
        change_plane_JV(p_Vid, PLANE_Y, NULL);
    } else 
        init_neighbors(p_Dec->p_Vid);
    if (mb_aff_frame_flag == 1) 
        set_loop_filter_functions_mbaff(p_Vid);
    else
        set_loop_filter_functions_normal(p_Vid);
}

void change_plane_JV(VideoParameters *p_Vid, int nplane, Slice *pSlice)
{
    p_Vid->mb_data     = p_Vid->mb_data_JV    [nplane];
    p_Vid->dec_picture = p_Vid->dec_picture_JV[nplane];
    p_Vid->siblock     = p_Vid->siblock_JV    [nplane];
    p_Vid->ipredmode   = p_Vid->ipredmode_JV  [nplane];
    p_Vid->intra_block = p_Vid->intra_block_JV[nplane];

    if (pSlice) {
        pSlice->mb_data     = p_Vid->mb_data_JV    [nplane];
        pSlice->dec_picture = p_Vid->dec_picture_JV[nplane];
        pSlice->siblock     = p_Vid->siblock_JV    [nplane];
        pSlice->ipredmode   = p_Vid->ipredmode_JV  [nplane];
        pSlice->intra_block = p_Vid->intra_block_JV[nplane];
    }
}

void pic_deblock(VideoParameters *p_Vid, StorablePicture *p)
{
    if (!p_Vid->iDeblockMode && (p_Vid->bDeblockEnable & (1<<p->used_for_reference))) {
        //deblocking for frame or field
        if (p_Vid->separate_colour_plane_flag != 0) {
            int nplane;
            int colour_plane_id = p_Vid->ppSliceList[0]->colour_plane_id;
            for (nplane = 0; nplane<MAX_PLANE; ++nplane) {
                p_Vid->ppSliceList[0]->colour_plane_id = nplane;
                change_plane_JV(p_Vid, nplane, NULL);
                DeblockPicture(p_Vid, p);
            }
            p_Vid->ppSliceList[0]->colour_plane_id = colour_plane_id;
            make_frame_picture_JV(p_Vid);
        } else
            DeblockPicture(p_Vid, p);
    } else {
        if (p_Vid->separate_colour_plane_flag != 0)
            make_frame_picture_JV(p_Vid);
    }
}
