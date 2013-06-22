
/*!
 *************************************************************************************
 * \file loop_filter_mbaff.c
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
#include "loop_filter.h"
#include "loop_filter_common.h"

static void edge_loop_luma_ver_MBAff   (ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p);
static void edge_loop_luma_hor_MBAff   (ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p);
static void edge_loop_chroma_ver_MBAff (imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p);
static void edge_loop_chroma_hor_MBAff (imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p);

void set_loop_filter_functions_mbaff(VideoParameters *p_Vid)
{
  p_Vid->EdgeLoopLumaVer   = edge_loop_luma_ver_MBAff;
  p_Vid->EdgeLoopLumaHor   = edge_loop_luma_hor_MBAff;
  p_Vid->EdgeLoopChromaVer = edge_loop_chroma_ver_MBAff;
  p_Vid->EdgeLoopChromaHor = edge_loop_chroma_hor_MBAff;
}


Macroblock* get_non_aff_neighbor_luma(Macroblock *mb, int xN, int yN)
{
  if (xN < 0)
    return(mb->mbleft);
  else if (yN < 0)
    return(mb->mbup);
  else
    return(mb);
}

Macroblock* get_non_aff_neighbor_chroma(Macroblock *mb, int xN, int yN, int block_width,int block_height)
{
  if (xN < 0) 
  {
    if (yN < block_height)
      return(mb->mbleft);
    else
      return(NULL);
  }
  else if (xN < block_width) 
  {
    if (yN < 0)
      return(mb->mbup);
    else if (yN < block_height)
      return(mb);
    else
      return(NULL);
  }
  else
    return(NULL);
}


/*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for MBAFF)
 *********************************************************************************************
 */
void get_strength_ver_MBAff(byte *Strength, Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p)
{
  //byte *Strength = MbQ->strength_ver[edge];
    short  blkP, blkQ, idx;

    int    StrValue;
    short  mb_x, mb_y;

    Macroblock *MbP;

    PixelPos pixP;
    VideoParameters *p_Vid = MbQ->p_Vid;
    BlockPos *PicPos = p_Vid->PicPos;

    if (p->slice_type == SP_SLICE || p->slice_type == SI_SLICE) {
        for (idx = 0; idx < MB_BLOCK_SIZE; ++idx) {
            getAffNeighbour(MbQ, edge - 1, idx, p_Vid->mb_size[IS_LUMA], &pixP);
            blkQ = (short) ((idx & 0xFFFC) + (edge >> 2));
            blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

            MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            MbQ->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field);

            Strength[idx] = (edge == 0) ? 4 : 3;
        }
    } else {
        getAffNeighbour(MbQ, edge - 1, 0, p_Vid->mb_size[IS_LUMA], &pixP);

        MbP = &(p_Vid->mb_data[pixP.mb_addr]);
        // Neighboring Frame MBs
        if (MbQ->mb_field == FALSE && MbP->mb_field == FALSE) {
            MbQ->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field); 
            if (MbQ->is_intra_block == TRUE || MbP->is_intra_block == TRUE) {
                //printf("idx %d %d %d %d %d\n", idx, pixP.x, pixP.y, pixP.pos_x, pixP.pos_y);
                // Start with Strength=3. or Strength=4 for Mb-edge
                StrValue = (edge == 0) ? 4 : 3;
                memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
            } else {
                get_mb_block_pos_mbaff (PicPos, MbQ->mbAddrX, &mb_x, &mb_y);
                for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
                    blkQ = (short) ((idx & 0xFFFC) + (edge >> 2));
                    blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

                    if (((MbQ->s_cbp[0].blk & i64_power2(blkQ)) != 0) || ((MbP->s_cbp[0].blk & i64_power2(blkP)) != 0))
                        StrValue = 2;
                    else if (edge && ((MbQ->mb_type == 1)  || (MbQ->mb_type == 2)))
                        StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
                    else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
                        int blk_y  = ((mb_y<<2) + (blkQ >> 2));
                        int blk_x  = ((mb_x<<2) + (blkQ  & 3));
                        int blk_y2 = (pixP.pos_y >> 2);
                        int blk_x2 = (pixP.pos_x >> 2);
                        PicMotionParams *mv_info_p = &p->mv_info[blk_y ][blk_x ];
                        PicMotionParams *mv_info_q = &p->mv_info[blk_y2][blk_x2];
                        StorablePicturePtr ref_p0 = mv_info_p->ref_pic[LIST_0];
                        StorablePicturePtr ref_q0 = mv_info_q->ref_pic[LIST_0];
                        StorablePicturePtr ref_p1 = mv_info_p->ref_pic[LIST_1];
                        StorablePicturePtr ref_q1 = mv_info_q->ref_pic[LIST_1];

                        if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1))||((ref_p0==ref_q1) && (ref_p1==ref_q0))) {
                            // L0 and L1 reference pictures of p0 are different; q0 as well
                            if (ref_p0 != ref_p1) {
                                // compare MV for the same reference picture
                                if (ref_p0 == ref_q0) {
                                    StrValue =  (byte) (
                                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit));
                                } else {
                                    StrValue =  (byte) (
                                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit));
                                }
                            } else { // L0 and L1 reference pictures of p0 are the same; q0 as well
                                StrValue = (byte) ((
                                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                                    &&(
                                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)));
                            }
                        } else
                            StrValue = 1;
                    }
                    *(int*)(Strength+idx) = StrValue * 0x01010101;
                    pixP.y += 4;
                    pixP.pos_y += 4;
                }
            }
        } else {
            for (idx = 0; idx < MB_BLOCK_SIZE; ++idx) {
                getAffNeighbour(MbQ, edge - 1, idx, p_Vid->mb_size[IS_LUMA], &pixP);
                blkQ = (short) ((idx & 0xFFFC) + (edge >> 2));
                blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

                MbP = &(p_Vid->mb_data[pixP.mb_addr]);
                MbQ->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field); 

                // Start with Strength=3. or Strength=4 for Mb-edge
                Strength[idx] = (edge == 0 && (((!p->mb_aff_frame_flag && (p->structure==FRAME)) ||
                    (p->mb_aff_frame_flag && !MbP->mb_field && !MbQ->mb_field)) ||
                    ((p->mb_aff_frame_flag || (p->structure!=FRAME))))) ? 4 : 3;

                if (MbQ->is_intra_block == FALSE && MbP->is_intra_block == FALSE) {
                    if (((MbQ->s_cbp[0].blk & i64_power2(blkQ)) != 0) || ((MbP->s_cbp[0].blk & i64_power2(blkP)) != 0))
                        Strength[idx] = 2 ;
                    else {
                        // if no coefs, but vector difference >= 1 set Strength=1
                        // if this is a mixed mode edge then one set of reference pictures will be frame and the
                        // other will be field
                        if (MbQ->mixedModeEdgeFlag)
                            Strength[idx] = 1;
                        else {
                            get_mb_block_pos_mbaff (PicPos, MbQ->mbAddrX, &mb_x, &mb_y);

                            int blk_y  = ((mb_y<<2) + (blkQ >> 2));
                            int blk_x  = ((mb_x<<2) + (blkQ  & 3));
                            int blk_y2 = (pixP.pos_y >> 2);
                            int blk_x2 = (pixP.pos_x >> 2);
                            PicMotionParams *mv_info_p = &p->mv_info[blk_y ][blk_x ];
                            PicMotionParams *mv_info_q = &p->mv_info[blk_y2][blk_x2];
                            StorablePicturePtr ref_p0 = mv_info_p->ref_pic[LIST_0];
                            StorablePicturePtr ref_q0 = mv_info_q->ref_pic[LIST_0];
                            StorablePicturePtr ref_p1 = mv_info_p->ref_pic[LIST_1];
                            StorablePicturePtr ref_q1 = mv_info_q->ref_pic[LIST_1];

                            if ( ((ref_p0==ref_q0) && (ref_p1==ref_q1))||((ref_p0==ref_q1) && (ref_p1==ref_q0))) {
                                Strength[idx]=0;
                                // L0 and L1 reference pictures of p0 are different; q0 as well
                                if (ref_p0 != ref_p1) {
                                    // compare MV for the same reference picture
                                    if (ref_p0==ref_q0) {
                                        Strength[idx] =  (byte) (
                                            compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                                            compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit));
                                    } else {
                                        Strength[idx] =  (byte) (
                                            compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                                            compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit));
                                    }
                                } else { // L0 and L1 reference pictures of p0 are the same; q0 as well
                                    Strength[idx] = (byte) ((
                                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                                        &&(
                                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)));
                                }
                            } else
                                Strength[idx] = 1;
                        }
                    }
                }
            }
        }
    }
}

/*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for MBAFF)
 *********************************************************************************************
 */
void get_strength_hor_MBAff(byte *Strength, Macroblock *MbQ, int edge, int mvlimit, StorablePicture *p)
{
    short  blkP, blkQ, idx;
    short  blk_x, blk_x2, blk_y, blk_y2 ;

    int    StrValue;
    int    xQ, yQ = (edge < MB_BLOCK_SIZE ? edge : 1);
    short  mb_x, mb_y;

    Macroblock *MbP;

    PixelPos pixP;
    VideoParameters *p_Vid = MbQ->p_Vid;
    BlockPos *PicPos = p_Vid->PicPos;

    if (p->slice_type == SP_SLICE || p->slice_type == SI_SLICE) {
        for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
            xQ = idx;
            getAffNeighbour(MbQ, xQ, yQ - 1, p_Vid->mb_size[IS_LUMA], &pixP);

            blkQ = (short) ((yQ & 0xFFFC) + (xQ >> 2));
            blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

            MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            MbQ->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field);

            StrValue = (edge == 0 && (!MbP->mb_field && !MbQ->mb_field)) ? 4 : 3;
      
            *(int*)(Strength+idx) = StrValue * 0x01010101;
        }
    } else {
        getAffNeighbour(MbQ, 0, yQ - 1, p_Vid->mb_size[IS_LUMA], &pixP);
        MbP = &(p_Vid->mb_data[pixP.mb_addr]);
        MbQ->mixedModeEdgeFlag = (byte) (MbQ->mb_field != MbP->mb_field); 

        // Set intra mode deblocking
        if (MbQ->is_intra_block == TRUE || MbP->is_intra_block == TRUE) {
            StrValue = (edge == 0 && (!MbP->mb_field && !MbQ->mb_field)) ? 4 : 3;
            memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
        } else {
            for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
                xQ = idx;    
                getAffNeighbour(MbQ, xQ, yQ - 1, p_Vid->mb_size[IS_LUMA], &pixP);

                blkQ = (short) ((yQ & 0xFFFC) + (xQ >> 2));
                blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

                if (((MbQ->s_cbp[0].blk & i64_power2(blkQ)) != 0) || ((MbP->s_cbp[0].blk & i64_power2(blkP)) != 0))
                    StrValue = 2;
                else {
                    // if no coefs, but vector difference >= 1 set Strength=1
                    // if this is a mixed mode edge then one set of reference pictures will be frame and the
                    // other will be field
                    if (MbQ->mixedModeEdgeFlag)
                        StrValue = 1;
                    else {
                        get_mb_block_pos_mbaff (PicPos, MbQ->mbAddrX, &mb_x, &mb_y);

                        blk_y  = (short) ((mb_y<<2) + (blkQ >> 2));
                        blk_x  = (short) ((mb_x<<2) + (blkQ  & 3));
                        blk_y2 = (short) (pixP.pos_y >> 2);
                        blk_x2 = (short) (pixP.pos_x >> 2);
                        PicMotionParams *mv_info_p = &p->mv_info[blk_y ][blk_x ];
                        PicMotionParams *mv_info_q = &p->mv_info[blk_y2][blk_x2];
                        StorablePicturePtr ref_p0 = mv_info_p->ref_pic[LIST_0];
                        StorablePicturePtr ref_q0 = mv_info_q->ref_pic[LIST_0];
                        StorablePicturePtr ref_p1 = mv_info_p->ref_pic[LIST_1];
                        StorablePicturePtr ref_q1 = mv_info_q->ref_pic[LIST_1];

                        if ((ref_p0 == ref_q0 && ref_p1 == ref_q1) || (ref_p0 == ref_q1 && ref_p1 == ref_q0)) {
                            StrValue = 0;
                            // L0 and L1 reference pictures of p0 are different; q0 as well
                            if (ref_p0 != ref_p1) {
                                // compare MV for the same reference picture
                                if (ref_p0==ref_q0) {
                                    StrValue =  (byte) (
                                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit));
                                } else {
                                    StrValue =  (byte) (
                                        compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                                        compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit));
                                }
                            } else { // L0 and L1 reference pictures of p0 are the same; q0 as well
                                StrValue = (byte) ((
                                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_0], mvlimit) ||
                                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_1], mvlimit))
                                    &&(
                                    compare_mvs(&mv_info_p->mv[LIST_0], &mv_info_q->mv[LIST_1], mvlimit) ||
                                    compare_mvs(&mv_info_p->mv[LIST_1], &mv_info_q->mv[LIST_0], mvlimit)));
                            }
                        } else
                            StrValue = 1;
                    }
                }
                *(int*)(Strength+idx) = StrValue * 0x01010101;
            }
        }
    }
}


static void deblock_mb_mbaff(bool verticalEdgeFlag, bool chromaEdgeFlag, bool chromaStyleFilteringFlag,
                       Macroblock *MbP, Macroblock *MbQ, byte *Strength,
                       ColorPlane pl, imgpel **Img, int edge, StorablePicture *p)
{
    VideoParameters *p_Vid = MbQ->p_Vid;

    int width = chromaEdgeFlag == 0 ? p->iLumaStride : p->iChromaStride;
    int PelNum = pl ? pelnum_cr[verticalEdgeFlag ? 0 : 1][p->chroma_format_idc] : MB_BLOCK_SIZE;
    int bitdepth_scale = pl ? p_Vid->bitdepth_scale[IS_CHROMA]
                            : p_Vid->bitdepth_scale[IS_LUMA];
    int *mb_size = p_Vid->mb_size[chromaEdgeFlag == 0 ? IS_LUMA : IS_CHROMA];
    PixelPos pixP, pixQ;

    int yQ = (edge < MB_BLOCK_SIZE ? edge : 1);

    if (MbQ->DFDisableIdc == 0) {
        for (int pel = 0; pel < PelNum; ++pel) {
            if (verticalEdgeFlag) {
                getAffNeighbour(MbQ, edge - 1, pel, mb_size, &pixP);     
                getAffNeighbour(MbQ, edge    , pel, mb_size, &pixQ);
            } else {
                getAffNeighbour(MbQ, pel, yQ - 1, mb_size, &pixP);     
                getAffNeighbour(MbQ, pel, yQ    , mb_size, &pixQ);
            }

            Macroblock *MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            int StrengthIdx = (PelNum == 8) ? ((MbQ->mb_field && !MbP->mb_field) ? pel << 1 : ((pel >> 1) << 2) + (pel & 0x01)) : pel;
            int bS = Strength[StrengthIdx];

            if (pixP.available && bS != 0) {
                imgpel *SrcPtrP = &(Img[pixP.pos_y][pixP.pos_x]);
                imgpel *SrcPtrQ = &(Img[pixQ.pos_y][pixQ.pos_x]);
                int incP, incQ;

                if (verticalEdgeFlag) {
                    incQ = 1;
                    incP = 1;
                } else {
                    incQ = ((MbP->mb_field && !MbQ->mb_field) ? 2 * width : width);
                    incP = ((MbQ->mb_field && !MbP->mb_field) ? 2 * width : width);
                }

                // Average QP of the two blocks
                int QP = pl? ((MbP->qpc[pl-1] + MbQ->qpc[pl-1] + 1) >> 1)
                            : (MbP->qp + MbQ->qp + 1) >> 1;
                int indexA = iClip3(0, MAX_QP, QP + MbQ->DFAlphaC0Offset);
                int indexB = iClip3(0, MAX_QP, QP + MbQ->DFBetaOffset);
                int Alpha  = ALPHA_TABLE[indexA] * bitdepth_scale;
                int Beta   = BETA_TABLE [indexB] * bitdepth_scale;

                if (bS == 4)
                    deblock_strong(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag);
                else if (bS > 0)
                    deblock_normal(SrcPtrP, SrcPtrQ, incP, incQ, Alpha, Beta, bS, chromaEdgeFlag, pl, p_Vid->bitdepth_luma, p_Vid->bitdepth_chroma, indexA);
            }
        }
    }
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Super MB Frame coded MBs
 *****************************************************************************************
 */
static void edge_loop_luma_ver_MBAff(ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p)
{
    deblock_mb_mbaff(1, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}

/*!
 *****************************************************************************************
 * \brief
 *    Filters 16 pel block edge of Super MB Frame coded MBs
 *****************************************************************************************
 */
static void edge_loop_luma_hor_MBAff(ColorPlane pl, imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, StorablePicture *p)
{
    deblock_mb_mbaff(0, 0, 0, MbQ, MbQ, Strength, pl, Img, edge, p);
}



/*!
*****************************************************************************************
* \brief
*    Filters chroma block edge for MBAFF types
*****************************************************************************************
 */
static void edge_loop_chroma_ver_MBAff(imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p)
{
    deblock_mb_mbaff(1, 1, 1, MbQ, MbQ, Strength, (ColorPlane)(uv+1), Img, edge, p);
}

/*!
*****************************************************************************************
* \brief
*    Filters chroma block edge for MBAFF types
*****************************************************************************************
 */
static void edge_loop_chroma_hor_MBAff(imgpel** Img, byte *Strength, Macroblock *MbQ, int edge, int uv, StorablePicture *p)
{
    deblock_mb_mbaff(0, 1, 1, MbQ, MbQ, Strength, (ColorPlane)(uv+1), Img, edge, p);
}
