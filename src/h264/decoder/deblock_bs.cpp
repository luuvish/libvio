#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "image.h"
#include "neighbour.h"
#include "deblock.h"


#define get_x_luma(x) (x & 15)
#define get_y_luma(y) (y & 15)
#define get_pos_x_luma(m,xx) (m->mb.x * 16 + (xx & 15))
#define get_pos_y_luma(m,yy) (m->mb.y * 16 + (yy & 15))


static inline int compare_mvs(const MotionVector *mv0, const MotionVector *mv1, int mvlimit)
{
    return (abs( mv0->mv_x - mv1->mv_x) >= 4) | (abs( mv0->mv_y - mv1->mv_y) >= mvlimit);
}


static mb_t* get_non_aff_neighbor_luma(mb_t *mb, int xN, int yN)
{
    if (xN < 0)
        return mb->mbleft;
    else if (yN < 0)
        return mb->mbup;
    else
        return mb;
}


  /*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for different Frame or Field types)
 *********************************************************************************************
 */
void get_strength_ver(mb_t *MbQ, int edge, int mvlimit, storable_picture *p)
{
    byte*     Strength = MbQ->strength_ver[edge];
    slice_t*  currSlice = MbQ->p_Slice;
    int       StrValue;
    BlockPos* PicPos = MbQ->p_Vid->PicPos;

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
    if (currSlice->slice_type == SP_SLICE || currSlice->slice_type == SI_SLICE ||
        MbQ->is_intra_block || MbP->is_intra_block) {
        // Set strength to either 3 or 4 regardless of pixel position
        StrValue = (edge == 0) ? 4 : 3;
        memset(Strength, (byte) StrValue, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && currSlice->slice_type == P_SLICE && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, BLOCK_SIZE * sizeof(byte));
        return;
    }
/*
    if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P16x8)) {
        for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
            int blkQ = idx + (edge);
            int blkP = idx + (get_x_luma(xQ) >> 2);
            if (MbQ->cbp_blks[0] & (((uint64_t)1 << blkQ) | ((uint64_t)1 << blkP)))
                StrValue = 2;
            else
                StrValue = 0; // if internal edge of certain types, then we already know StrValue should be 0
            Strength[idx >> 2] = StrValue;
        }
        return;
    }
*/
    BlockPos mb = PicPos[MbQ->mbAddrX];
    mb.x <<= BLOCK_SHIFT;
    mb.y <<= BLOCK_SHIFT;

    for (int idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
        int blkQ = idx + (edge);
        int blkP = idx + (get_x_luma(xQ) >> 2);
        if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 ||
            (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
            int blk_x  = mb.x + (blkQ  & 3);
            int blk_y  = mb.y + (blkQ >> 2);
            int blk_x2 = (get_pos_x_luma(neighbor, xQ)      ) >> 2;
            int blk_y2 = (get_pos_y_luma(neighbor,  0) + idx) >> 2;
            pic_motion_params *mv_info_p = &p->mv_info[blk_y ][blk_x ];            
            pic_motion_params *mv_info_q = &p->mv_info[blk_y2][blk_x2];            
            storable_picture* ref_p0 = mv_info_p->ref_pic[LIST_0];
            storable_picture* ref_q0 = mv_info_q->ref_pic[LIST_0];            
            storable_picture* ref_p1 = mv_info_p->ref_pic[LIST_1];
            storable_picture* ref_q1 = mv_info_q->ref_pic[LIST_1];

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
void get_strength_hor(mb_t *MbQ, int edge, int mvlimit, storable_picture *p)
{  
    byte*     Strength = MbQ->strength_hor[edge];
    int       StrValue;
    slice_t*  currSlice = MbQ->p_Slice;
    BlockPos* PicPos = MbQ->p_Vid->PicPos;

    int yQ = (edge < BLOCK_SIZE ? (edge << 2) - 1: 0);
    mb_t* neighbor = get_non_aff_neighbor_luma(MbQ, 0, yQ);
    mb_t* MbP = edge ? MbQ : neighbor;

    if (currSlice->slice_type == SP_SLICE || currSlice->slice_type == SI_SLICE ||
        MbQ->is_intra_block || MbP->is_intra_block) {
        // Set strength to either 3 or 4 regardless of pixel position
        StrValue = (edge == 0 && !currSlice->field_pic_flag) ? 4 : 3;
        memset(Strength, (byte) StrValue, BLOCK_SIZE * sizeof(byte));
        return;
    }

    if (edge && currSlice->slice_type == P_SLICE && MbQ->mb_type == PSKIP) {
        memset(Strength, 0, BLOCK_SIZE * sizeof(byte));
        return;
    }
/*
    if (edge && (MbQ->mb_type == P16x16 || MbQ->mb_type == P8x16)) {
        for (int idx = 0; idx < BLOCK_SIZE; idx++) {
            int blkQ = (yQ + 1) + idx;
            int blkP = (get_y_luma(yQ) & 0xFFFC) + idx;
            if ((MbQ->cbp_blks[0] & (((uint64_t)1 << blkQ) | ((uint64_t)1 << blkP))) != 0)
                StrValue = 2;
            else
                StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
            *(int*)(Strength + idx) = StrValue;
        }
        return;
    }
*/
    BlockPos mb = PicPos[MbQ->mbAddrX];
    mb.x <<= 2;
    mb.y <<= 2;

    for (int idx = 0; idx < BLOCK_SIZE; idx++) {
        int blkQ = (yQ + 1) + idx;
        int blkP = (get_y_luma(yQ) & 0xFFFC) + idx;

        if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 || (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
            StrValue = 2;
        else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
            int blk_x  = mb.x + (blkQ  & 3);
            int blk_y  = mb.y + (blkQ >> 2);
            int blk_x2 = (get_pos_x_luma(neighbor,  0) >> 2) + idx;
            int blk_y2 = (get_pos_y_luma(neighbor, yQ) >> 2);
            pic_motion_params *mv_info_p = &p->mv_info[blk_y ][blk_x ];
            pic_motion_params *mv_info_q = &p->mv_info[blk_y2][blk_x2];
            storable_picture* ref_p0 = mv_info_p->ref_pic[LIST_0];
            storable_picture* ref_q0 = mv_info_q->ref_pic[LIST_0];
            storable_picture* ref_p1 = mv_info_p->ref_pic[LIST_1];
            storable_picture* ref_q1 = mv_info_q->ref_pic[LIST_1];            

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

/*!
 *********************************************************************************************
 * \brief
 *    returns a buffer of 16 Strength values for one stripe in a mb (for MBAFF)
 *********************************************************************************************
 */
void get_strength_ver_MBAff(byte *Strength, mb_t *MbQ, int edge, int mvlimit, storable_picture *p)
{
  //byte *Strength = MbQ->strength_ver[edge];
    short  blkP, blkQ, idx;

    int    StrValue;
    short  mb_x, mb_y;

    mb_t *MbP;

    PixelPos pixP;
    VideoParameters *p_Vid = MbQ->p_Vid;
    slice_t *currSlice = MbQ->p_Slice;
    BlockPos *PicPos = p_Vid->PicPos;

    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    if (p->slice_type == SP_SLICE || p->slice_type == SI_SLICE) {
        for (idx = 0; idx < MB_BLOCK_SIZE; ++idx) {
            getAffNeighbour(MbQ, edge - 1, idx, mb_size, &pixP);
            blkQ = (short) ((idx & 0xFFFC) + (edge >> 2));
            blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

            MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag);

            Strength[idx] = (edge == 0) ? 4 : 3;
        }
    } else {
        getAffNeighbour(MbQ, edge - 1, 0, mb_size, &pixP);

        MbP = &(p_Vid->mb_data[pixP.mb_addr]);
        // Neighboring Frame MBs
        if (!MbQ->mb_field_decoding_flag && !MbP->mb_field_decoding_flag) {
            MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag); 
            if (MbQ->is_intra_block || MbP->is_intra_block) {
                //printf("idx %d %d %d %d %d\n", idx, pixP.x, pixP.y, pixP.pos_x, pixP.pos_y);
                // Start with Strength=3. or Strength=4 for Mb-edge
                StrValue = (edge == 0) ? 4 : 3;
                memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
            } else {
                get_mb_block_pos_mbaff (PicPos, MbQ->mbAddrX, &mb_x, &mb_y);
                for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
                    blkQ = (short) ((idx & 0xFFFC) + (edge >> 2));
                    blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

                    if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 || (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
                        StrValue = 2;
                    else if (edge && ((MbQ->mb_type == 1)  || (MbQ->mb_type == 2)))
                        StrValue = 0; // if internal edge of certain types, we already know StrValue should be 0
                    else { // for everything else, if no coefs, but vector difference >= 1 set Strength=1
                        int blk_y  = ((mb_y<<2) + (blkQ >> 2));
                        int blk_x  = ((mb_x<<2) + (blkQ  & 3));
                        int blk_y2 = (pixP.pos_y >> 2);
                        int blk_x2 = (pixP.pos_x >> 2);
                        pic_motion_params *mv_info_p = &p->mv_info[blk_y ][blk_x ];
                        pic_motion_params *mv_info_q = &p->mv_info[blk_y2][blk_x2];
                        storable_picture* ref_p0 = mv_info_p->ref_pic[LIST_0];
                        storable_picture* ref_q0 = mv_info_q->ref_pic[LIST_0];
                        storable_picture* ref_p1 = mv_info_p->ref_pic[LIST_1];
                        storable_picture* ref_q1 = mv_info_q->ref_pic[LIST_1];

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
                getAffNeighbour(MbQ, edge - 1, idx, mb_size, &pixP);
                blkQ = (short) ((idx & 0xFFFC) + (edge >> 2));
                blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

                MbP = &(p_Vid->mb_data[pixP.mb_addr]);
                MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag); 

                // Start with Strength=3. or Strength=4 for Mb-edge
                Strength[idx] = (edge == 0 && (((!p->mb_aff_frame_flag && !currSlice->field_pic_flag) ||
                    (p->mb_aff_frame_flag && !MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag)) ||
                    ((p->mb_aff_frame_flag || currSlice->field_pic_flag)))) ? 4 : 3;

                if (!MbQ->is_intra_block && !MbP->is_intra_block) {
                    if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 || (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
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
                            pic_motion_params *mv_info_p = &p->mv_info[blk_y ][blk_x ];
                            pic_motion_params *mv_info_q = &p->mv_info[blk_y2][blk_x2];
                            storable_picture* ref_p0 = mv_info_p->ref_pic[LIST_0];
                            storable_picture* ref_q0 = mv_info_q->ref_pic[LIST_0];
                            storable_picture* ref_p1 = mv_info_p->ref_pic[LIST_1];
                            storable_picture* ref_q1 = mv_info_q->ref_pic[LIST_1];

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
void get_strength_hor_MBAff(byte *Strength, mb_t *MbQ, int edge, int mvlimit, storable_picture *p)
{
    short  blkP, blkQ, idx;
    short  blk_x, blk_x2, blk_y, blk_y2 ;

    int    StrValue;
    int    xQ, yQ = (edge < MB_BLOCK_SIZE ? edge : 1);
    short  mb_x, mb_y;

    mb_t *MbP;

    PixelPos pixP;
    VideoParameters *p_Vid = MbQ->p_Vid;
    BlockPos *PicPos = p_Vid->PicPos;

    int mb_size[2] = { MB_BLOCK_SIZE, MB_BLOCK_SIZE };

    if (p->slice_type == SP_SLICE || p->slice_type == SI_SLICE) {
        for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
            xQ = idx;
            getAffNeighbour(MbQ, xQ, yQ - 1, mb_size, &pixP);

            blkQ = (short) ((yQ & 0xFFFC) + (xQ >> 2));
            blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

            MbP = &(p_Vid->mb_data[pixP.mb_addr]);
            MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag);

            StrValue = (edge == 0 && (!MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag)) ? 4 : 3;
      
            *(int*)(Strength+idx) = StrValue * 0x01010101;
        }
    } else {
        getAffNeighbour(MbQ, 0, yQ - 1, mb_size, &pixP);
        MbP = &(p_Vid->mb_data[pixP.mb_addr]);
        MbQ->mixedModeEdgeFlag = (MbQ->mb_field_decoding_flag != MbP->mb_field_decoding_flag); 

        // Set intra mode deblocking
        if (MbQ->is_intra_block || MbP->is_intra_block) {
            StrValue = (edge == 0 && (!MbP->mb_field_decoding_flag && !MbQ->mb_field_decoding_flag)) ? 4 : 3;
            memset(Strength, (byte) StrValue, MB_BLOCK_SIZE * sizeof(byte));
        } else {
            for (idx = 0; idx < MB_BLOCK_SIZE; idx += BLOCK_SIZE) {
                xQ = idx;    
                getAffNeighbour(MbQ, xQ, yQ - 1, mb_size, &pixP);

                blkQ = (short) ((yQ & 0xFFFC) + (xQ >> 2));
                blkP = (short) ((pixP.y & 0xFFFC) + (pixP.x >> 2));

                if ((MbQ->cbp_blks[0] & ((uint64_t)1 << blkQ)) != 0 || (MbP->cbp_blks[0] & ((uint64_t)1 << blkP)) != 0)
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
                        pic_motion_params *mv_info_p = &p->mv_info[blk_y ][blk_x ];
                        pic_motion_params *mv_info_q = &p->mv_info[blk_y2][blk_x2];
                        storable_picture* ref_p0 = mv_info_p->ref_pic[LIST_0];
                        storable_picture* ref_q0 = mv_info_q->ref_pic[LIST_0];
                        storable_picture* ref_p1 = mv_info_p->ref_pic[LIST_1];
                        storable_picture* ref_q1 = mv_info_q->ref_pic[LIST_1];

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

