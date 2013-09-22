#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "transform.h"
#include "inter_prediction.h"
#include "intra_prediction.h"
#include "neighbour.h"
#include "mv_prediction.h"


#define INVALIDINDEX  (-135792468)


static const uint8_t decode_block_scan[16] = { 0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15 };

static inline int RSD(int x)
{
    return ((x&2)?(x|1):(x&(~1)));
}

//! used to control block sizes : Not used/16x16/16x8/8x16/8x8/8x4/4x8/4x4
static const int BLOCK_STEP[8][2] = {
    {0, 0}, {4, 4}, {4, 2}, {2, 4},
    {2, 2}, {2, 1}, {1, 2}, {1, 1}
};

static void update_direct_mv_info_temporal(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    int step_h0 = BLOCK_STEP [currMB->mb_type][0];
    int step_v0 = BLOCK_STEP [currMB->mb_type][1];

    int j;
    int j6;
    int j4, i4;
    storable_picture *dec_picture = currSlice->dec_picture;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    bool has_direct = (currMB->SubMbType[0] == 0) | (currMB->SubMbType[1] == 0) |
                      (currMB->SubMbType[2] == 0) | (currMB->SubMbType[3] == 0);

    int block_y_aff;
    if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
        block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
    else
        block_y_aff = currMB->mb.y * 4;

    if (!has_direct)
        return;

    int mv_scale = 0;

    for (int k = 0; k < 4; ++k) { // Scan all blocks
        if (currMB->SubMbType[k] == 0) {
            currMB->SubMbPredMode[k] = 2;
            for (int j0 = 2 * (k >> 1); j0 < 2 * (k >> 1) + 2; j0 += step_v0) {
                for (int i0 = currMB->mb.x * 4 + 2 * (k & 0x01);
                     i0 < currMB->mb.x * 4 + 2 * (k & 0x01) + 2; i0 += step_h0) {
                    pic_motion_params *colocated;
                    int mapped_idx = -1, iref;

                    colocated = sps->direct_8x8_inference_flag ? 
                                &list1[0]->mv_info[RSD(block_y_aff + j0)][RSD(i0)] :
                                &list1[0]->mv_info[block_y_aff + j0][i0];
                    if (currSlice->MbaffFrameFlag) {
                        assert(sps->direct_8x8_inference_flag);
                        if (!currMB->mb_field_decoding_flag &&
                            ((currSlice->listX[LIST_1][0]->iCodingType == FRAME_MB_PAIR_CODING &&
                              currSlice->listX[LIST_1][0]->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
                             (currSlice->listX[LIST_1][0]->iCodingType == FIELD_CODING))) {
                            if (abs(dec_picture->poc - currSlice->listX[LIST_1+4][0]->poc) >
                                abs(dec_picture->poc - currSlice->listX[LIST_1+2][0]->poc)) {
                                colocated = sps->direct_8x8_inference_flag ? 
                                            &currSlice->listX[LIST_1+2][0]->mv_info[RSD(block_y_aff + j0)>>1][RSD(i0)] :
                                            &currSlice->listX[LIST_1+2][0]->mv_info[(block_y_aff + j0)>>1][i0];
                            } else {
                                colocated = sps->direct_8x8_inference_flag ? 
                                            &currSlice->listX[LIST_1+4][0]->mv_info[RSD(block_y_aff + j0)>>1][RSD(i0)] :
                                            &currSlice->listX[LIST_1+4][0]->mv_info[(block_y_aff + j0)>>1][i0];
                            }
                        }
                    } else if (!sps->frame_mbs_only_flag && !currSlice->field_pic_flag &&
                               currSlice->listX[LIST_1][0]->iCodingType != FRAME_CODING) {
                        if (abs(dec_picture->poc - list1[0]->bottom_field->poc) >
                            abs(dec_picture->poc - list1[0]->top_field->poc)) {
                            colocated = sps->direct_8x8_inference_flag ? 
                                        &list1[0]->top_field->mv_info[RSD(block_y_aff + j0)>>1][RSD(i0)] :
                                        &list1[0]->top_field->mv_info[(block_y_aff + j0)>>1][i0];
                        } else {
                            colocated = sps->direct_8x8_inference_flag ? 
                                        &list1[0]->bottom_field->mv_info[RSD(block_y_aff + j0)>>1][RSD(i0)] :
                                        &list1[0]->bottom_field->mv_info[(block_y_aff + j0)>>1][i0];
                        }
                    } else if (!sps->frame_mbs_only_flag &&
                                currSlice->field_pic_flag &&
                                currSlice->structure != list1[0]->structure &&
                                list1[0]->coded_frame) {
                        if (!currSlice->bottom_field_flag) {
                            colocated = sps->direct_8x8_inference_flag ? 
                                        &list1[0]->frame->top_field->mv_info[RSD(block_y_aff + j0)][RSD(i0)] :
                                        &list1[0]->frame->top_field->mv_info[block_y_aff + j0][i0];
                        } else {
                            colocated = sps->direct_8x8_inference_flag ? 
                                        &list1[0]->frame->bottom_field->mv_info[RSD(block_y_aff + j0)][RSD(i0)] :
                                        &list1[0]->frame->bottom_field->mv_info[block_y_aff + j0][i0];
                        }
                    }

                    int refList = colocated->ref_idx[LIST_0 ]== -1 ? LIST_1 : LIST_0;
                    int ref_idx = colocated->ref_idx[refList];

                    if (ref_idx == -1) {
                        for (j4 = currMB->mb.y * 4 + j0; j4 < currMB->mb.y * 4 + j0 + step_v0; ++j4) {
                            for (i4 = i0; i4 < i0 + step_h0; ++i4) {
                                pic_motion_params *mv_info = &dec_picture->mv_info[j4][i4];
                                mv_info->ref_pic[LIST_0] = list0[0];
                                mv_info->ref_pic[LIST_1] = list1[0];
                                mv_info->mv     [LIST_0] = zero_mv;
                                mv_info->mv     [LIST_1] = zero_mv;
                                mv_info->ref_idx[LIST_0] = 0;
                                mv_info->ref_idx[LIST_1] = 0;
                            }
                        }
                    } else {
                        if ((currSlice->MbaffFrameFlag &&
                             ((currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure==FRAME) || 
                             (!currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure!=FRAME))) ||
                            (!currSlice->MbaffFrameFlag &&
                             ((currSlice->field_pic_flag==0 && colocated->ref_pic[refList]->structure!=FRAME) ||
                             (currSlice->field_pic_flag==1 && colocated->ref_pic[refList]->structure==FRAME)))) {
                            //! Frame with field co-located
                            for (iref = 0; iref < min<int>(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]); ++iref) {
                                if (currSlice->listX[LIST_0 + list_offset][iref]->top_field == colocated->ref_pic[refList] ||
                                    currSlice->listX[LIST_0 + list_offset][iref]->bottom_field == colocated->ref_pic[refList] ||
                                    currSlice->listX[LIST_0 + list_offset][iref]->frame == colocated->ref_pic[refList]) {
                                    if (currSlice->field_pic_flag &&
                                        (currSlice->listX[LIST_0 + list_offset][iref]->structure != currSlice->structure)) {
                                        mapped_idx = INVALIDINDEX;
                                    } else {
                                        mapped_idx = iref;
                                        break;
                                    }
                                } else //! invalid index. Default to zero even though this case should not happen
                                    mapped_idx = INVALIDINDEX;
                            }
                        } else {
                            for (iref = 0; iref < min<int>(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]); ++iref) {
                                if (currSlice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                                    mapped_idx=iref;
                                    break;
                                } else //! invalid index. Default to zero even though this case should not happen
                                    mapped_idx = INVALIDINDEX;
                            }
                        }

                        if (mapped_idx != INVALIDINDEX) {
                            for (j = j0; j < j0 + step_v0; ++j) {
                                j4 = currMB->mb.y * 4 + j;
                                j6 = block_y_aff + j;

                                for (i4 = i0; i4 < i0 + step_h0; ++i4) {
                                    pic_motion_params *colocated = sps->direct_8x8_inference_flag ? 
                                        &list1[0]->mv_info[RSD(j6)][RSD(i4)] : &list1[0]->mv_info[j6][i4];
                                    pic_motion_params *mv_info = &dec_picture->mv_info[j4][i4];
                                    int mv_y;
                                    if (currSlice->MbaffFrameFlag) {
                                        if (!currMB->mb_field_decoding_flag &&
                                            ((currSlice->listX[LIST_1][0]->iCodingType==FRAME_MB_PAIR_CODING &&
                                              currSlice->listX[LIST_1][0]->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
                                            (currSlice->listX[LIST_1][0]->iCodingType==FIELD_CODING))) {
                                            if (abs(dec_picture->poc - currSlice->listX[LIST_1+4][0]->poc) >
                                                abs(dec_picture->poc - currSlice->listX[LIST_1+2][0]->poc)) {
                                                colocated = sps->direct_8x8_inference_flag ? 
                                                            &currSlice->listX[LIST_1+2][0]->mv_info[RSD(j6)>>1][RSD(i4)] :
                                                            &currSlice->listX[LIST_1+2][0]->mv_info[j6>>1][i4];
                                            } else {
                                                colocated = sps->direct_8x8_inference_flag ? 
                                                            &currSlice->listX[LIST_1+4][0]->mv_info[RSD(j6)>>1][RSD(i4)] :
                                                            &currSlice->listX[LIST_1+4][0]->mv_info[j6>>1][i4];
                                            }
                                        }
                                    } else if (!sps->frame_mbs_only_flag &&
                                               !currSlice->field_pic_flag &&
                                                currSlice->listX[LIST_1][0]->iCodingType != FRAME_CODING) {
                                        if (abs(dec_picture->poc - list1[0]->bottom_field->poc) >
                                            abs(dec_picture->poc - list1[0]->top_field->poc)) {
                                            colocated = sps->direct_8x8_inference_flag ? 
                                                        &list1[0]->top_field->mv_info[RSD(j6)>>1][RSD(i4)] :
                                                        &list1[0]->top_field->mv_info[(j6)>>1][i4];
                                        } else {
                                            colocated = sps->direct_8x8_inference_flag ? 
                                                        &list1[0]->bottom_field->mv_info[RSD(j6)>>1][RSD(i4)] :
                                                        &list1[0]->bottom_field->mv_info[(j6)>>1][i4];
                                        }
                                    } else if (!sps->frame_mbs_only_flag &&
                                               currSlice->field_pic_flag &&
                                               currSlice->structure != list1[0]->structure &&
                                               list1[0]->coded_frame) {
                                        if (!currSlice->bottom_field_flag) {
                                            colocated = sps->direct_8x8_inference_flag ? 
                                                        &list1[0]->frame->top_field->mv_info[RSD(j6)][RSD(i4)] :
                                                        &list1[0]->frame->top_field->mv_info[j6][i4];
                                        } else {
                                            colocated = sps->direct_8x8_inference_flag ? 
                                                        &list1[0]->frame->bottom_field->mv_info[RSD(j6)][RSD(i4)] :
                                                        &list1[0]->frame->bottom_field->mv_info[j6][i4];
                                        }
                                    }

                                    mv_y = colocated->mv[refList].mv_y; 
                                    if ((currSlice->MbaffFrameFlag &&
                                         !currMB->mb_field_decoding_flag &&
                                         colocated->ref_pic[refList]->structure!=FRAME) ||
                                        (!currSlice->MbaffFrameFlag &&
                                         currSlice->field_pic_flag == 0 &&
                                         colocated->ref_pic[refList]->structure!=FRAME))
                                        mv_y *= 2;
                                    else if ((currSlice->MbaffFrameFlag &&
                                              currMB->mb_field_decoding_flag &&
                                              colocated->ref_pic[refList]->structure==FRAME) ||
                                             (!currSlice->MbaffFrameFlag &&
                                              currSlice->field_pic_flag == 1 &&
                                              colocated->ref_pic[refList]->structure==FRAME))
                                        mv_y /= 2;

                                    mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];

                                    mv_info->ref_idx[LIST_0] = (char) mapped_idx;
                                    mv_info->ref_idx[LIST_1] = 0;
                                    mv_info->ref_pic[LIST_0] = list0[mapped_idx];
                                    mv_info->ref_pic[LIST_1] = list1[0];

                                    if (mv_scale == 9999 || currSlice->listX[LIST_0+list_offset][mapped_idx]->is_long_term) {
                                        mv_info->mv[LIST_0].mv_x = colocated->mv[refList].mv_x;
                                        mv_info->mv[LIST_0].mv_y = (short) mv_y;
                                        mv_info->mv[LIST_1] = zero_mv;
                                    } else {
                                        mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * colocated->mv[refList].mv_x + 128 ) >> 8);
                                        mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * mv_y + 128 ) >> 8);
                                        mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - colocated->mv[refList].mv_x);
                                        mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - mv_y);
                                    }
                                }
                            }
                        } else if (INVALIDINDEX == mapped_idx) {
                            error("temporal direct error: colocated block has ref that is unavailable",-1111);
                        }
                    }
                }
            }
        }
    }
}


static inline void update_neighbor_mvs(pic_motion_params **motion, const pic_motion_params *mv_info, int i4)
{
    (*motion++)[i4 + 1] = *mv_info;
    (*motion  )[i4    ] = *mv_info;
    (*motion  )[i4 + 1] = *mv_info;
}


static int get_colocated_info_4x4(mb_t *currMB, storable_picture *list1, int i, int j)
{
    if (list1->is_long_term)
        return 1;

    pic_motion_params *fs = &list1->mv_info[j][i];

    int moving =
        !((fs->ref_idx[LIST_0] == 0 &&
           abs(fs->mv[LIST_0].mv_x) >> 1 == 0 && abs(fs->mv[LIST_0].mv_y) >> 1 == 0) ||
          (fs->ref_idx[LIST_0] == -1 && fs->ref_idx[LIST_1] == 0 &&
           abs(fs->mv[LIST_1].mv_x) >> 1 == 0 && abs(fs->mv[LIST_1].mv_y) >> 1 == 0));
    return moving;  
}

static int get_colocated_info_8x8(mb_t *currMB, storable_picture *list1, int i, int j)
{
    if (list1->is_long_term)
        return 1;

    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;

    int jj = RSD(j);
    int ii = RSD(i);
    pic_motion_params *fs = &list1->mv_info[jj][ii];

    if (currSlice->MbaffFrameFlag ||
        (!sps->frame_mbs_only_flag &&
         ((!currSlice->field_pic_flag && list1->iCodingType == FIELD_CODING) ||
          (currSlice->structure != list1->structure && list1->coded_frame)))) {
        if (currSlice->field_pic_flag && currSlice->structure != list1->structure && list1->coded_frame) {
            if (!currSlice->bottom_field_flag)
                fs = &list1->top_field->mv_info[jj][ii];
            else
                fs = &list1->bottom_field->mv_info[jj][ii];
        } else if ((!currMB->mb_field_decoding_flag &&
                    (list1->iCodingType == FIELD_CODING || list1->motion.mb_field_decoding_flag[currMB->mbAddrX]))) {
            if (abs(currSlice->dec_picture->poc - list1->bottom_field->poc) >
                abs(currSlice->dec_picture->poc - list1->top_field->poc))
                fs = &list1->top_field->mv_info[jj >> 1][ii];
            else
                fs = &list1->bottom_field->mv_info[jj >> 1][ii];
        }
    } else if (sps->separate_colour_plane_flag && sps->chroma_format_idc == YUV444)
        fs = &list1->JVmv_info[currMB->p_Slice->colour_plane_id][jj][ii];

    int moving =
        !((fs->ref_idx[LIST_0] == 0 &&
           abs(fs->mv[LIST_0].mv_x) >> 1 == 0 && abs(fs->mv[LIST_0].mv_y) >> 1 == 0) ||
          (fs->ref_idx[LIST_0] == -1 && fs->ref_idx[LIST_1] == 0 &&
           abs(fs->mv[LIST_1].mv_x) >> 1 == 0 && abs(fs->mv[LIST_1].mv_y) >> 1 == 0));
    return moving;  
}


static void update_direct_mv_info_spatial_8x8(mb_t *currMB)
{
    bool has_direct = (currMB->SubMbType[0] == 0) | (currMB->SubMbType[1] == 0) |
                      (currMB->SubMbType[2] == 0) | (currMB->SubMbType[3] == 0);
    if (!has_direct)
        return;

    slice_t *currSlice = currMB->p_Slice;
    storable_picture *dec_picture = currSlice->dec_picture;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    char l0_rFrame, l1_rFrame;
    MotionVector pmvl0, pmvl1;
    currSlice->inter_prediction.prepare_direct_params(currMB, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    storable_picture *ref_pic_l[2];
    char             ref_idx_l[2];
    if (l0_rFrame < 0 && l1_rFrame < 0) {
        ref_pic_l[LIST_0] = list0[0];
        ref_pic_l[LIST_1] = list1[0];
        ref_idx_l[LIST_0] = 0;
        ref_idx_l[LIST_1] = 0;
    } else {
        ref_pic_l[LIST_0] = l0_rFrame == -1 ? NULL : list0[(short)l0_rFrame];
        ref_pic_l[LIST_1] = l1_rFrame == -1 ? NULL : list1[(short)l1_rFrame];
        ref_idx_l[LIST_0] = l0_rFrame;
        ref_idx_l[LIST_1] = l1_rFrame;
    }

    for (int block8x8 = 0; block8x8 < 4; block8x8++) {
        if (currMB->SubMbType[block8x8] == 0) {
            int i = (block8x8 % 2) * 2;
            int j = (block8x8 / 2) * 2;

            pic_motion_params *mv_info = &dec_picture->mv_info[currMB->mb.y * 4 + j][currMB->mb.x * 4 + i];
            mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
            mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
            mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
            mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

            if (l0_rFrame < 0 && l1_rFrame < 0) {
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = zero_mv;
            } else {
                int block_y_aff;
                if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
                    block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
                else
                    block_y_aff = currMB->mb.y * 4;
                bool is_not_moving = (get_colocated_info_8x8(currMB, list1[0], currMB->mb.x * 4 + i, block_y_aff + j) == 0);
                mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? zero_mv : pmvl0;
                mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? zero_mv : pmvl1;
            }

            update_neighbor_mvs(&dec_picture->mv_info[currMB->mb.y * 4 + j], mv_info, currMB->mb.x * 4 + i);              
        }
    }
}

static void update_direct_mv_info_spatial_4x4(mb_t *currMB)
{
    bool has_direct = (currMB->SubMbType[0] == 0) | (currMB->SubMbType[1] == 0) |
                      (currMB->SubMbType[2] == 0) | (currMB->SubMbType[3] == 0);
    if (!has_direct)
        return;

    slice_t *currSlice = currMB->p_Slice;
    storable_picture *dec_picture = currMB->p_Vid->dec_picture;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    char l0_rFrame, l1_rFrame;
    MotionVector pmvl0, pmvl1;
    currSlice->inter_prediction.prepare_direct_params(currMB, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    storable_picture *ref_pic_l[2];
    char              ref_idx_l[2];
    if (l0_rFrame < 0 && l1_rFrame < 0) {
        ref_pic_l[LIST_0] = list0[0];
        ref_pic_l[LIST_1] = list1[0];
        ref_idx_l[LIST_0] = 0;
        ref_idx_l[LIST_1] = 0;
    } else {
        ref_pic_l[LIST_0] = l0_rFrame == -1 ? NULL : list0[(short)l0_rFrame];
        ref_pic_l[LIST_1] = l1_rFrame == -1 ? NULL : list1[(short)l1_rFrame];
        ref_idx_l[LIST_0] = l0_rFrame;
        ref_idx_l[LIST_1] = l1_rFrame;
    }

    for (int block4x4 = 0; block4x4 < 16; block4x4++) {
        if (currMB->SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);

            pic_motion_params *mv_info = &dec_picture->mv_info[currMB->mb.y * 4 + j][currMB->mb.x * 4 + i];
            mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
            mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
            mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
            mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

            if (l0_rFrame < 0 && l1_rFrame < 0) {
                mv_info->mv[LIST_0] = zero_mv;
                mv_info->mv[LIST_1] = zero_mv;
            } else {
                int block_y_aff;
                if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
                    block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
                else
                    block_y_aff = currMB->mb.y * 4;
                bool is_not_moving = (get_colocated_info_4x4(currMB, list1[0], currMB->mb.x * 4 + i, block_y_aff + j) == 0);
                mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? zero_mv : pmvl0;
                mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? zero_mv : pmvl1;
            }
        }
    }
}


void inter_prediction_t::update_direct_mv_info(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    sps_t *sps = currSlice->active_sps;
    if (!currSlice->direct_spatial_mv_pred_flag)
        update_direct_mv_info_temporal(currMB);
    else if (sps->direct_8x8_inference_flag)
        update_direct_mv_info_spatial_8x8(currMB);
    else
        update_direct_mv_info_spatial_4x4(currMB);
}


void inter_prediction_t::get_direct8x8temporal(mb_t *currMB, int block8x8)
{
    slice_t *currSlice = currMB->p_Slice;
    storable_picture* dec_picture = currSlice->dec_picture;
    sps_t *sps = currSlice->active_sps;
  
    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    int block_y_aff;
    if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
        block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
    else
        block_y_aff = currMB->mb.y * 4;

    for (int block4x4 = block8x8 * 4; block4x4 < block8x8 * 4 + 4; block4x4++) {
        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
        int i4 = currMB->mb.x * 4 + i;
        int j4 = currMB->mb.y * 4 + j;
        int j6 = block_y_aff + j;
        pic_motion_params *mv_info = &dec_picture->mv_info[j4][i4];
        pic_motion_params *colocated = &list1[0]->mv_info[RSD(j6)][RSD(i4)];

        if (sps->separate_colour_plane_flag && sps->chroma_format_idc == YUV444)
            colocated = &list1[0]->JVmv_info[currSlice->colour_plane_id][RSD(j6)][RSD(i4)];

        if (currSlice->MbaffFrameFlag) {
            if (!currMB->mb_field_decoding_flag &&
                ((currSlice->listX[LIST_1][0]->iCodingType == FRAME_MB_PAIR_CODING &&
                  currSlice->listX[LIST_1][0]->motion.mb_field_decoding_flag[currMB->mbAddrX]) ||
                 currSlice->listX[LIST_1][0]->iCodingType == FIELD_CODING)) {
                if (abs(dec_picture->poc - currSlice->listX[LIST_1+4][0]->poc) >
                    abs(dec_picture->poc - currSlice->listX[LIST_1+2][0]->poc))
                    colocated = &currSlice->listX[LIST_1+2][0]->mv_info[RSD(j6)>>1][RSD(i4)];
                else
                    colocated = &currSlice->listX[LIST_1+4][0]->mv_info[RSD(j6)>>1][RSD(i4)];
            }
        } else if (!sps->frame_mbs_only_flag && !currSlice->field_pic_flag &&
                   currSlice->listX[LIST_1][0]->iCodingType != FRAME_CODING) {
            if (abs(dec_picture->poc - list1[0]->bottom_field->poc) >
                abs(dec_picture->poc - list1[0]->top_field->poc) )
                colocated = &list1[0]->top_field->mv_info[RSD(j6)>>1][RSD(i4)];
            else
                colocated = &list1[0]->bottom_field->mv_info[RSD(j6)>>1][RSD(i4)];
        } else if (!sps->frame_mbs_only_flag && currSlice->field_pic_flag &&
                   currSlice->structure != list1[0]->structure && list1[0]->coded_frame) {
            if (!currSlice->bottom_field_flag)
                colocated = &list1[0]->frame->top_field->mv_info[RSD(j6)][RSD(i4)];
            else
                colocated = &list1[0]->frame->bottom_field->mv_info[RSD(j6)][RSD(i4)];
        }

        int   refList = colocated->ref_idx[LIST_0] == -1 ? LIST_1 : LIST_0;
        short ref_idx = colocated->ref_idx[refList];

        if (ref_idx == -1) { // co-located is intra mode
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;

            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = 0;
        } else { // co-located skip or inter mode
            int mapped_idx = 0;
            int iref;
            if ((currSlice->MbaffFrameFlag &&
                 ((currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure == FRAME) || 
                  (!currMB->mb_field_decoding_flag && colocated->ref_pic[refList]->structure != FRAME))) ||
                (!currSlice->MbaffFrameFlag &&
                 ((currSlice->field_pic_flag == 0 && colocated->ref_pic[refList]->structure != FRAME) ||
                  (currSlice->field_pic_flag == 1 && colocated->ref_pic[refList]->structure == FRAME)))) {
                for (iref = 0; iref < min<int>(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]);iref++) {
                    if (currSlice->listX[LIST_0 + list_offset][iref]->top_field == colocated->ref_pic[refList] || 
                        currSlice->listX[LIST_0 + list_offset][iref]->bottom_field == colocated->ref_pic[refList] ||
                        currSlice->listX[LIST_0 + list_offset][iref]->frame == colocated->ref_pic[refList]) {
                        if (currSlice->field_pic_flag == 1 &&
                            currSlice->listX[LIST_0 + list_offset][iref]->structure != currSlice->structure)
                            mapped_idx = INVALIDINDEX;
                        else {
                            mapped_idx = iref;
                            break;
                        }
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx = INVALIDINDEX;
                }
            } else {
                for (iref = 0; iref < min<int>(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]);iref++) {
                    if (currSlice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                        mapped_idx = iref;            
                        break;
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx = INVALIDINDEX;
                }
            }

            if (INVALIDINDEX != mapped_idx) {
                int mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];
                int mv_y = colocated->mv[refList].mv_y; 
                if ((currSlice->MbaffFrameFlag && !currMB->mb_field_decoding_flag &&
                     colocated->ref_pic[refList]->structure!=FRAME) ||
                    (!currSlice->MbaffFrameFlag && currSlice->field_pic_flag == 0 &&
                     colocated->ref_pic[refList]->structure != FRAME))
                    mv_y *= 2;
                else if ((currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag &&
                          colocated->ref_pic[refList]->structure==FRAME) ||
                         (!currSlice->MbaffFrameFlag && currSlice->field_pic_flag == 1 &&
                          colocated->ref_pic[refList]->structure==FRAME))
                    mv_y /= 2;

                //! In such case, an array is needed for each different reference.
                if (mv_scale == 9999 || currSlice->listX[LIST_0 + list_offset][mapped_idx]->is_long_term) {
                    mv_info->mv[LIST_0].mv_x = colocated->mv[refList].mv_x;
                    mv_info->mv[LIST_0].mv_y = (short) mv_y;
                    mv_info->mv[LIST_1] = zero_mv;
                } else {
                    mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * colocated->mv[refList].mv_x + 128 ) >> 8);
                    mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * mv_y/*colocated->mv[refList].mv_y*/ + 128 ) >> 8);

                    mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - colocated->mv[refList].mv_x);
                    mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - mv_y/*colocated->mv[refList].mv_y*/);
                }

                mv_info->ref_idx[LIST_0] = (char) mapped_idx; //colocated->ref_idx[refList];
                mv_info->ref_idx[LIST_1] = 0;
            } else if (INVALIDINDEX == mapped_idx)
                error("temporal direct error: colocated block has ref that is unavailable",-1111);
        }
        // store reference picture ID determined by direct mode
        mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
        mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
    }
}

void inter_prediction_t::get_direct4x4temporal(mb_t *currMB, int block8x8)
{
    slice_t *currSlice = currMB->p_Slice;
    storable_picture* dec_picture = currSlice->dec_picture;
    sps_t *sps = currSlice->active_sps;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    int block_y_aff;
    if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
        block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
    else
        block_y_aff = currMB->mb.y * 4;

    for (int block4x4 = block8x8 * 4; block4x4 < block8x8 * 4 + 4; block4x4++) {
        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
        int i4 = currMB->mb.x * 4 + i;
        int j4 = currMB->mb.y * 4 + j;
        int j6 = block_y_aff + j;
        pic_motion_params *mv_info = &dec_picture->mv_info[j4][i4];
        pic_motion_params *colocated = &list1[0]->mv_info[j6][i4];
        if (sps->separate_colour_plane_flag && sps->chroma_format_idc == YUV444)
            colocated = &list1[0]->JVmv_info[currMB->p_Slice->colour_plane_id][RSD(j6)][RSD(i4)];

        int   refList = colocated->ref_idx[LIST_0] == -1 ? LIST_1 : LIST_0;
        short ref_idx = colocated->ref_idx[refList];

        if (ref_idx == -1) { // co-located is intra mode
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;

            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = 0;
        } else { // co-located skip or inter mode
            int mapped_idx = 0;
            int iref;

            for (iref = 0; iref < min<int>(currSlice->num_ref_idx_l0_active_minus1+1, currSlice->listXsize[LIST_0 + list_offset]); iref++) {
                if (currSlice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                    mapped_idx = iref;
                    break;
                } else //! invalid index. Default to zero even though this case should not happen
                    mapped_idx = INVALIDINDEX;
            }
            if (INVALIDINDEX == mapped_idx)
                error("temporal direct error: colocated block has ref that is unavailable",-1111);

            int mv_scale = currSlice->mvscale[LIST_0 + list_offset][mapped_idx];

            //! In such case, an array is needed for each different reference.
            if (mv_scale == 9999 || currSlice->listX[LIST_0 + list_offset][mapped_idx]->is_long_term) {
                mv_info->mv[LIST_0] = colocated->mv[refList];
                mv_info->mv[LIST_1] = zero_mv;
            } else {
                mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * colocated->mv[refList].mv_x + 128 ) >> 8);
                mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * colocated->mv[refList].mv_y + 128 ) >> 8);

                mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - colocated->mv[refList].mv_x);
                mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - colocated->mv[refList].mv_y);
            }

            mv_info->ref_idx[LIST_0] = (char) mapped_idx;
            mv_info->ref_idx[LIST_1] = 0;
        }
        // store reference picture ID determined by direct mode
        mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
        mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
    }
}

void inter_prediction_t::get_direct8x8spatial(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    storable_picture* dec_picture = currSlice->dec_picture;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    prepare_direct_params(currMB, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    int pred_dir = 0;
    if (l0_rFrame < 0 && l1_rFrame < 0)
        pred_dir = 2;
    else
        pred_dir = l1_rFrame == -1 ? 0 : l0_rFrame == -1 ? 1 : 2;

    storable_picture *ref_pic_l[2];
    char             ref_idx_l[2];
    if (l0_rFrame < 0 && l1_rFrame < 0) {
        ref_pic_l[LIST_0] = list0[0];
        ref_pic_l[LIST_1] = list1[0];
        ref_idx_l[LIST_0] = 0;
        ref_idx_l[LIST_1] = 0;
    } else {
        ref_pic_l[LIST_0] = l0_rFrame == -1 ? NULL : list0[(short)l0_rFrame];
        ref_pic_l[LIST_1] = l1_rFrame == -1 ? NULL : list1[(short)l1_rFrame];
        ref_idx_l[LIST_0] = l0_rFrame;
        ref_idx_l[LIST_1] = l1_rFrame;
    }

    for (int block8x8 = 0; block8x8 < 4; block8x8++) {
        int i = (block8x8 % 2) * 2;
        int j = (block8x8 / 2) * 2;
        currMB->SubMbPredMode[block8x8] = pred_dir;

        pic_motion_params *mv_info = &dec_picture->mv_info[currMB->mb.y * 4 + j][currMB->mb.x * 4 + i];
        mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
        mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
        mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
        mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

        if (l0_rFrame == 0 || l1_rFrame == 0) {
            int block_y_aff;
            if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
                block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
            else
                block_y_aff = currMB->mb.y * 4;
            bool is_not_moving = (get_colocated_info_8x8(currMB, list1[0], currMB->mb.x * 4 + i, block_y_aff + j) == 0);
            mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? zero_mv : pmvl0;
            mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? zero_mv : pmvl1;
        } else if (l0_rFrame < 0 && l1_rFrame < 0) {
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;
        } else {
            mv_info->mv[LIST_0] = l0_rFrame == -1 ? zero_mv : pmvl0;
            mv_info->mv[LIST_1] = l1_rFrame == -1 ? zero_mv : pmvl1;
        }

        update_neighbor_mvs(&dec_picture->mv_info[currMB->mb.y * 4 + j], mv_info, currMB->mb.x * 4 + i);
    }
}

void inter_prediction_t::get_direct4x4spatial(mb_t *currMB)
{
    slice_t *currSlice = currMB->p_Slice;
    storable_picture* dec_picture = currSlice->dec_picture;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    char l0_rFrame = -1, l1_rFrame = -1;
    MotionVector pmvl0 = zero_mv, pmvl1 = zero_mv;
    prepare_direct_params(currMB, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    int pred_dir = 0;
    if (l0_rFrame < 0 && l1_rFrame < 0)
        pred_dir = 2;
    else
        pred_dir = l1_rFrame == -1 ? 0 : l0_rFrame == -1 ? 1 : 2;

    storable_picture *ref_pic_l[2];
    char             ref_idx_l[2];
    if (l0_rFrame < 0 && l1_rFrame < 0) {
        ref_pic_l[LIST_0] = list0[0];
        ref_pic_l[LIST_1] = list1[0];
        ref_idx_l[LIST_0] = 0;
        ref_idx_l[LIST_1] = 0;
    } else {
        ref_pic_l[LIST_0] = l0_rFrame == -1 ? NULL : list0[(short)l0_rFrame];
        ref_pic_l[LIST_1] = l1_rFrame == -1 ? NULL : list1[(short)l1_rFrame];
        ref_idx_l[LIST_0] = l0_rFrame;
        ref_idx_l[LIST_1] = l1_rFrame;
    }

    for (int block4x4 = 0; block4x4 < 16; block4x4++) {
        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
        currMB->SubMbPredMode[block4x4 / 4] = pred_dir;

        pic_motion_params *mv_info = &dec_picture->mv_info[currMB->mb.y * 4 + j][currMB->mb.x * 4 + i];
        mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
        mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
        mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
        mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

        if (l0_rFrame < 0 && l1_rFrame < 0) {
            mv_info->mv[LIST_0] = zero_mv;
            mv_info->mv[LIST_1] = zero_mv;
        } else {
            int block_y_aff;
            if (currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag)
                block_y_aff = (currMB->mb.y * 4 - 4 * (currMB->mbAddrX % 2)) / 2;
            else
                block_y_aff = currMB->mb.y * 4;
            bool is_not_moving = (get_colocated_info_4x4(currMB, list1[0], currMB->mb.x * 4 + i, block_y_aff + j) == 0);
            mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? zero_mv : pmvl0;
            mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? zero_mv : pmvl1;
        }
    }
}

int inter_prediction_t::get_inter8x8(mb_t *currMB, int block8x8)
{
    slice_t *currSlice = currMB->p_Slice;
    storable_picture* dec_picture = currSlice->dec_picture;

    int list_offset = currSlice->MbaffFrameFlag && currMB->mb_field_decoding_flag ?
                      currMB->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture **list0 = currSlice->listX[LIST_0 + list_offset];
    storable_picture **list1 = currSlice->listX[LIST_1 + list_offset];

    int mv_mode  = currMB->SubMbType    [block8x8];
    int pred_dir = currMB->SubMbPredMode[block8x8];

    if (mv_mode == 0) {
        int k_start = (block8x8 << 2);

        for (int k = k_start; k < k_start + BLOCK_MULTIPLE; k ++) {
            int i  =  (decode_block_scan[k] & 3);
            int j  = ((decode_block_scan[k] >> 2) & 3);
            int i4 = currMB->mb.x * 4 + i;
            int j4 = currMB->mb.y * 4 + j;
            pic_motion_params *mv_info = &dec_picture->mv_info[j4][i4];

            if (currSlice->direct_spatial_mv_pred_flag) {
                if (mv_info->ref_idx[LIST_1] == -1)
                    pred_dir = 0;
                else if (mv_info->ref_idx[LIST_0] == -1)
                    pred_dir = 1;
                else
                    pred_dir = 2;
            } else {
                mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
                mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
            }
        }
    }

    return pred_dir;
}


static inline void set_direct_references(mb_t *currMB, const PixelPos *mb, char *l0_rFrame, char *l1_rFrame, pic_motion_params **mv_info)
{
    slice_t *currSlice = currMB->p_Slice;
    mb_t *mb_data = currMB->p_Vid->mb_data;

    if (!mb->available) {
        *l0_rFrame = -1;
        *l1_rFrame = -1;
        return;
    }

    char *ref_idx = mv_info[mb->pos_y][mb->pos_x].ref_idx;
    if (!currSlice->MbaffFrameFlag ||
        (currMB->mb_field_decoding_flag == mb_data[mb->mb_addr].mb_field_decoding_flag)) {
        *l0_rFrame = ref_idx[LIST_0];
        *l1_rFrame = ref_idx[LIST_1];
    } else if (currMB->mb_field_decoding_flag) {
        *l0_rFrame = (ref_idx[LIST_0] < 0) ? ref_idx[LIST_0] : ref_idx[LIST_0] * 2;
        *l1_rFrame = (ref_idx[LIST_1] < 0) ? ref_idx[LIST_1] : ref_idx[LIST_1] * 2;
    } else {
        *l0_rFrame = (ref_idx[LIST_0] >> 1);
        *l1_rFrame = (ref_idx[LIST_1] >> 1);
    }
}

void inter_prediction_t::prepare_direct_params(mb_t *currMB, MotionVector *pmvl0, MotionVector *pmvl1, char *l0_rFrame, char *l1_rFrame)
{
    slice_t *currSlice = currMB->p_Slice;
    storable_picture* dec_picture = currSlice->dec_picture;

    PixelPos mb[4];
    get_neighbors(currMB, mb, 0, 0, 16);

    char l0_refA, l0_refB, l0_refC;
    char l1_refA, l1_refB, l1_refC;
    pic_motion_params **mv_info = dec_picture->mv_info;
    set_direct_references(currMB, &mb[0], &l0_refA, &l1_refA, mv_info);
    set_direct_references(currMB, &mb[1], &l0_refB, &l1_refB, mv_info);
    set_direct_references(currMB, &mb[2], &l0_refC, &l1_refC, mv_info);

    *l0_rFrame = (char) min(min((unsigned char) l0_refA, (unsigned char) l0_refB), (unsigned char) l0_refC);
    *l1_rFrame = (char) min(min((unsigned char) l1_refA, (unsigned char) l1_refB), (unsigned char) l1_refC);

    if (*l0_rFrame >=0)
        GetMVPredictor(currMB, mb, pmvl0, *l0_rFrame, mv_info, LIST_0, 0, 0, 16, 16);
    if (*l1_rFrame >=0)
        GetMVPredictor(currMB, mb, pmvl1, *l1_rFrame, mv_info, LIST_1, 0, 0, 16, 16);
}
