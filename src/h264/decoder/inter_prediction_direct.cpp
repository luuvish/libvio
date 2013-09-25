#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "inter_prediction.h"
#include "neighbour.h"
#include "mv_prediction.h"


#define INVALIDINDEX  (-135792468)


static inline int RSD(int x)
{
    return (x & 2) ? (x | 1) : (x & ~1);
}


void inter_prediction_t::get_direct_temporal(mb_t* mb, bool dir)
{
    bool has_direct = (mb->SubMbType[0] == 0) | (mb->SubMbType[1] == 0) |
                      (mb->SubMbType[2] == 0) | (mb->SubMbType[3] == 0);
    if (!has_direct)
        return;

    slice_t* slice = mb->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;
    sps_t* sps = slice->active_sps;

    int list_offset = slice->MbaffFrameFlag && mb->mb_field_decoding_flag ?
                      mb->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice->listX[LIST_0 + list_offset];
    storable_picture** list1 = slice->listX[LIST_1 + list_offset];

    int block_y_aff;
    if (slice->MbaffFrameFlag && mb->mb_field_decoding_flag)
        block_y_aff = (mb->mb.y * 4 - 4 * (mb->mbAddrX % 2)) / 2;
    else
        block_y_aff = mb->mb.y * 4;

    for (int block4x4 = 0; block4x4 < 16; ++block4x4) {
        if (mb->SubMbType[block4x4 / 4] != 0)
            continue;
        if (dir)
            mb->SubMbPredMode[block4x4 / 4] = 2;

        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
        int i4 = mb->mb.x * 4 + i;
        int j4 = mb->mb.y * 4 + j;
        int j6 = block_y_aff + j;
        auto mv_info = &dec_picture->mv_info[j4][i4];
        pic_motion_params *colocated;
        if (sps->direct_8x8_inference_flag)
            colocated = &list1[0]->mv_info[RSD(j6)][RSD(i4)];
        else
            colocated = &list1[0]->mv_info[j6][i4];
        if (sps->separate_colour_plane_flag && sps->chroma_format_idc == YUV444)
            colocated = &list1[0]->JVmv_info[slice->colour_plane_id][RSD(j6)][RSD(i4)];

        if (sps->direct_8x8_inference_flag) {
            if (slice->MbaffFrameFlag) {
                if (!mb->mb_field_decoding_flag &&
                    ((slice->listX[LIST_1][0]->iCodingType == FRAME_MB_PAIR_CODING &&
                      slice->listX[LIST_1][0]->motion.mb_field_decoding_flag[mb->mbAddrX]) ||
                     (slice->listX[LIST_1][0]->iCodingType == FIELD_CODING))) {
                    if (abs(dec_picture->poc - slice->listX[LIST_1+4][0]->poc) >
                        abs(dec_picture->poc - slice->listX[LIST_1+2][0]->poc))
                        colocated = &slice->listX[LIST_1+2][0]->mv_info[RSD(j6)>>1][RSD(i4)];
                    else
                        colocated = &slice->listX[LIST_1+4][0]->mv_info[RSD(j6)>>1][RSD(i4)];
                }
            } else if (!sps->frame_mbs_only_flag && !slice->field_pic_flag &&
                       slice->listX[LIST_1][0]->iCodingType != FRAME_CODING) {
                if (abs(dec_picture->poc - list1[0]->bottom_field->poc) >
                    abs(dec_picture->poc - list1[0]->top_field->poc) )
                    colocated = &list1[0]->top_field->mv_info[RSD(j6)>>1][RSD(i4)];
                else
                    colocated = &list1[0]->bottom_field->mv_info[RSD(j6)>>1][RSD(i4)];
            } else if (!sps->frame_mbs_only_flag && slice->field_pic_flag &&
                       slice->structure != list1[0]->structure && list1[0]->coded_frame) {
                if (!slice->bottom_field_flag)
                    colocated = &list1[0]->frame->top_field->mv_info[RSD(j6)][RSD(i4)];
                else
                    colocated = &list1[0]->frame->bottom_field->mv_info[RSD(j6)][RSD(i4)];
            }
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

            if (sps->direct_8x8_inference_flag &&
                ((slice->MbaffFrameFlag &&
                  ((mb->mb_field_decoding_flag && colocated->ref_pic[refList]->structure == FRAME) || 
                   (!mb->mb_field_decoding_flag && colocated->ref_pic[refList]->structure != FRAME))) ||
                 (!slice->MbaffFrameFlag &&
                  ((slice->field_pic_flag == 0 && colocated->ref_pic[refList]->structure != FRAME) ||
                   (slice->field_pic_flag == 1 && colocated->ref_pic[refList]->structure == FRAME))))) {
                for (int iref = 0; iref < min<int>(slice->num_ref_idx_l0_active_minus1+1, slice->listXsize[LIST_0 + list_offset]);iref++) {
                    if (slice->listX[LIST_0 + list_offset][iref]->top_field == colocated->ref_pic[refList] || 
                        slice->listX[LIST_0 + list_offset][iref]->bottom_field == colocated->ref_pic[refList] ||
                        slice->listX[LIST_0 + list_offset][iref]->frame == colocated->ref_pic[refList]) {
                        if (slice->field_pic_flag == 1 &&
                            slice->listX[LIST_0 + list_offset][iref]->structure != slice->structure)
                            mapped_idx = INVALIDINDEX;
                        else {
                            mapped_idx = iref;
                            break;
                        }
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx = INVALIDINDEX;
                }
            } else {
                for (int iref = 0; iref < min<int>(slice->num_ref_idx_l0_active_minus1+1, slice->listXsize[LIST_0 + list_offset]);iref++) {
                    if (slice->listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                        mapped_idx = iref;            
                        break;
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx = INVALIDINDEX;
                }
            }

            if (INVALIDINDEX == mapped_idx)
                error("temporal direct error: colocated block has ref that is unavailable",-1111);

            int mv_scale = slice->mvscale[LIST_0 + list_offset][mapped_idx];
            int mv_x = colocated->mv[refList].mv_x;
            int mv_y = colocated->mv[refList].mv_y;

            if (sps->direct_8x8_inference_flag) {
                if ((slice->MbaffFrameFlag && !mb->mb_field_decoding_flag &&
                     colocated->ref_pic[refList]->structure!=FRAME) ||
                    (!slice->MbaffFrameFlag && slice->field_pic_flag == 0 &&
                     colocated->ref_pic[refList]->structure != FRAME))
                    mv_y *= 2;
                else if ((slice->MbaffFrameFlag && mb->mb_field_decoding_flag &&
                          colocated->ref_pic[refList]->structure==FRAME) ||
                         (!slice->MbaffFrameFlag && slice->field_pic_flag == 1 &&
                          colocated->ref_pic[refList]->structure==FRAME))
                    mv_y /= 2;
            }

            //! In such case, an array is needed for each different reference.
            if (mv_scale == 9999 || slice->listX[LIST_0 + list_offset][mapped_idx]->is_long_term) {
                mv_info->mv[LIST_0].mv_x = (short) mv_x;
                mv_info->mv[LIST_0].mv_y = (short) mv_y;
                mv_info->mv[LIST_1] = zero_mv;
            } else {
                mv_info->mv[LIST_0].mv_x = (short) ((mv_scale * mv_x + 128 ) >> 8);
                mv_info->mv[LIST_0].mv_y = (short) ((mv_scale * mv_y + 128 ) >> 8);

                mv_info->mv[LIST_1].mv_x = (short) (mv_info->mv[LIST_0].mv_x - mv_x);
                mv_info->mv[LIST_1].mv_y = (short) (mv_info->mv[LIST_0].mv_y - mv_y);
            }

            mv_info->ref_idx[LIST_0] = (char) mapped_idx;
            mv_info->ref_idx[LIST_1] = 0;
        }
        // store reference picture ID determined by direct mode
        mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
        mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
    }
}

int inter_prediction_t::get_colocated_info(mb_t* mb, storable_picture* list1, int i, int j)
{
    if (list1->is_long_term)
        return 1;

    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;

    pic_motion_params* fs;
    if (sps->direct_8x8_inference_flag)
        fs = &list1->mv_info[RSD(j)][RSD(i)];
    else
        fs = &list1->mv_info[j][i];

    if (sps->direct_8x8_inference_flag) {
        if (slice->MbaffFrameFlag ||
            (!sps->frame_mbs_only_flag &&
             ((!slice->field_pic_flag && list1->iCodingType == FIELD_CODING) ||
              (slice->structure != list1->structure && list1->coded_frame)))) {
            if (slice->field_pic_flag && slice->structure != list1->structure && list1->coded_frame) {
                if (!slice->bottom_field_flag)
                    fs = &list1->top_field->mv_info[RSD(j)][RSD(i)];
                else
                    fs = &list1->bottom_field->mv_info[RSD(j)][RSD(i)];
            } else if ((!mb->mb_field_decoding_flag &&
                        (list1->iCodingType == FIELD_CODING || list1->motion.mb_field_decoding_flag[mb->mbAddrX]))) {
                if (abs(slice->dec_picture->poc - list1->bottom_field->poc) >
                    abs(slice->dec_picture->poc - list1->top_field->poc))
                    fs = &list1->top_field->mv_info[RSD(j) >> 1][RSD(i)];
                else
                    fs = &list1->bottom_field->mv_info[RSD(j) >> 1][RSD(i)];
            }
        } else if (sps->separate_colour_plane_flag && sps->chroma_format_idc == YUV444)
            fs = &list1->JVmv_info[slice->colour_plane_id][RSD(j)][RSD(i)];
    }

    int moving =
        !((fs->ref_idx[LIST_0] == 0 &&
           abs(fs->mv[LIST_0].mv_x) >> 1 == 0 && abs(fs->mv[LIST_0].mv_y) >> 1 == 0) ||
          (fs->ref_idx[LIST_0] == -1 && fs->ref_idx[LIST_1] == 0 &&
           abs(fs->mv[LIST_1].mv_x) >> 1 == 0 && abs(fs->mv[LIST_1].mv_y) >> 1 == 0));
    return moving;  
}

inline void inter_prediction_t::set_direct_references(mb_t* mb, const PixelPos* pix, char* l0_rFrame, char* l1_rFrame, pic_motion_params** mv_info)
{
    slice_t* slice = mb->p_Slice;

    if (!pix->available) {
        *l0_rFrame = -1;
        *l1_rFrame = -1;
        return;
    }

    char* ref_idx = mv_info[pix->pos_y][pix->pos_x].ref_idx;
    if (!slice->MbaffFrameFlag ||
        (mb->mb_field_decoding_flag == slice->mb_data[pix->mb_addr].mb_field_decoding_flag)) {
        *l0_rFrame = ref_idx[LIST_0];
        *l1_rFrame = ref_idx[LIST_1];
    } else if (mb->mb_field_decoding_flag) {
        *l0_rFrame = (ref_idx[LIST_0] < 0) ? ref_idx[LIST_0] : ref_idx[LIST_0] * 2;
        *l1_rFrame = (ref_idx[LIST_1] < 0) ? ref_idx[LIST_1] : ref_idx[LIST_1] * 2;
    } else {
        *l0_rFrame = (ref_idx[LIST_0] >> 1);
        *l1_rFrame = (ref_idx[LIST_1] >> 1);
    }
}

void inter_prediction_t::prepare_direct_params(mb_t* mb, MotionVector* pmvl0, MotionVector* pmvl1, char* l0_rFrame, char* l1_rFrame)
{
    slice_t* slice = mb->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;

    PixelPos pix[4];
    get_neighbors(mb, pix, 0, 0, 16);

    char l0_refA, l0_refB, l0_refC;
    char l1_refA, l1_refB, l1_refC;
    auto mv_info = dec_picture->mv_info;
    this->set_direct_references(mb, &pix[0], &l0_refA, &l1_refA, mv_info);
    this->set_direct_references(mb, &pix[1], &l0_refB, &l1_refB, mv_info);
    this->set_direct_references(mb, &pix[2], &l0_refC, &l1_refC, mv_info);

    *l0_rFrame = (char) min(min((unsigned char) l0_refA, (unsigned char) l0_refB), (unsigned char) l0_refC);
    *l1_rFrame = (char) min(min((unsigned char) l1_refA, (unsigned char) l1_refB), (unsigned char) l1_refC);

    if (*l0_rFrame >=0)
        GetMVPredictor(mb, pix, pmvl0, *l0_rFrame, mv_info, LIST_0, 0, 0, 16, 16);
    if (*l1_rFrame >=0)
        GetMVPredictor(mb, pix, pmvl1, *l1_rFrame, mv_info, LIST_1, 0, 0, 16, 16);
}

void inter_prediction_t::get_direct_spatial(mb_t* mb, bool dir)
{
    bool has_direct = (mb->SubMbType[0] == 0) | (mb->SubMbType[1] == 0) |
                      (mb->SubMbType[2] == 0) | (mb->SubMbType[3] == 0);
    if (!has_direct)
        return;

    slice_t* slice = mb->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;

    int list_offset = slice->MbaffFrameFlag && mb->mb_field_decoding_flag ?
                      mb->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice->listX[LIST_0 + list_offset];
    storable_picture** list1 = slice->listX[LIST_1 + list_offset];

    char l0_rFrame, l1_rFrame;
    MotionVector pmvl0, pmvl1;
    prepare_direct_params(mb, &pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

    int pred_dir = 0;
    if (l0_rFrame < 0 && l1_rFrame < 0)
        pred_dir = 2;
    else
        pred_dir = l1_rFrame == -1 ? 0 : l0_rFrame == -1 ? 1 : 2;

    storable_picture* ref_pic_l[2];
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

    int blks = sps->direct_8x8_inference_flag ? 4 : 1;
    for (int block4x4 = 0; block4x4 < 16; block4x4 += blks) {
        if (mb->SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
            if (dir)
                mb->SubMbPredMode[block4x4 / 4] = pred_dir;

            auto mv_info = &dec_picture->mv_info[mb->mb.y * 4 + j][mb->mb.x * 4 + i];
            mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
            mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
            mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
            mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

            int block_y_aff;
            if (slice->MbaffFrameFlag && mb->mb_field_decoding_flag)
                block_y_aff = (mb->mb.y * 4 - 4 * (mb->mbAddrX % 2)) / 2;
            else
                block_y_aff = mb->mb.y * 4;
            bool is_not_moving = (get_colocated_info(mb, list1[0], mb->mb.x * 4 + i, block_y_aff + j) == 0);
            mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? zero_mv : pmvl0;
            mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? zero_mv : pmvl1;

            if (sps->direct_8x8_inference_flag) {
                dec_picture->mv_info[mb->mb.y * 4 + j + 0][mb->mb.x * 4 + i + 1] = *mv_info;
                dec_picture->mv_info[mb->mb.y * 4 + j + 1][mb->mb.x * 4 + i + 0] = *mv_info;
                dec_picture->mv_info[mb->mb.y * 4 + j + 1][mb->mb.x * 4 + i + 1] = *mv_info;
            }
        }
    }
}

int inter_prediction_t::get_inter8x8(mb_t* mb, int block8x8)
{
    slice_t* slice = mb->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;

    int list_offset = slice->MbaffFrameFlag && mb->mb_field_decoding_flag ?
                      mb->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice->listX[LIST_0 + list_offset];
    storable_picture** list1 = slice->listX[LIST_1 + list_offset];

    int pred_dir = mb->SubMbPredMode[block8x8];

    for (int block4x4 = block8x8 * 4; block4x4 < block8x8 * 4 + 4; ++block4x4) {
        if (mb->SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);

            auto mv_info = &dec_picture->mv_info[mb->mb.y * 4 + j][mb->mb.x * 4 + i];
            if (slice->direct_spatial_mv_pred_flag) {
                if (mv_info->ref_idx[LIST_1] == -1)
                    pred_dir = 0;
                else if (mv_info->ref_idx[LIST_0] == -1)
                    pred_dir = 1;
                else
                    pred_dir = 2;
                //mb->SubMbPredMode[block8x8] = pred_dir;
            } else {
                mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
                mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
            }
        }
    }

    return pred_dir;
}
