#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "neighbour.h"


#define INVALIDINDEX  (-135792468)


static inline int RSD(int x)
{
    return (x & 2) ? (x | 1) : (x & ~1);
}


enum {
    MVPRED_MEDIAN = 0,
    MVPRED_L      = 1,
    MVPRED_U      = 2,
    MVPRED_UR     = 3
};


static inline int imedian(int a, int b, int c)
{
    if (a > b) { // a > b
        if (b > c)
            return b; // a > b > c
        else if (a > c)
            return c; // a > c > b
        else
            return a; // c > a > b
    } else { // b > a
        if (a > c)
            return a; // b > a > c
        else if (b > c)
            return c; // b > c > a
        else
            return b; // c > b > a
    }
}

void macroblock_t::GetMotionVectorPredictorMBAFF(
    PixelPos* block, MotionVector* pmv,
    short ref_frame, pic_motion_params** mv_info,
    int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y)
{
    slice_t* slice = this->p_Slice;
    int mv_a, mv_b, mv_c, pred_vec = 0;
    int mvPredType, rFrameL, rFrameU, rFrameUR;
    int hv;

    mvPredType = MVPRED_MEDIAN;

    if (this->mb_field_decoding_flag) {
        rFrameL  = block[0].available
            ? (slice->mb_data[block[0].mb_addr].mb_field_decoding_flag
            ? mv_info[block[0].pos_y][block[0].pos_x].ref_idx[list]
            : mv_info[block[0].pos_y][block[0].pos_x].ref_idx[list] * 2) : -1;
        rFrameU  = block[1].available
            ? (slice->mb_data[block[1].mb_addr].mb_field_decoding_flag
            ? mv_info[block[1].pos_y][block[1].pos_x].ref_idx[list]
            : mv_info[block[1].pos_y][block[1].pos_x].ref_idx[list] * 2) : -1;
        rFrameUR = block[2].available
            ? (slice->mb_data[block[2].mb_addr].mb_field_decoding_flag
            ? mv_info[block[2].pos_y][block[2].pos_x].ref_idx[list]
            : mv_info[block[2].pos_y][block[2].pos_x].ref_idx[list] * 2) : -1;
    } else {
        rFrameL = block[0].available
            ? (slice->mb_data[block[0].mb_addr].mb_field_decoding_flag
            ? mv_info[block[0].pos_y][block[0].pos_x].ref_idx[list] >>1
            : mv_info[block[0].pos_y][block[0].pos_x].ref_idx[list]) : -1;
        rFrameU  = block[1].available
            ? (slice->mb_data[block[1].mb_addr].mb_field_decoding_flag
            ? mv_info[block[1].pos_y][block[1].pos_x].ref_idx[list] >>1
            : mv_info[block[1].pos_y][block[1].pos_x].ref_idx[list]) : -1;
        rFrameUR = block[2].available
            ? (slice->mb_data[block[2].mb_addr].mb_field_decoding_flag
            ? mv_info[block[2].pos_y][block[2].pos_x].ref_idx[list] >>1
            : mv_info[block[2].pos_y][block[2].pos_x].ref_idx[list]) : -1;
    }

    /* Prediction if only one of the neighbors uses the reference frame
     *  we are checking
     */
    if (rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       
        mvPredType = MVPRED_L;
    else if (rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  
        mvPredType = MVPRED_U;
    else if (rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  
        mvPredType = MVPRED_UR;
    // Directional predictions
    if (blockshape_x == 8 && blockshape_y == 16) {
        if (mb_x == 0) {
            if (rFrameL == ref_frame)
                mvPredType = MVPRED_L;
        } else {
            if (rFrameUR == ref_frame)
                mvPredType = MVPRED_UR;
        }
    } else if (blockshape_x == 16 && blockshape_y == 8) {
        if (mb_y == 0) {
            if (rFrameU == ref_frame)
                mvPredType = MVPRED_U;
        } else {
            if (rFrameL == ref_frame)
                mvPredType = MVPRED_L;
        }
    }

    for (hv = 0; hv < 2; hv++) {
        if (hv == 0) {
            mv_a = block[0].available ? mv_info[block[0].pos_y][block[0].pos_x].mv[list].mv_x : 0;
            mv_b = block[1].available ? mv_info[block[1].pos_y][block[1].pos_x].mv[list].mv_x : 0;
            mv_c = block[2].available ? mv_info[block[2].pos_y][block[2].pos_x].mv[list].mv_x : 0;
        } else {
            if (this->mb_field_decoding_flag) {
                mv_a = block[0].available
                    ? slice->mb_data[block[0].mb_addr].mb_field_decoding_flag
                    ? mv_info[block[0].pos_y][block[0].pos_x].mv[list].mv_y
                    : mv_info[block[0].pos_y][block[0].pos_x].mv[list].mv_y / 2
                    : 0;
                mv_b = block[1].available
                    ? slice->mb_data[block[1].mb_addr].mb_field_decoding_flag
                    ? mv_info[block[1].pos_y][block[1].pos_x].mv[list].mv_y
                    : mv_info[block[1].pos_y][block[1].pos_x].mv[list].mv_y / 2
                    : 0;
                mv_c = block[2].available
                    ? slice->mb_data[block[2].mb_addr].mb_field_decoding_flag
                    ? mv_info[block[2].pos_y][block[2].pos_x].mv[list].mv_y
                    : mv_info[block[2].pos_y][block[2].pos_x].mv[list].mv_y / 2
                    : 0;
            } else {
                mv_a = block[0].available
                    ? slice->mb_data[block[0].mb_addr].mb_field_decoding_flag
                    ? mv_info[block[0].pos_y][block[0].pos_x].mv[list].mv_y * 2
                    : mv_info[block[0].pos_y][block[0].pos_x].mv[list].mv_y
                    : 0;
                mv_b = block[1].available
                    ? slice->mb_data[block[1].mb_addr].mb_field_decoding_flag
                    ? mv_info[block[1].pos_y][block[1].pos_x].mv[list].mv_y * 2
                    : mv_info[block[1].pos_y][block[1].pos_x].mv[list].mv_y
                    : 0;
                mv_c = block[2].available
                    ? slice->mb_data[block[2].mb_addr].mb_field_decoding_flag
                    ? mv_info[block[2].pos_y][block[2].pos_x].mv[list].mv_y * 2
                    : mv_info[block[2].pos_y][block[2].pos_x].mv[list].mv_y
                    : 0;
            }
        }

        switch (mvPredType) {
        case MVPRED_MEDIAN:
            if (!(block[1].available || block[2].available))
                pred_vec = mv_a;
            else
                pred_vec = imedian(mv_a, mv_b, mv_c);
            break;
        case MVPRED_L:
            pred_vec = mv_a;
            break;
        case MVPRED_U:
            pred_vec = mv_b;
            break;
        case MVPRED_UR:
            pred_vec = mv_c;
            break;
        default:
            break;
        }

        if (hv == 0)
            pmv->mv_x = (short) pred_vec;
        else
            pmv->mv_y = (short) pred_vec;
    }
}

void macroblock_t::GetMotionVectorPredictorNormal(
    PixelPos* block, MotionVector* pmv,
    short ref_frame, pic_motion_params** mv_info,
    int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y)
{
    int mvPredType = MVPRED_MEDIAN;

    int rFrameL    = block[0].available ? mv_info[block[0].pos_y][block[0].pos_x].ref_idx[list] : -1;
    int rFrameU    = block[1].available ? mv_info[block[1].pos_y][block[1].pos_x].ref_idx[list] : -1;
    int rFrameUR   = block[2].available ? mv_info[block[2].pos_y][block[2].pos_x].ref_idx[list] : -1;

    /* Prediction if only one of the neighbors uses the reference frame
     *  we are checking
     */
    if (rFrameL == ref_frame && rFrameU != ref_frame && rFrameUR != ref_frame)       
        mvPredType = MVPRED_L;
    else if(rFrameL != ref_frame && rFrameU == ref_frame && rFrameUR != ref_frame)  
        mvPredType = MVPRED_U;
    else if(rFrameL != ref_frame && rFrameU != ref_frame && rFrameUR == ref_frame)  
        mvPredType = MVPRED_UR;

    // Directional predictions
    if (blockshape_x == 8 && blockshape_y == 16) {
        if (mb_x == 0) {
            if (rFrameL == ref_frame)
                mvPredType = MVPRED_L;
        } else {
            if (rFrameUR == ref_frame)
                mvPredType = MVPRED_UR;
        }
    } else if (blockshape_x == 16 && blockshape_y == 8) {
        if (mb_y == 0) {
            if (rFrameU == ref_frame)
                mvPredType = MVPRED_U;
        } else {
            if (rFrameL == ref_frame)
                mvPredType = MVPRED_L;
        }
    }

    switch (mvPredType) {
    case MVPRED_MEDIAN:
        if (!(block[1].available || block[2].available)) {
            if (block[0].available)
                *pmv = mv_info[block[0].pos_y][block[0].pos_x].mv[list];
            else
                *pmv = zero_mv;
        } else {
            MotionVector *mv_a = block[0].available ? &mv_info[block[0].pos_y][block[0].pos_x].mv[list] : (MotionVector *) &zero_mv;
            MotionVector *mv_b = block[1].available ? &mv_info[block[1].pos_y][block[1].pos_x].mv[list] : (MotionVector *) &zero_mv;
            MotionVector *mv_c = block[2].available ? &mv_info[block[2].pos_y][block[2].pos_x].mv[list] : (MotionVector *) &zero_mv;

            pmv->mv_x = (short) imedian(mv_a->mv_x, mv_b->mv_x, mv_c->mv_x);
            pmv->mv_y = (short) imedian(mv_a->mv_y, mv_b->mv_y, mv_c->mv_y);
        }
        break;
    case MVPRED_L:
        if (block[0].available)
            *pmv = mv_info[block[0].pos_y][block[0].pos_x].mv[list];
        else
            *pmv = zero_mv;
        break;
    case MVPRED_U:
        if (block[1].available)
            *pmv = mv_info[block[1].pos_y][block[1].pos_x].mv[list];
        else
            *pmv = zero_mv;
        break;
    case MVPRED_UR:
        if (block[2].available)
            *pmv = mv_info[block[2].pos_y][block[2].pos_x].mv[list];
        else
            *pmv = zero_mv;
        break;
    default:
        break;
    }
}

void macroblock_t::GetMVPredictor(PixelPos* block, MotionVector* pmv,
                                  short ref_frame, pic_motion_params** mv_info,
                                  int list, int mb_x, int mb_y, int blockshape_x, int blockshape_y)
{
    slice_t* slice = this->p_Slice;

    if (slice->MbaffFrameFlag)
        this->GetMotionVectorPredictorMBAFF(block, pmv, ref_frame, mv_info,
                                            list, mb_x, mb_y, blockshape_x, blockshape_y);
    else
        this->GetMotionVectorPredictorNormal(block, pmv, ref_frame, mv_info,
                                             list, mb_x, mb_y, blockshape_x, blockshape_y);
}


void macroblock_t::get_direct_temporal(bool dir)
{
    bool has_direct = (this->SubMbType[0] == 0) | (this->SubMbType[1] == 0) |
                      (this->SubMbType[2] == 0) | (this->SubMbType[3] == 0);
    if (!has_direct)
        return;

    slice_t* slice = this->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;
    sps_t* sps = slice->active_sps;

    int list_offset = slice->MbaffFrameFlag && this->mb_field_decoding_flag ?
                      this->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice->listX[LIST_0 + list_offset];
    storable_picture** list1 = slice->listX[LIST_1 + list_offset];

    int block_y_aff;
    if (slice->MbaffFrameFlag && this->mb_field_decoding_flag)
        block_y_aff = (this->mb.y * 4 - 4 * (this->mbAddrX % 2)) / 2;
    else
        block_y_aff = this->mb.y * 4;

    for (int block4x4 = 0; block4x4 < 16; ++block4x4) {
        if (this->SubMbType[block4x4 / 4] != 0)
            continue;
        if (dir)
            this->SubMbPredMode[block4x4 / 4] = 2;

        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
        int i4 = this->mb.x * 4 + i;
        int j4 = this->mb.y * 4 + j;
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
                if (!this->mb_field_decoding_flag &&
                    ((slice->listX[LIST_1][0]->iCodingType == FRAME_MB_PAIR_CODING &&
                      slice->listX[LIST_1][0]->motion.mb_field_decoding_flag[this->mbAddrX]) ||
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
                  ((this->mb_field_decoding_flag && colocated->ref_pic[refList]->structure == FRAME) || 
                   (!this->mb_field_decoding_flag && colocated->ref_pic[refList]->structure != FRAME))) ||
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
                if ((slice->MbaffFrameFlag && !this->mb_field_decoding_flag &&
                     colocated->ref_pic[refList]->structure!=FRAME) ||
                    (!slice->MbaffFrameFlag && slice->field_pic_flag == 0 &&
                     colocated->ref_pic[refList]->structure != FRAME))
                    mv_y *= 2;
                else if ((slice->MbaffFrameFlag && this->mb_field_decoding_flag &&
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

int macroblock_t::get_colocated_info(storable_picture* list1, int i, int j)
{
    if (list1->is_long_term)
        return 1;

    slice_t* slice = this->p_Slice;
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
            } else if ((!this->mb_field_decoding_flag &&
                        (list1->iCodingType == FIELD_CODING || list1->motion.mb_field_decoding_flag[this->mbAddrX]))) {
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

inline void macroblock_t::set_direct_references(const PixelPos* pix, char* l0_rFrame, char* l1_rFrame, pic_motion_params** mv_info)
{
    slice_t* slice = this->p_Slice;

    if (!pix->available) {
        *l0_rFrame = -1;
        *l1_rFrame = -1;
        return;
    }

    char* ref_idx = mv_info[pix->pos_y][pix->pos_x].ref_idx;
    if (!slice->MbaffFrameFlag ||
        (this->mb_field_decoding_flag == slice->mb_data[pix->mb_addr].mb_field_decoding_flag)) {
        *l0_rFrame = ref_idx[LIST_0];
        *l1_rFrame = ref_idx[LIST_1];
    } else if (this->mb_field_decoding_flag) {
        *l0_rFrame = (ref_idx[LIST_0] < 0) ? ref_idx[LIST_0] : ref_idx[LIST_0] * 2;
        *l1_rFrame = (ref_idx[LIST_1] < 0) ? ref_idx[LIST_1] : ref_idx[LIST_1] * 2;
    } else {
        *l0_rFrame = (ref_idx[LIST_0] >> 1);
        *l1_rFrame = (ref_idx[LIST_1] >> 1);
    }
}

void macroblock_t::prepare_direct_params(MotionVector* pmvl0, MotionVector* pmvl1, char* l0_rFrame, char* l1_rFrame)
{
    slice_t* slice = this->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;

    PixelPos pix[4];
    get_neighbors(this, pix, 0, 0, 16);

    char l0_refA, l0_refB, l0_refC;
    char l1_refA, l1_refB, l1_refC;
    auto mv_info = dec_picture->mv_info;
    this->set_direct_references(&pix[0], &l0_refA, &l1_refA, mv_info);
    this->set_direct_references(&pix[1], &l0_refB, &l1_refB, mv_info);
    this->set_direct_references(&pix[2], &l0_refC, &l1_refC, mv_info);

    *l0_rFrame = (char) min(min((unsigned char) l0_refA, (unsigned char) l0_refB), (unsigned char) l0_refC);
    *l1_rFrame = (char) min(min((unsigned char) l1_refA, (unsigned char) l1_refB), (unsigned char) l1_refC);

    if (*l0_rFrame >=0)
        this->GetMVPredictor(pix, pmvl0, *l0_rFrame, mv_info, LIST_0, 0, 0, 16, 16);
    if (*l1_rFrame >=0)
        this->GetMVPredictor(pix, pmvl1, *l1_rFrame, mv_info, LIST_1, 0, 0, 16, 16);
}

void macroblock_t::get_direct_spatial(bool dir)
{
    bool has_direct = (this->SubMbType[0] == 0) | (this->SubMbType[1] == 0) |
                      (this->SubMbType[2] == 0) | (this->SubMbType[3] == 0);
    if (!has_direct)
        return;

    slice_t* slice = this->p_Slice;
    sps_t* sps = slice->active_sps;
    storable_picture* dec_picture = slice->dec_picture;

    int list_offset = slice->MbaffFrameFlag && this->mb_field_decoding_flag ?
                      this->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice->listX[LIST_0 + list_offset];
    storable_picture** list1 = slice->listX[LIST_1 + list_offset];

    char l0_rFrame, l1_rFrame;
    MotionVector pmvl0, pmvl1;
    this->prepare_direct_params(&pmvl0, &pmvl1, &l0_rFrame, &l1_rFrame);

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
        if (this->SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
            if (dir)
                this->SubMbPredMode[block4x4 / 4] = pred_dir;

            auto mv_info = &dec_picture->mv_info[this->mb.y * 4 + j][this->mb.x * 4 + i];
            mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
            mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
            mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
            mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

            int block_y_aff;
            if (slice->MbaffFrameFlag && this->mb_field_decoding_flag)
                block_y_aff = (this->mb.y * 4 - 4 * (this->mbAddrX % 2)) / 2;
            else
                block_y_aff = this->mb.y * 4;
            bool is_not_moving = (this->get_colocated_info(list1[0], this->mb.x * 4 + i, block_y_aff + j) == 0);
            mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? zero_mv : pmvl0;
            mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? zero_mv : pmvl1;

            if (sps->direct_8x8_inference_flag) {
                dec_picture->mv_info[this->mb.y * 4 + j + 0][this->mb.x * 4 + i + 1] = *mv_info;
                dec_picture->mv_info[this->mb.y * 4 + j + 1][this->mb.x * 4 + i + 0] = *mv_info;
                dec_picture->mv_info[this->mb.y * 4 + j + 1][this->mb.x * 4 + i + 1] = *mv_info;
            }
        }
    }
}

int macroblock_t::get_inter8x8(int block8x8)
{
    slice_t* slice = this->p_Slice;
    storable_picture* dec_picture = slice->dec_picture;

    int list_offset = slice->MbaffFrameFlag && this->mb_field_decoding_flag ?
                      this->mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice->listX[LIST_0 + list_offset];
    storable_picture** list1 = slice->listX[LIST_1 + list_offset];

    int pred_dir = this->SubMbPredMode[block8x8];

    for (int block4x4 = block8x8 * 4; block4x4 < block8x8 * 4 + 4; ++block4x4) {
        if (this->SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);

            auto mv_info = &dec_picture->mv_info[this->mb.y * 4 + j][this->mb.x * 4 + i];
            if (slice->direct_spatial_mv_pred_flag) {
                if (mv_info->ref_idx[LIST_1] == -1)
                    pred_dir = 0;
                else if (mv_info->ref_idx[LIST_0] == -1)
                    pred_dir = 1;
                else
                    pred_dir = 2;
                //this->SubMbPredMode[block8x8] = pred_dir;
            } else {
                mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
                mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
            }
        }
    }

    return pred_dir;
}
