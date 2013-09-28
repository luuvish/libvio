#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "dpb.h"
#include "neighbour.h"


namespace vio  {
namespace h264 {


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


mv_t Parser::Macroblock::GetMVPredictor(char ref_frame, int list, int i, int j, int step_h4, int step_v4)
{
    nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4 - 1      , j * 4    });
    nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4          , j * 4 - 1});
    nb_t nbC = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4 + step_h4, j * 4 - 1});
    nb_t nbD = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4 - 1      , j * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbC.mb = nbC.mb && nbC.mb->slice_nr == mb.slice_nr ? nbC.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    if (j > 0) {
        if (i < 2) { // first column of 8x8 blocks
            if (j == 2) {
                if (step_h4 == 16)
                    nbC.mb = nullptr;
            } else if (i * 4 + step_h4 == 8)
                nbC.mb = nullptr;
        } else if (i * 4 + step_h4 == 16)
            nbC.mb = nullptr;
    }

    if (!nbC.mb)
        nbC = nbD;

    char refIdxA = -1;
    mv_t mvA = {0, 0};
    if (nbA.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbA.y / 4][nbA.x / 4];
        refIdxA = mv_info->ref_idx[list];
        mvA     = mv_info->mv     [list];
        if (mb.mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag) {
            if (refIdxA >= 0)
                refIdxA *= 2;
            mvA.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag) {
            refIdxA >>= 1;
            mvA.mv_y *= 2;
        }
    }
    char refIdxB = -1;
    mv_t mvB = {0, 0};
    if (nbB.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbB.y / 4][nbB.x / 4];
        refIdxB = mv_info->ref_idx[list];
        mvB     = mv_info->mv     [list];
        if (mb.mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag) {
            if (refIdxB >= 0)
                refIdxB *= 2;
            mvB.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag) {
            refIdxB >>= 1;
            mvB.mv_y *= 2;
        }
    }
    char refIdxC = -1;
    mv_t mvC = {0, 0};
    if (nbC.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbC.y / 4][nbC.x / 4];
        refIdxC = mv_info->ref_idx[list];
        mvC     = mv_info->mv     [list];
        if (mb.mb_field_decoding_flag && !nbC.mb->mb_field_decoding_flag) {
            if (refIdxC >= 0)
                refIdxC *= 2;
            mvC.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbC.mb->mb_field_decoding_flag) {
            refIdxC >>= 1;
            mvC.mv_y *= 2;
        }
    }

    /* Prediction if only one of the neighbors uses the reference frame
     *  we are checking
     */
    int mvPredType = MVPRED_MEDIAN;
    if (refIdxA == ref_frame && refIdxB != ref_frame && refIdxC != ref_frame)       
        mvPredType = MVPRED_L;
    else if (refIdxA != ref_frame && refIdxB == ref_frame && refIdxC != ref_frame)  
        mvPredType = MVPRED_U;
    else if (refIdxA != ref_frame && refIdxB != ref_frame && refIdxC == ref_frame)  
        mvPredType = MVPRED_UR;

    // Directional predictions
    if (step_h4 == 8 && step_v4 == 16) {
        if (i == 0) {
            if (refIdxA == ref_frame)
                mvPredType = MVPRED_L;
        } else {
            if (refIdxC == ref_frame)
                mvPredType = MVPRED_UR;
        }
    } else if (step_h4 == 16 && step_v4 == 8) {
        if (j == 0) {
            if (refIdxB == ref_frame)
                mvPredType = MVPRED_U;
        } else {
            if (refIdxA == ref_frame)
                mvPredType = MVPRED_L;
        }
    }

    mv_t pred_mv;
    switch (mvPredType) {
    case MVPRED_MEDIAN:
        if (!(nbB.mb || nbC.mb))
            pred_mv = mvA;
        else {
            pred_mv.mv_x = median(mvA.mv_x, mvB.mv_x, mvC.mv_x);
            pred_mv.mv_y = median(mvA.mv_y, mvB.mv_y, mvC.mv_y);
        }
        break;
    case MVPRED_L:
        pred_mv = mvA;
        break;
    case MVPRED_U:
        pred_mv = mvB;
        break;
    case MVPRED_UR:
        pred_mv = mvC;
        break;
    default:
        break;
    }

    return pred_mv;
}

mv_t Parser::Macroblock::GetMVPredictor2(char& ref_frame, int list, int i, int j, int step_h4, int step_v4)
{
    nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4 - 1      , j * 4    });
    nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4          , j * 4 - 1});
    nb_t nbC = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4 + step_h4, j * 4 - 1});
    nb_t nbD = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {i * 4 - 1      , j * 4 - 1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbC.mb = nbC.mb && nbC.mb->slice_nr == mb.slice_nr ? nbC.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    if (j > 0) {
        if (i < 2) { // first column of 8x8 blocks
            if (j == 2) {
                if (step_h4 == 16)
                    nbC.mb = nullptr;
            } else if (i * 4 + step_h4 == 8)
                nbC.mb = nullptr;
        } else if (i * 4 + step_h4 == 16)
            nbC.mb = nullptr;
    }

    if (!nbC.mb)
        nbC = nbD;

    char refIdxA = -1;
    mv_t mvA = {0, 0};
    if (nbA.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbA.y / 4][nbA.x / 4];
        refIdxA = mv_info->ref_idx[list];
        mvA     = mv_info->mv     [list];
        if (mb.mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag) {
            if (refIdxA >= 0)
                refIdxA *= 2;
            mvA.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag) {
            refIdxA >>= 1;
            mvA.mv_y *= 2;
        }
    }
    char refIdxB = -1;
    mv_t mvB = {0, 0};
    if (nbB.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbB.y / 4][nbB.x / 4];
        refIdxB = mv_info->ref_idx[list];
        mvB     = mv_info->mv     [list];
        if (mb.mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag) {
            if (refIdxB >= 0)
                refIdxB *= 2;
            mvB.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag) {
            refIdxB >>= 1;
            mvB.mv_y *= 2;
        }
    }
    char refIdxC = -1;
    mv_t mvC = {0, 0};
    if (nbC.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbC.y / 4][nbC.x / 4];
        refIdxC = mv_info->ref_idx[list];
        mvC     = mv_info->mv     [list];
        if (mb.mb_field_decoding_flag && !nbC.mb->mb_field_decoding_flag) {
            if (refIdxC >= 0)
                refIdxC *= 2;
            mvC.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbC.mb->mb_field_decoding_flag) {
            refIdxC >>= 1;
            mvC.mv_y *= 2;
        }
    }

    ref_frame = (char) min(min((unsigned char) refIdxA, (unsigned char) refIdxB), (unsigned char) refIdxC);

    /* Prediction if only one of the neighbors uses the reference frame
     *  we are checking
     */
    int mvPredType = MVPRED_MEDIAN;
    if (refIdxA == ref_frame && refIdxB != ref_frame && refIdxC != ref_frame)       
        mvPredType = MVPRED_L;
    else if (refIdxA != ref_frame && refIdxB == ref_frame && refIdxC != ref_frame)  
        mvPredType = MVPRED_U;
    else if (refIdxA != ref_frame && refIdxB != ref_frame && refIdxC == ref_frame)  
        mvPredType = MVPRED_UR;

    // Directional predictions
    if (step_h4 == 8 && step_v4 == 16) {
        if (i == 0) {
            if (refIdxA == ref_frame)
                mvPredType = MVPRED_L;
        } else {
            if (refIdxC == ref_frame)
                mvPredType = MVPRED_UR;
        }
    } else if (step_h4 == 16 && step_v4 == 8) {
        if (j == 0) {
            if (refIdxB == ref_frame)
                mvPredType = MVPRED_U;
        } else {
            if (refIdxA == ref_frame)
                mvPredType = MVPRED_L;
        }
    }

    mv_t pred_mv;
    switch (mvPredType) {
    case MVPRED_MEDIAN:
        if (!(nbB.mb || nbC.mb))
            pred_mv = mvA;
        else {
            pred_mv.mv_x = median(mvA.mv_x, mvB.mv_x, mvC.mv_x);
            pred_mv.mv_y = median(mvA.mv_y, mvB.mv_y, mvC.mv_y);
        }
        break;
    case MVPRED_L:
        pred_mv = mvA;
        break;
    case MVPRED_U:
        pred_mv = mvB;
        break;
    case MVPRED_UR:
        pred_mv = mvC;
        break;
    default:
        break;
    }

    return pred_mv;
}

void Parser::Macroblock::skip_macroblock()
{
    int list_offset = slice.MbaffFrameFlag && mb.mb_field_decoding_flag ?
                      mb.mbAddrX % 2 ? 4 : 2 : 0;

    nb_t nbA = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {-1,  0});
    nb_t nbB = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, { 0, -1});
    nb_t nbC = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {16, -1});
    nb_t nbD = slice.neighbour.get_neighbour(&slice, false, mb.mbAddrX, {-1, -1});
    nbA.mb = nbA.mb && nbA.mb->slice_nr == mb.slice_nr ? nbA.mb : nullptr;
    nbB.mb = nbB.mb && nbB.mb->slice_nr == mb.slice_nr ? nbB.mb : nullptr;
    nbC.mb = nbC.mb && nbC.mb->slice_nr == mb.slice_nr ? nbC.mb : nullptr;
    nbD.mb = nbD.mb && nbD.mb->slice_nr == mb.slice_nr ? nbD.mb : nullptr;

    if (!nbC.mb)
        nbC = nbD;

    int refIdxA = 0;
    mv_t mvA = {0, 0};
    if (nbA.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbA.y / 4][nbA.x / 4];
        refIdxA = mv_info->ref_idx[LIST_0];
        mvA     = mv_info->mv     [LIST_0];
        if (mb.mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag) {
            refIdxA  *= 2;
            mvA.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag) {
            refIdxA >>= 1;
            mvA.mv_y *= 2;
        }
    }

    int refIdxB = 0;
    mv_t mvB = {0, 0};
    if (nbB.mb) {
        auto mv_info = &slice.dec_picture->mv_info[nbB.y / 4][nbB.x / 4];
        refIdxB = mv_info->ref_idx[LIST_0];
        mvB     = mv_info->mv     [LIST_0];
        if (mb.mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag) {
            refIdxB  *= 2;
            mvB.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag) {
            refIdxB >>= 1;
            mvB.mv_y *= 2;
        }
    }

    mb.CodedBlockPatternLuma   = 0;
    mb.CodedBlockPatternChroma = 0;
    if (!pps.entropy_coding_mode_flag)
        memset(mb.nz_coeff, 0, 3 * 16 * sizeof(uint8_t));

    bool zeroMotionA = !nbA.mb || (refIdxA == 0 && mvA.mv_x == 0 && mvA.mv_y == 0);
    bool zeroMotionB = !nbB.mb || (refIdxB == 0 && mvB.mv_x == 0 && mvB.mv_y == 0);
    mv_t pred_mv;
    if (zeroMotionA || zeroMotionB)
        pred_mv = {0, 0};
    else {
        int refIdxA = -1;
        mv_t mvA = {0, 0};
        if (nbA.mb) {
            auto mv_info = &slice.dec_picture->mv_info[nbA.y / 4][nbA.x / 4];
            refIdxA = mv_info->ref_idx[LIST_0];
            mvA     = mv_info->mv     [LIST_0];
            if (mb.mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag) {
                refIdxA  *= 2;
                mvA.mv_y /= 2;
            }
            if (!mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag) {
                refIdxA >>= 1;
                mvA.mv_y *= 2;
            }
        }
        int refIdxB = -1;
        mv_t mvB = {0, 0};
        if (nbB.mb) {
            auto mv_info = &slice.dec_picture->mv_info[nbB.y / 4][nbB.x / 4];
            refIdxB = mv_info->ref_idx[LIST_0];
            mvB     = mv_info->mv     [LIST_0];
            if (mb.mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag) {
                refIdxB  *= 2;
                mvB.mv_y /= 2;
            }
            if (!mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag) {
                refIdxB >>= 1;
                mvB.mv_y *= 2;
            }
        }
        int refIdxC = -1;
        mv_t mvC = {0, 0};
        if (nbC.mb) {
            auto mv_info = &slice.dec_picture->mv_info[nbC.y / 4][nbC.x / 4];
            refIdxC = mv_info->ref_idx[LIST_0];
            mvC     = mv_info->mv     [LIST_0];
            if (mb.mb_field_decoding_flag && !nbC.mb->mb_field_decoding_flag) {
                refIdxC *= 2;
                mvC.mv_y /= 2;
            }
            if (!mb.mb_field_decoding_flag && nbC.mb->mb_field_decoding_flag) {
                refIdxC >>= 1;
                mvC.mv_y *= 2;
            }
        }

        /* Prediction if only one of the neighbors uses the reference frame
         *  we are checking
         */
        short ref_frame = 0;
        int mvPredType = MVPRED_MEDIAN;
        if (refIdxA == ref_frame && refIdxB != ref_frame && refIdxC != ref_frame)       
            mvPredType = MVPRED_L;
        else if (refIdxA != ref_frame && refIdxB == ref_frame && refIdxC != ref_frame)  
            mvPredType = MVPRED_U;
        else if (refIdxA != ref_frame && refIdxB != ref_frame && refIdxC == ref_frame)  
            mvPredType = MVPRED_UR;

        switch (mvPredType) {
        case MVPRED_MEDIAN:
            if (!(nbB.mb || nbC.mb))
                pred_mv = mvA;
            else {
                pred_mv.mv_x = median(mvA.mv_x, mvB.mv_x, mvC.mv_x);
                pred_mv.mv_y = median(mvA.mv_y, mvB.mv_y, mvC.mv_y);
            }
            break;
        case MVPRED_L:
            pred_mv = mvA;
            break;
        case MVPRED_U:
            pred_mv = mvB;
            break;
        case MVPRED_UR:
            pred_mv = mvC;
            break;
        default:
            break;
        }
    }

    auto mv_info = slice.dec_picture->mv_info;
    storable_picture* cur_pic = slice.listX[list_offset][0];

    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            auto mv = &mv_info[mb.mb.y * 4 + y][mb.mb.x * 4 + x];
            mv->ref_pic[LIST_0] = cur_pic;
            mv->ref_idx[LIST_0] = 0;
            mv->mv     [LIST_0] = pred_mv;
        }
    }
}

void Parser::Macroblock::get_direct_temporal(bool dir)
{
    bool has_direct = (mb.SubMbType[0] == 0) | (mb.SubMbType[1] == 0) |
                      (mb.SubMbType[2] == 0) | (mb.SubMbType[3] == 0);
    if (!has_direct)
        return;

    int list_offset = slice.MbaffFrameFlag && mb.mb_field_decoding_flag ?
                      mb.mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice.listX[LIST_0 + list_offset];
    storable_picture** list1 = slice.listX[LIST_1 + list_offset];

    int block_y_aff;
    if (slice.MbaffFrameFlag && mb.mb_field_decoding_flag)
        block_y_aff = (mb.mb.y * 4 - 4 * (mb.mbAddrX % 2)) / 2;
    else
        block_y_aff = mb.mb.y * 4;

    for (int block4x4 = 0; block4x4 < 16; ++block4x4) {
        if (mb.SubMbType[block4x4 / 4] != 0)
            continue;
        if (dir)
            mb.SubMbPredMode[block4x4 / 4] = 2;

        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
        int i4 = mb.mb.x * 4 + i;
        int j4 = mb.mb.y * 4 + j;
        int j6 = block_y_aff + j;
        auto mv_info = &slice.dec_picture->mv_info[j4][i4];
        pic_motion_params *colocated;
        if (sps.direct_8x8_inference_flag)
            colocated = &list1[0]->mv_info[RSD(j6)][RSD(i4)];
        else
            colocated = &list1[0]->mv_info[j6][i4];
        if (sps.separate_colour_plane_flag && sps.chroma_format_idc == YUV444)
            colocated = &list1[0]->JVmv_info[slice.colour_plane_id][RSD(j6)][RSD(i4)];

        if (sps.direct_8x8_inference_flag) {
            if (slice.MbaffFrameFlag) {
                if (!mb.mb_field_decoding_flag &&
                    ((slice.listX[LIST_1][0]->iCodingType == FRAME_MB_PAIR_CODING &&
                      slice.listX[LIST_1][0]->motion.mb_field_decoding_flag[mb.mbAddrX]) ||
                     (slice.listX[LIST_1][0]->iCodingType == FIELD_CODING))) {
                    if (abs(slice.dec_picture->poc - slice.listX[LIST_1+4][0]->poc) >
                        abs(slice.dec_picture->poc - slice.listX[LIST_1+2][0]->poc))
                        colocated = &slice.listX[LIST_1+2][0]->mv_info[RSD(j6)>>1][RSD(i4)];
                    else
                        colocated = &slice.listX[LIST_1+4][0]->mv_info[RSD(j6)>>1][RSD(i4)];
                }
            } else if (!sps.frame_mbs_only_flag && !slice.field_pic_flag &&
                       slice.listX[LIST_1][0]->iCodingType != FRAME_CODING) {
                if (abs(slice.dec_picture->poc - list1[0]->bottom_field->poc) >
                    abs(slice.dec_picture->poc - list1[0]->top_field->poc) )
                    colocated = &list1[0]->top_field->mv_info[RSD(j6)>>1][RSD(i4)];
                else
                    colocated = &list1[0]->bottom_field->mv_info[RSD(j6)>>1][RSD(i4)];
            } else if (!sps.frame_mbs_only_flag && slice.field_pic_flag &&
                       slice.structure != list1[0]->structure && list1[0]->coded_frame) {
                if (!slice.bottom_field_flag)
                    colocated = &list1[0]->frame->top_field->mv_info[RSD(j6)][RSD(i4)];
                else
                    colocated = &list1[0]->frame->bottom_field->mv_info[RSD(j6)][RSD(i4)];
            }
        }

        int   refList = colocated->ref_idx[LIST_0] == -1 ? LIST_1 : LIST_0;
        short ref_idx = colocated->ref_idx[refList];

        if (ref_idx == -1) { // co-located is intra mode
            mv_info->mv[LIST_0] = {0, 0};
            mv_info->mv[LIST_1] = {0, 0};

            mv_info->ref_idx[LIST_0] = 0;
            mv_info->ref_idx[LIST_1] = 0;
        } else { // co-located skip or inter mode
            int mapped_idx = 0;

            if (sps.direct_8x8_inference_flag &&
                ((slice.MbaffFrameFlag &&
                  ((mb.mb_field_decoding_flag && colocated->ref_pic[refList]->structure == FRAME) || 
                   (!mb.mb_field_decoding_flag && colocated->ref_pic[refList]->structure != FRAME))) ||
                 (!slice.MbaffFrameFlag &&
                  ((slice.field_pic_flag == 0 && colocated->ref_pic[refList]->structure != FRAME) ||
                   (slice.field_pic_flag == 1 && colocated->ref_pic[refList]->structure == FRAME))))) {
                for (int iref = 0; iref < min<int>(slice.num_ref_idx_l0_active_minus1+1, slice.listXsize[LIST_0 + list_offset]);iref++) {
                    if (slice.listX[LIST_0 + list_offset][iref]->top_field == colocated->ref_pic[refList] || 
                        slice.listX[LIST_0 + list_offset][iref]->bottom_field == colocated->ref_pic[refList] ||
                        slice.listX[LIST_0 + list_offset][iref]->frame == colocated->ref_pic[refList]) {
                        if (slice.field_pic_flag == 1 &&
                            slice.listX[LIST_0 + list_offset][iref]->structure != slice.structure)
                            mapped_idx = INVALIDINDEX;
                        else {
                            mapped_idx = iref;
                            break;
                        }
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx = INVALIDINDEX;
                }
            } else {
                for (int iref = 0; iref < min<int>(slice.num_ref_idx_l0_active_minus1+1, slice.listXsize[LIST_0 + list_offset]);iref++) {
                    if (slice.listX[LIST_0 + list_offset][iref] == colocated->ref_pic[refList]) {
                        mapped_idx = iref;            
                        break;
                    } else //! invalid index. Default to zero even though this case should not happen
                        mapped_idx = INVALIDINDEX;
                }
            }

            if (INVALIDINDEX == mapped_idx)
                error("temporal direct error: colocated block has ref that is unavailable",-1111);

            int mv_scale = slice.mvscale[LIST_0 + list_offset][mapped_idx];
            int mv_x = colocated->mv[refList].mv_x;
            int mv_y = colocated->mv[refList].mv_y;

            if (sps.direct_8x8_inference_flag) {
                if ((slice.MbaffFrameFlag && !mb.mb_field_decoding_flag &&
                     colocated->ref_pic[refList]->structure!=FRAME) ||
                    (!slice.MbaffFrameFlag && slice.field_pic_flag == 0 &&
                     colocated->ref_pic[refList]->structure != FRAME))
                    mv_y *= 2;
                else if ((slice.MbaffFrameFlag && mb.mb_field_decoding_flag &&
                          colocated->ref_pic[refList]->structure==FRAME) ||
                         (!slice.MbaffFrameFlag && slice.field_pic_flag == 1 &&
                          colocated->ref_pic[refList]->structure==FRAME))
                    mv_y /= 2;
            }

            //! In such case, an array is needed for each different reference.
            if (mv_scale == 9999 || slice.listX[LIST_0 + list_offset][mapped_idx]->is_long_term) {
                mv_info->mv[LIST_0].mv_x = (short) mv_x;
                mv_info->mv[LIST_0].mv_y = (short) mv_y;
                mv_info->mv[LIST_1] = {0, 0};
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

int Parser::Macroblock::get_colocated_info(storable_picture* list1, int i, int j)
{
    if (list1->is_long_term)
        return 1;

    pic_motion_params* fs;
    if (sps.direct_8x8_inference_flag)
        fs = &list1->mv_info[RSD(j)][RSD(i)];
    else
        fs = &list1->mv_info[j][i];

    if (sps.direct_8x8_inference_flag) {
        if (slice.MbaffFrameFlag ||
            (!sps.frame_mbs_only_flag &&
             ((!slice.field_pic_flag && list1->iCodingType == FIELD_CODING) ||
              (slice.structure != list1->structure && list1->coded_frame)))) {
            if (slice.field_pic_flag && slice.structure != list1->structure && list1->coded_frame) {
                if (!slice.bottom_field_flag)
                    fs = &list1->top_field->mv_info[RSD(j)][RSD(i)];
                else
                    fs = &list1->bottom_field->mv_info[RSD(j)][RSD(i)];
            } else if ((!mb.mb_field_decoding_flag &&
                        (list1->iCodingType == FIELD_CODING || list1->motion.mb_field_decoding_flag[mb.mbAddrX]))) {
                if (abs(slice.dec_picture->poc - list1->bottom_field->poc) >
                    abs(slice.dec_picture->poc - list1->top_field->poc))
                    fs = &list1->top_field->mv_info[RSD(j) >> 1][RSD(i)];
                else
                    fs = &list1->bottom_field->mv_info[RSD(j) >> 1][RSD(i)];
            }
        } else if (sps.separate_colour_plane_flag && sps.chroma_format_idc == YUV444)
            fs = &list1->JVmv_info[slice.colour_plane_id][RSD(j)][RSD(i)];
    }

    int moving =
        !((fs->ref_idx[LIST_0] == 0 &&
           abs(fs->mv[LIST_0].mv_x) >> 1 == 0 && abs(fs->mv[LIST_0].mv_y) >> 1 == 0) ||
          (fs->ref_idx[LIST_0] == -1 && fs->ref_idx[LIST_1] == 0 &&
           abs(fs->mv[LIST_1].mv_x) >> 1 == 0 && abs(fs->mv[LIST_1].mv_y) >> 1 == 0));
    return moving;  
}

void Parser::Macroblock::get_direct_spatial(bool dir)
{
    bool has_direct = (mb.SubMbType[0] == 0) | (mb.SubMbType[1] == 0) |
                      (mb.SubMbType[2] == 0) | (mb.SubMbType[3] == 0);
    if (!has_direct)
        return;

    int list_offset = slice.MbaffFrameFlag && mb.mb_field_decoding_flag ?
                      mb.mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice.listX[LIST_0 + list_offset];
    storable_picture** list1 = slice.listX[LIST_1 + list_offset];

    char l0_rFrame, l1_rFrame;
    mv_t pmvl0 = this->GetMVPredictor2(l0_rFrame, LIST_0, 0, 0, 16, 16);
    mv_t pmvl1 = this->GetMVPredictor2(l1_rFrame, LIST_1, 0, 0, 16, 16);

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

    int blks = sps.direct_8x8_inference_flag ? 4 : 1;
    for (int block4x4 = 0; block4x4 < 16; block4x4 += blks) {
        if (mb.SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);
            if (dir)
                mb.SubMbPredMode[block4x4 / 4] = pred_dir;

            auto mv_info = &slice.dec_picture->mv_info[mb.mb.y * 4 + j][mb.mb.x * 4 + i];
            mv_info->ref_pic[LIST_0] = ref_pic_l[LIST_0];
            mv_info->ref_pic[LIST_1] = ref_pic_l[LIST_1];
            mv_info->ref_idx[LIST_0] = ref_idx_l[LIST_0];
            mv_info->ref_idx[LIST_1] = ref_idx_l[LIST_1];

            int block_y_aff;
            if (slice.MbaffFrameFlag && mb.mb_field_decoding_flag)
                block_y_aff = (mb.mb.y * 4 - 4 * (mb.mbAddrX % 2)) / 2;
            else
                block_y_aff = mb.mb.y * 4;
            bool is_not_moving = (this->get_colocated_info(list1[0], mb.mb.x * 4 + i, block_y_aff + j) == 0);
            mv_info->mv[LIST_0] = l0_rFrame == -1 || (l0_rFrame == 0 && is_not_moving) ? mv_t{0, 0} : pmvl0;
            mv_info->mv[LIST_1] = l1_rFrame == -1 || (l1_rFrame == 0 && is_not_moving) ? mv_t{0, 0} : pmvl1;

            if (sps.direct_8x8_inference_flag) {
                slice.dec_picture->mv_info[mb.mb.y * 4 + j + 0][mb.mb.x * 4 + i + 1] = *mv_info;
                slice.dec_picture->mv_info[mb.mb.y * 4 + j + 1][mb.mb.x * 4 + i + 0] = *mv_info;
                slice.dec_picture->mv_info[mb.mb.y * 4 + j + 1][mb.mb.x * 4 + i + 1] = *mv_info;
            }
        }
    }
}

int Parser::get_inter8x8(mb_t& mb, int block8x8)
{
    slice_t& slice = *mb.p_Slice;

    int list_offset = slice.MbaffFrameFlag && mb.mb_field_decoding_flag ?
                      mb.mbAddrX % 2 ? 4 : 2 : 0;
    storable_picture** list0 = slice.listX[LIST_0 + list_offset];
    storable_picture** list1 = slice.listX[LIST_1 + list_offset];

    int pred_dir = mb.SubMbPredMode[block8x8];

    for (int block4x4 = block8x8 * 4; block4x4 < block8x8 * 4 + 4; ++block4x4) {
        if (mb.SubMbType[block4x4 / 4] == 0) {
            int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
            int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);

            auto mv_info = &slice.dec_picture->mv_info[mb.mb.y * 4 + j][mb.mb.x * 4 + i];
            if (slice.direct_spatial_mv_pred_flag) {
                if (mv_info->ref_idx[LIST_1] == -1)
                    pred_dir = 0;
                else if (mv_info->ref_idx[LIST_0] == -1)
                    pred_dir = 1;
                else
                    pred_dir = 2;
                //mb.SubMbPredMode[block8x8] = pred_dir;
            } else {
                mv_info->ref_pic[LIST_0] = list0[(short)mv_info->ref_idx[LIST_0]];
                mv_info->ref_pic[LIST_1] = list1[(short)mv_info->ref_idx[LIST_1]];
            }
        }
    }

    return pred_dir;
}

    
}
}
