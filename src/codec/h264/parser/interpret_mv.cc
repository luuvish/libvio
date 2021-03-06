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


struct nb_mv_t {
    bool available;
    bool predFlagL;
    int  refIdxL;
    mv_t mvL;
};

void neighbour_mv(mb_t& mb, nb_mv_t nb_mv[3], int listSuffixFlag, int i, int j, int step_h4, int step_v4)
{
    slice_t& slice = *mb.p_Slice;

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

    int refIdxLXA = -1;
    mv_t mvLXA = {0, 0};
    if (nbA.mb && !nbA.mb->is_intra_block) {
        auto mv_info = &slice.dec_picture->mv_info[nbA.y / 4][nbA.x / 4];
        refIdxLXA = mv_info->ref_idx[listSuffixFlag];
        mvLXA     = mv_info->mv     [listSuffixFlag];
        if (mb.mb_field_decoding_flag && !nbA.mb->mb_field_decoding_flag) {
            if (refIdxLXA >= 0)
                refIdxLXA  *= 2;
            mvLXA.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbA.mb->mb_field_decoding_flag) {
            refIdxLXA >>= 1;
            mvLXA.mv_y *= 2;
        }
    }

    int refIdxLXB = -1;
    mv_t mvLXB = {0, 0};
    if (nbB.mb && !nbB.mb->is_intra_block) {
        auto mv_info = &slice.dec_picture->mv_info[nbB.y / 4][nbB.x / 4];
        refIdxLXB = mv_info->ref_idx[listSuffixFlag];
        mvLXB     = mv_info->mv     [listSuffixFlag];
        if (mb.mb_field_decoding_flag && !nbB.mb->mb_field_decoding_flag) {
            if (refIdxLXB >= 0)
                refIdxLXB  *= 2;
            mvLXB.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbB.mb->mb_field_decoding_flag) {
            refIdxLXB >>= 1;
            mvLXB.mv_y *= 2;
        }
    }

    int refIdxLXC = -1;
    mv_t mvLXC = {0, 0};
    if (nbC.mb && !nbC.mb->is_intra_block) {
        auto mv_info = &slice.dec_picture->mv_info[nbC.y / 4][nbC.x / 4];
        refIdxLXC = mv_info->ref_idx[listSuffixFlag];
        mvLXC     = mv_info->mv     [listSuffixFlag];
        if (mb.mb_field_decoding_flag && !nbC.mb->mb_field_decoding_flag) {
            if (refIdxLXC >= 0)
                refIdxLXC  *= 2;
            mvLXC.mv_y /= 2;
        }
        if (!mb.mb_field_decoding_flag && nbC.mb->mb_field_decoding_flag) {
            refIdxLXC >>= 1;
            mvLXC.mv_y *= 2;
        }
    }

    nb_mv[0].available = nbA.mb ? true : false;
    nb_mv[0].refIdxL   = refIdxLXA;
    nb_mv[0].mvL       = mvLXA;
    nb_mv[1].available = nbB.mb ? true : false;
    nb_mv[1].refIdxL   = refIdxLXB;
    nb_mv[1].mvL       = mvLXB;
    nb_mv[2].available = nbC.mb ? true : false;
    nb_mv[2].refIdxL   = refIdxLXC;
    nb_mv[2].mvL       = mvLXC;
}

mv_t predict_mv(nb_mv_t nb_mv[3], int refIdxLX, int i, int j, int step_h4, int step_v4)
{
    const int&  refIdxLXA = nb_mv[0].refIdxL;
    const int&  refIdxLXB = nb_mv[1].refIdxL;
    const int&  refIdxLXC = nb_mv[2].refIdxL;
    const mv_t& mvLXA     = nb_mv[0].mvL;
    const mv_t& mvLXB     = nb_mv[1].mvL;
    const mv_t& mvLXC     = nb_mv[2].mvL;

    mv_t mvpLX;

    if      (step_h4 == 8 && step_v4 == 16 && i == 0 && refIdxLXA == refIdxLX)
        mvpLX = mvLXA;
    else if (step_h4 == 8 && step_v4 == 16 && i != 0 && refIdxLXC == refIdxLX)
        mvpLX = mvLXC;
    else if (step_h4 == 16 && step_v4 == 8 && j == 0 && refIdxLXB == refIdxLX)
        mvpLX = mvLXB;
    else if (step_h4 == 16 && step_v4 == 8 && j != 0 && refIdxLXA == refIdxLX)
        mvpLX = mvLXA;
    else if (nb_mv[0].available && !nb_mv[1].available && !nb_mv[2].available)
        mvpLX = mvLXA;
    else if (refIdxLXA == refIdxLX && refIdxLXB != refIdxLX && refIdxLXC != refIdxLX)
        mvpLX = mvLXA;
    else if (refIdxLXA != refIdxLX && refIdxLXB == refIdxLX && refIdxLXC != refIdxLX)
        mvpLX = mvLXB;
    else if (refIdxLXA != refIdxLX && refIdxLXB != refIdxLX && refIdxLXC == refIdxLX)
        mvpLX = mvLXC;
    else {
        mvpLX.mv_x = median(mvLXA.mv_x, mvLXB.mv_x, mvLXC.mv_x);
        mvpLX.mv_y = median(mvLXA.mv_y, mvLXB.mv_y, mvLXC.mv_y);
    }

    return mvpLX;
}

mv_t Parser::Macroblock::GetMVPredictor(char refIdxLX, int list, int i, int j, int step_h4, int step_v4)
{
    nb_mv_t nb_mv[3];
    neighbour_mv(mb, nb_mv, list, i, j, step_h4, step_v4);
    return predict_mv(nb_mv, refIdxLX, i, j, step_h4, step_v4);
}

void Parser::Macroblock::skip_macroblock()
{
    int list = LIST_0;
    int refIdxLX = 0;

    nb_mv_t nb_mv[3];
    neighbour_mv(mb, nb_mv, list, 0, 0, 16, 16);
    int refIdxLXA = nb_mv[0].refIdxL;
    int refIdxLXB = nb_mv[1].refIdxL;
    const mv_t& mvLXA = nb_mv[0].mvL;
    const mv_t& mvLXB = nb_mv[1].mvL;

    mv_t mvpLX = {0, 0};
    if (!(!nb_mv[0].available || (refIdxLXA == 0 && mvLXA == mv_t{0, 0}) ||
          !nb_mv[1].available || (refIdxLXB == 0 && mvLXB == mv_t{0, 0})))
        mvpLX = predict_mv(nb_mv, refIdxLX, 0, 0, 16, 16);

    storable_picture* ref_pic = get_ref_pic(mb, slice.RefPicList[0], refIdxLX);

    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            auto mv_info = &slice.dec_picture->mv_info[mb.mb.y * 4 + y][mb.mb.x * 4 + x];
            mv_info->ref_pic[list] = ref_pic;
            mv_info->ref_idx[list] = refIdxLX;
            mv_info->mv     [list] = mvpLX;
        }
    }

    mb.CodedBlockPatternLuma   = 0;
    mb.CodedBlockPatternChroma = 0;
    if (!pps.entropy_coding_mode_flag)
        memset(mb.nz_coeff, 0, 3 * 16 * sizeof(uint8_t));
}

pic_motion_params* get_colocated(mb_t& mb, int i, int j)
{
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;

    storable_picture* ref_pic = get_ref_pic(mb, slice.RefPicList[LIST_1], 0);

    int block_y_aff;
    if (shr.MbaffFrameFlag && mb.mb_field_decoding_flag)
        block_y_aff = (mb.mb.y * 4 - 4 * (mb.mbAddrX % 2)) / 2;
    else
        block_y_aff = mb.mb.y * 4;

    int i4 = mb.mb.x * 4 + i;
    int j4 = block_y_aff + j;

    pic_motion_params* colocated;
    if (sps.direct_8x8_inference_flag)
        colocated = &ref_pic->mv_info[RSD(j4)][RSD(i4)];
    else
        colocated = &ref_pic->mv_info[j4][i4];

    if (sps.direct_8x8_inference_flag) {
        if (shr.MbaffFrameFlag) {
            if (!mb.mb_field_decoding_flag &&
                (ref_pic->slice.iCodingType == FIELD_CODING || ref_pic->motion.mb_field_decoding_flag[mb.mbAddrX])) {
                if (abs(slice.dec_picture->poc - ref_pic->bottom_field->poc) >
                    abs(slice.dec_picture->poc - ref_pic->top_field->poc))
                    colocated = &ref_pic->top_field->mv_info[RSD(j4) >> 1][RSD(i4)];
                else
                    colocated = &ref_pic->bottom_field->mv_info[RSD(j4) >> 1][RSD(i4)];
            }
        } else if (!sps.frame_mbs_only_flag && !shr.field_pic_flag && ref_pic->slice.iCodingType == FIELD_CODING) {
            if (abs(slice.dec_picture->poc - ref_pic->bottom_field->poc) >
                abs(slice.dec_picture->poc - ref_pic->top_field->poc) )
                colocated = &ref_pic->top_field->mv_info[RSD(j4) >> 1][RSD(i4)];
            else
                colocated = &ref_pic->bottom_field->mv_info[RSD(j4) >> 1][RSD(i4)];
        } else if (shr.field_pic_flag && ref_pic->slice.iCodingType != FIELD_CODING) {
            if (!shr.bottom_field_flag)
                colocated = &ref_pic->frame->top_field->mv_info[RSD(j4)][RSD(i4)];
            else
                colocated = &ref_pic->frame->bottom_field->mv_info[RSD(j4)][RSD(i4)];
        }
    }

    return colocated;
}

static int MapColToList0(mb_t& mb, pic_motion_params* colocated)
{
    slice_t& slice = *mb.p_Slice;
    sps_t& sps = *slice.active_sps;
    shr_t& shr = slice.header;

    int num_ref_list = min<int>(shr.num_ref_idx_l0_active_minus1 + 1,
                                slice.RefPicSize[LIST_0] * (1 + (shr.MbaffFrameFlag && mb.mb_field_decoding_flag)));

    int  refList = colocated->ref_idx[LIST_0] == -1 ? LIST_1 : LIST_0;
    auto ref_pic = colocated->ref_pic[refList];

    bool direct_8x8 = sps.direct_8x8_inference_flag && (
        (shr.MbaffFrameFlag && !mb.mb_field_decoding_flag && ref_pic->slice.structure != FRAME) ||
        (!shr.MbaffFrameFlag && !shr.field_pic_flag && ref_pic->slice.structure != FRAME) ||
        (shr.MbaffFrameFlag && mb.mb_field_decoding_flag && ref_pic->slice.structure == FRAME) ||
        (!shr.MbaffFrameFlag && shr.field_pic_flag && ref_pic->slice.structure == FRAME));

    int mapped_idx = 0;
    for (int iref = 0; iref < num_ref_list; ++iref) {
        storable_picture* ref_pic_i = get_ref_pic(mb, slice.RefPicList[LIST_0], iref);
        if (direct_8x8) {
            bool same_pic = ref_pic_i->top_field == ref_pic || 
                            ref_pic_i->bottom_field == ref_pic ||
                            ref_pic_i->frame == ref_pic;

            if (same_pic && (!shr.field_pic_flag || ref_pic_i->slice.structure == shr.structure)) {
                mapped_idx = iref;
                break;
            } else
                mapped_idx = INVALIDINDEX;
        } else {
            if (ref_pic_i == ref_pic) {
                mapped_idx = iref;            
                break;
            } else
                mapped_idx = INVALIDINDEX;
        }
    }
    if (INVALIDINDEX == mapped_idx)
        error(-1111, "temporal direct error: colocated block has ref that is unavailable");

    return mapped_idx;
}

static int DistScaleFactor(mb_t& mb, int ref_idx)
{
    slice_t& slice = *mb.p_Slice;
    shr_t& shr = slice.header;

    storable_picture* dec_picture = slice.p_Vid->dec_picture;
    storable_picture* ref_pic0 = get_ref_pic(mb, slice.RefPicList[LIST_0], ref_idx);
    storable_picture* ref_pic1 = get_ref_pic(mb, slice.RefPicList[LIST_1], 0);

    if (ref_pic0->is_long_term)
        return 9999;

    int curPoc = (!shr.MbaffFrameFlag || !mb.mb_field_decoding_flag) ? dec_picture->poc :
                 (mb.mbAddrX % 2 == 0) ? dec_picture->top_poc : dec_picture->bottom_poc;

    int tb = clip3(-128, 127, curPoc - ref_pic0->poc);
    int td = clip3(-128, 127, ref_pic1->poc - ref_pic0->poc);
    if (td == 0)
        return 9999;
    int tx = (16384 + abs(td / 2)) / td;

    int DistScaleFactor = clip3(-1024, 1023, (tb * tx + 32) >> 6);

    return DistScaleFactor;
}

void Parser::Macroblock::get_direct_temporal()
{
    bool has_direct = (mb.SubMbType[0] == 0) | (mb.SubMbType[1] == 0) |
                      (mb.SubMbType[2] == 0) | (mb.SubMbType[3] == 0);
    if (!has_direct)
        return;

    shr_t& shr = slice.header;

    for (int block4x4 = 0; block4x4 < 16; ++block4x4) {
        if (mb.SubMbType[block4x4 / 4] != 0)
            continue;
        mb.SubMbPredMode[block4x4 / 4] = 2;

        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);

        auto colocated = vio::h264::get_colocated(mb, i, j);
        int  refList = colocated->ref_idx[LIST_0] == -1 ? LIST_1 : LIST_0;
        int  ref_idx = colocated->ref_idx[refList];

        auto mv_info = &slice.dec_picture->mv_info[mb.mb.y * 4 + j][mb.mb.x * 4 + i];
        if (ref_idx == -1) { // co-located is intra mode
            mv_info->ref_idx[LIST_0] = 0;
            mv_info->mv[LIST_0] = {0, 0};
            mv_info->mv[LIST_1] = {0, 0};
        } else { // co-located skip or inter mode
            auto ref_pic = colocated->ref_pic[refList];
            mv_t mvCol   = colocated->mv     [refList];
            if (sps.direct_8x8_inference_flag) {
                if ((shr.MbaffFrameFlag && !mb.mb_field_decoding_flag && ref_pic->slice.structure != FRAME) ||
                    (!shr.MbaffFrameFlag && !shr.field_pic_flag && ref_pic->slice.structure != FRAME))
                    mvCol.mv_y *= 2;
                if ((shr.MbaffFrameFlag && mb.mb_field_decoding_flag && ref_pic->slice.structure == FRAME) ||
                    (!shr.MbaffFrameFlag && shr.field_pic_flag && ref_pic->slice.structure == FRAME))
                    mvCol.mv_y /= 2;
            }

            int mapped_idx = MapColToList0(mb, colocated);
            int mv_scale = DistScaleFactor(mb, mapped_idx);
            mv_info->ref_idx[LIST_0] = (char) mapped_idx;
            //! In such case, an array is needed for each different reference.
            if (mv_scale == 9999) {
                mv_info->mv[LIST_0] = mvCol;
                mv_info->mv[LIST_1] = {0, 0};
            } else {
                mv_info->mv[LIST_0].mv_x = (mv_scale * mvCol.mv_x + 128) >> 8;
                mv_info->mv[LIST_0].mv_y = (mv_scale * mvCol.mv_y + 128) >> 8;
                mv_info->mv[LIST_1] = mv_info->mv[LIST_0] - mvCol;
            }
        }
        // store reference picture ID determined by direct mode
        mv_info->ref_idx[LIST_1] = 0;
        mv_info->ref_pic[LIST_0] = get_ref_pic(mb, slice.RefPicList[LIST_0], (short)mv_info->ref_idx[LIST_0]);
        mv_info->ref_pic[LIST_1] = get_ref_pic(mb, slice.RefPicList[LIST_1], (short)mv_info->ref_idx[LIST_1]);
    }
}

void Parser::Macroblock::get_direct_spatial()
{
    bool has_direct = (mb.SubMbType[0] == 0) | (mb.SubMbType[1] == 0) |
                      (mb.SubMbType[2] == 0) | (mb.SubMbType[3] == 0);
    if (!has_direct)
        return;

    nb_mv_t nb_mv_l0[3];
    nb_mv_t nb_mv_l1[3];
    neighbour_mv(mb, nb_mv_l0, LIST_0, 0, 0, 16, 16);
    neighbour_mv(mb, nb_mv_l1, LIST_1, 0, 0, 16, 16);
    int refIdxL0A = nb_mv_l0[0].refIdxL;
    int refIdxL0B = nb_mv_l0[1].refIdxL;
    int refIdxL0C = nb_mv_l0[2].refIdxL;
    int refIdxL1A = nb_mv_l1[0].refIdxL;
    int refIdxL1B = nb_mv_l1[1].refIdxL;
    int refIdxL1C = nb_mv_l1[2].refIdxL;

    int refIdxL0 = min_positive(refIdxL0A, refIdxL0B, refIdxL0C);
    int refIdxL1 = min_positive(refIdxL1A, refIdxL1B, refIdxL1C);

    bool directZeroPredictionFlag = refIdxL0 < 0 && refIdxL1 < 0;
    if (directZeroPredictionFlag) {
        refIdxL0 = 0;
        refIdxL1 = 0;
    }

    mv_t pmvl0 = predict_mv(nb_mv_l0, refIdxL0, 0, 0, 16, 16);
    mv_t pmvl1 = predict_mv(nb_mv_l1, refIdxL1, 0, 0, 16, 16);

    int blks = sps.direct_8x8_inference_flag ? 4 : 1;
    for (int block4x4 = 0; block4x4 < 16; block4x4 += blks) {
        if (mb.SubMbType[block4x4 / 4] != 0)
            continue;
        mb.SubMbPredMode[block4x4 / 4] = refIdxL1 < 0 ? 0 : refIdxL0 < 0 ? 1 : 2;

        int i = ((block4x4 / 4) % 2) * 2 + ((block4x4 % 4) % 2);
        int j = ((block4x4 / 4) / 2) * 2 + ((block4x4 % 4) / 2);

        bool colZeroFlag = false;
        if (!get_ref_pic(mb, slice.RefPicList[LIST_1], 0)->is_long_term) {
            auto colocated = vio::h264::get_colocated(mb, i, j);
            colZeroFlag =
                (colocated->ref_idx[LIST_0] == 0 &&
                   abs(colocated->mv[LIST_0].mv_x) >> 1 == 0 && abs(colocated->mv[LIST_0].mv_y) >> 1 == 0) ||
                (colocated->ref_idx[LIST_0] == -1 && colocated->ref_idx[LIST_1] == 0 &&
                   abs(colocated->mv[LIST_1].mv_x) >> 1 == 0 && abs(colocated->mv[LIST_1].mv_y) >> 1 == 0);
        }

        auto mv_info = &slice.dec_picture->mv_info[mb.mb.y * 4 + j][mb.mb.x * 4 + i];
        mv_info->ref_pic[LIST_0] = refIdxL0 == -1 ? nullptr : get_ref_pic(mb, slice.RefPicList[LIST_0], refIdxL0);
        mv_info->ref_pic[LIST_1] = refIdxL1 == -1 ? nullptr : get_ref_pic(mb, slice.RefPicList[LIST_1], refIdxL1);
        mv_info->ref_idx[LIST_0] = refIdxL0;
        mv_info->ref_idx[LIST_1] = refIdxL1;
        mv_info->mv[LIST_0] = (directZeroPredictionFlag || refIdxL0 < 0 || (refIdxL0 == 0 && colZeroFlag)) ? mv_t{0, 0} : pmvl0;
        mv_info->mv[LIST_1] = (directZeroPredictionFlag || refIdxL1 < 0 || (refIdxL1 == 0 && colZeroFlag)) ? mv_t{0, 0} : pmvl1;

        if (sps.direct_8x8_inference_flag) {
            slice.dec_picture->mv_info[mb.mb.y * 4 + j + 0][mb.mb.x * 4 + i + 1] = *mv_info;
            slice.dec_picture->mv_info[mb.mb.y * 4 + j + 1][mb.mb.x * 4 + i + 0] = *mv_info;
            slice.dec_picture->mv_info[mb.mb.y * 4 + j + 1][mb.mb.x * 4 + i + 1] = *mv_info;
        }
    }
}

    
}
}
