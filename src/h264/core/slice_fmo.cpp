#include "global.h"
#include "slice.h"


void slice_t::FmoGenerateType0MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    int i = 0;

    do {
        for (int iGroup = 0;
             iGroup <= pps.num_slice_groups_minus1 && i < sps.PicSizeInMapUnits;
             i += pps.slice_groups[iGroup++].run_length_minus1 + 1) {
            for (int j = 0; j <= pps.slice_groups[iGroup].run_length_minus1 && i + j < sps.PicSizeInMapUnits; ++j)
                mapUnitToSliceGroupMap[i + j] = iGroup;
        }
    } while (i < sps.PicSizeInMapUnits);
}

void slice_t::FmoGenerateType1MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i) {
        mapUnitToSliceGroupMap[i] =
            ((i % sps.PicWidthInMbs) +
             (((i / sps.PicWidthInMbs) * (pps.num_slice_groups_minus1 + 1)) / 2))
            % (pps.num_slice_groups_minus1 + 1);
    }
}

void slice_t::FmoGenerateType2MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        mapUnitToSliceGroupMap[i] = pps.num_slice_groups_minus1;

    for (int iGroup = pps.num_slice_groups_minus1 - 1 ; iGroup >= 0; --iGroup) {
        int xTopLeft     = pps.slice_groups[iGroup].top_left     % sps.PicWidthInMbs;
        int yTopLeft     = pps.slice_groups[iGroup].top_left     / sps.PicWidthInMbs;
        int xBottomRight = pps.slice_groups[iGroup].bottom_right % sps.PicWidthInMbs;
        int yBottomRight = pps.slice_groups[iGroup].bottom_right / sps.PicWidthInMbs;
        for (int y = yTopLeft; y <= yBottomRight; ++y) {
            for (int x = xTopLeft; x <= xBottomRight; ++x)
                mapUnitToSliceGroupMap[y * sps.PicWidthInMbs + x] = iGroup;
        }
    }
}

void slice_t::FmoGenerateType3MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        mapUnitToSliceGroupMap[i] = 1;

    int x = (sps.PicWidthInMbs       - pps.slice_group_change_direction_flag) / 2;
    int y = (sps.PicHeightInMapUnits - pps.slice_group_change_direction_flag) / 2;

    int leftBound   = x;
    int topBound    = y;
    int rightBound  = x;
    int bottomBound = y;

    int xDir = pps.slice_group_change_direction_flag - 1;
    int yDir = pps.slice_group_change_direction_flag;

    int mapUnitVacant = 1;
    for (int k = 0; k < shr.MapUnitsInSliceGroup0; k += mapUnitVacant) {
        mapUnitVacant = mapUnitToSliceGroupMap[y * sps.PicWidthInMbs + x] == 1;
        if (mapUnitVacant)
            mapUnitToSliceGroupMap[y * sps.PicWidthInMbs + x] = 0;

        if (xDir == -1 && x == leftBound) {
            leftBound = max<int>(leftBound - 1, 0);
            x = leftBound;
            xDir = 0;
            yDir = 2 * pps.slice_group_change_direction_flag - 1;
        } else if (xDir == 1 && x == rightBound) {
            rightBound = min<int>(rightBound + 1, sps.PicWidthInMbs - 1);
            x = rightBound;
            xDir = 0;
            yDir = 1 - 2 * pps.slice_group_change_direction_flag;
        } else if (yDir == -1 && y == topBound) {
            topBound = max<int>(topBound - 1, 0);
            y = topBound;
            xDir = 1 - 2 * pps.slice_group_change_direction_flag;
            yDir = 0;
        } else if (yDir == 1 && y == bottomBound) {
            bottomBound = min<int>(bottomBound + 1, sps.PicHeightInMapUnits - 1);
            y = bottomBound;
            xDir = 2 * pps.slice_group_change_direction_flag - 1;
            yDir = 0;
        } else {
            x += xDir;
            y += yDir;
        }
    }
}

void slice_t::FmoGenerateType4MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    uint32_t sizeOfUpperLeftGroup = pps.slice_group_change_direction_flag ?
        (sps.PicSizeInMapUnits - shr.MapUnitsInSliceGroup0) : shr.MapUnitsInSliceGroup0;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i) {
        if (i < sizeOfUpperLeftGroup)
            mapUnitToSliceGroupMap[i] = pps.slice_group_change_direction_flag;
        else
            mapUnitToSliceGroupMap[i] = 1 - pps.slice_group_change_direction_flag;
    }
}

void slice_t::FmoGenerateType5MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    uint32_t sizeOfUpperLeftGroup = pps.slice_group_change_direction_flag ?
        (sps.PicSizeInMapUnits - shr.MapUnitsInSliceGroup0) : shr.MapUnitsInSliceGroup0;

    int k = 0;
    for (int j = 0; j < sps.PicWidthInMbs; ++j) {
        for (int i = 0; i < sps.PicHeightInMapUnits; ++i) {
            if (k++ < sizeOfUpperLeftGroup)
                mapUnitToSliceGroupMap[i * sps.PicWidthInMbs + j] = pps.slice_group_change_direction_flag;
            else
                mapUnitToSliceGroupMap[i * sps.PicWidthInMbs + j] = 1 - pps.slice_group_change_direction_flag;
        }
    }
}

void slice_t::FmoGenerateType6MapUnitMap(slice_t::uint8_v& mapUnitToSliceGroupMap)
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps; 

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        mapUnitToSliceGroupMap[i] = pps.slice_group_id[i];
}

void slice_t::init_slice_group_map()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    slice_t::uint8_v mapUnitToSliceGroupMap = slice_t::uint8_v(sps.PicSizeInMapUnits);

    if (pps.num_slice_groups_minus1 == 0) {
        for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
            mapUnitToSliceGroupMap[i] = 0;
    } else {
        switch (pps.slice_group_map_type) {
        case 0:
            this->FmoGenerateType0MapUnitMap(mapUnitToSliceGroupMap);
            break;
        case 1:
            this->FmoGenerateType1MapUnitMap(mapUnitToSliceGroupMap);
            break;
        case 2:
            this->FmoGenerateType2MapUnitMap(mapUnitToSliceGroupMap);
            break;
        case 3:
            this->FmoGenerateType3MapUnitMap(mapUnitToSliceGroupMap);
            break;
        case 4:
            this->FmoGenerateType4MapUnitMap(mapUnitToSliceGroupMap);
            break;
        case 5:
            this->FmoGenerateType5MapUnitMap(mapUnitToSliceGroupMap);
            break;
        case 6:
            this->FmoGenerateType6MapUnitMap(mapUnitToSliceGroupMap);
            break;
        }
    }

    this->MbToSliceGroupMap = slice_t::uint8_v(shr.PicSizeInMbs);

    if (sps.frame_mbs_only_flag || shr.field_pic_flag) {
        for (int i = 0; i < shr.PicSizeInMbs; ++i)
            this->MbToSliceGroupMap[i] = mapUnitToSliceGroupMap[i];
    } else if (shr.MbaffFrameFlag) {
        for (int i = 0; i < shr.PicSizeInMbs; ++i)
            this->MbToSliceGroupMap[i] = mapUnitToSliceGroupMap[i / 2];
    } else {
        for (int i = 0; i < shr.PicSizeInMbs; ++i)
            this->MbToSliceGroupMap[i] =
                mapUnitToSliceGroupMap[(i / (2 * sps.PicWidthInMbs)) * sps.PicWidthInMbs
                                        + (i % sps.PicWidthInMbs)];
    }
}

int slice_t::NextMbAddress(int n)
{
    shr_t& shr = this->header;

    int i = n + 1;
    while (i < shr.PicSizeInMbs && this->MbToSliceGroupMap[i] != this->MbToSliceGroupMap[n])
        ++i;

    return (i < shr.PicSizeInMbs) ? i : -1;
}


void slice_t::decode_poc()
{
    VideoParameters* p_Vid = this->p_Vid;
    sps_t& sps = *this->active_sps;
    shr_t& shr = this->header;

    switch (sps.pic_order_cnt_type) {
    case 0:
        if (this->IdrPicFlag) {
            p_Vid->prevPicOrderCntMsb = 0;
            p_Vid->prevPicOrderCntLsb = 0;
        } else if (p_Vid->last_has_mmco_5) {
            p_Vid->prevPicOrderCntMsb = 0;
            p_Vid->prevPicOrderCntLsb = !p_Vid->last_pic_bottom_field ? shr.TopFieldOrderCnt : 0;
        }

        if ((shr.pic_order_cnt_lsb < p_Vid->prevPicOrderCntLsb) &&
            (p_Vid->prevPicOrderCntLsb - shr.pic_order_cnt_lsb >= sps.MaxPicOrderCntLsb / 2))
            shr.PicOrderCntMsb = p_Vid->prevPicOrderCntMsb + sps.MaxPicOrderCntLsb;
        else if ((shr.pic_order_cnt_lsb > p_Vid->prevPicOrderCntLsb) &&
                 (shr.pic_order_cnt_lsb - p_Vid->prevPicOrderCntLsb > sps.MaxPicOrderCntLsb / 2))
            shr.PicOrderCntMsb = p_Vid->prevPicOrderCntMsb - sps.MaxPicOrderCntLsb;
        else
            shr.PicOrderCntMsb = p_Vid->prevPicOrderCntMsb;

        shr.TopFieldOrderCnt    = 0;
        shr.BottomFieldOrderCnt = 0;
        if (!shr.field_pic_flag || !shr.bottom_field_flag)
            shr.TopFieldOrderCnt = shr.PicOrderCntMsb + shr.pic_order_cnt_lsb;
        if (!shr.field_pic_flag)
            shr.BottomFieldOrderCnt = shr.TopFieldOrderCnt + shr.delta_pic_order_cnt_bottom;
        else if (shr.bottom_field_flag)
            shr.BottomFieldOrderCnt = shr.PicOrderCntMsb + shr.pic_order_cnt_lsb;

        if (this->nal_ref_idc) {
            p_Vid->prevPicOrderCntMsb = shr.PicOrderCntMsb;
            p_Vid->prevPicOrderCntLsb = shr.pic_order_cnt_lsb;
        }
        break;

    case 1:
        if (p_Vid->last_has_mmco_5) {
            p_Vid->prevFrameNum       = 0;
            p_Vid->prevFrameNumOffset = 0;
        }

        if (this->IdrPicFlag)
            shr.FrameNumOffset = 0;
        else if (p_Vid->prevFrameNum > shr.frame_num)
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset + sps.MaxFrameNum;
        else
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset;

        int32_t absFrameNum;
        if (sps.num_ref_frames_in_pic_order_cnt_cycle != 0)
            absFrameNum = shr.FrameNumOffset + shr.frame_num;
        else
            absFrameNum = 0;
        if (this->nal_ref_idc == 0 && absFrameNum > 0)
            --absFrameNum;

        int32_t expectedPicOrderCnt;
        if (absFrameNum > 0) {
            int32_t picOrderCntCycleCnt        = (absFrameNum - 1) / sps.num_ref_frames_in_pic_order_cnt_cycle;
            int32_t frameNumInPicOrderCntCycle = (absFrameNum - 1) % sps.num_ref_frames_in_pic_order_cnt_cycle;
            expectedPicOrderCnt = picOrderCntCycleCnt * sps.ExpectedDeltaPerPicOrderCntCycle;
            for (int i = 0; i <= frameNumInPicOrderCntCycle; i++)
                expectedPicOrderCnt += sps.offset_for_ref_frame[i];
        } else
            expectedPicOrderCnt = 0;
        if (this->nal_ref_idc == 0)
            expectedPicOrderCnt += sps.offset_for_non_ref_pic;

        shr.TopFieldOrderCnt    = 0;
        shr.BottomFieldOrderCnt = 0;
        if (!shr.field_pic_flag || !shr.bottom_field_flag)
            shr.TopFieldOrderCnt = expectedPicOrderCnt + shr.delta_pic_order_cnt[0];
        if (!shr.field_pic_flag)
            shr.BottomFieldOrderCnt = shr.TopFieldOrderCnt + sps.offset_for_top_to_bottom_field + shr.delta_pic_order_cnt[1];
        else if (shr.bottom_field_flag)
            shr.BottomFieldOrderCnt = expectedPicOrderCnt + sps.offset_for_top_to_bottom_field + shr.delta_pic_order_cnt[0];

        p_Vid->prevFrameNum       = shr.frame_num;
        p_Vid->prevFrameNumOffset = shr.FrameNumOffset;
        break;

    case 2:
        if (p_Vid->last_has_mmco_5) {
            p_Vid->prevFrameNum       = 0;
            p_Vid->prevFrameNumOffset = 0;
        }

        if (this->IdrPicFlag)
            shr.FrameNumOffset = 0;
        else if (p_Vid->prevFrameNum > shr.frame_num)
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset + sps.MaxFrameNum;
        else
            shr.FrameNumOffset = p_Vid->prevFrameNumOffset;

        int tempPicOrderCnt;
        if (this->IdrPicFlag)
            tempPicOrderCnt = 0;
        else if (this->nal_ref_idc == 0)
            tempPicOrderCnt = 2 * (shr.FrameNumOffset + shr.frame_num) - 1;
        else
            tempPicOrderCnt = 2 * (shr.FrameNumOffset + shr.frame_num);

        shr.TopFieldOrderCnt    = 0;
        shr.BottomFieldOrderCnt = 0;
        if (!shr.field_pic_flag || !shr.bottom_field_flag)
            shr.TopFieldOrderCnt = tempPicOrderCnt;
        if (!shr.field_pic_flag || shr.bottom_field_flag)
            shr.BottomFieldOrderCnt = tempPicOrderCnt;

        p_Vid->prevFrameNum       = shr.frame_num;
        p_Vid->prevFrameNumOffset = shr.FrameNumOffset;
        break;

    default:
        assert(false);
        break;
    }

    if (!shr.field_pic_flag)
        shr.PicOrderCnt = min(shr.TopFieldOrderCnt, shr.BottomFieldOrderCnt);
    else if (!shr.bottom_field_flag)
        shr.PicOrderCnt = shr.TopFieldOrderCnt;
    else
        shr.PicOrderCnt = shr.BottomFieldOrderCnt;
}
