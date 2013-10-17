#include "global.h"
#include "slice.h"


void slice_t::FmoGenerateType0MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;
    int i = 0;

    do {
        for (int iGroup = 0;
             iGroup <= pps.num_slice_groups_minus1 && i < sps.PicSizeInMapUnits;
             i += pps.slice_groups[iGroup++].run_length_minus1 + 1) {
            for (int j = 0; j <= pps.slice_groups[iGroup].run_length_minus1 && i + j < sps.PicSizeInMapUnits; ++j)
                shr.MapUnitToSliceGroupMap[i + j] = iGroup;
        }
    } while (i < sps.PicSizeInMapUnits);
}

void slice_t::FmoGenerateType1MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        shr.MapUnitToSliceGroupMap[i] =
            ((i % sps.PicWidthInMbs) +
             (((i / sps.PicWidthInMbs) * (pps.num_slice_groups_minus1 + 1)) / 2))
            % (pps.num_slice_groups_minus1 + 1);
}

void slice_t::FmoGenerateType2MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        shr.MapUnitToSliceGroupMap[i] = pps.num_slice_groups_minus1;

    for (int iGroup = pps.num_slice_groups_minus1 - 1 ; iGroup >= 0; --iGroup) {
        int xTopLeft     = pps.slice_groups[iGroup].top_left     % sps.PicWidthInMbs;
        int yTopLeft     = pps.slice_groups[iGroup].top_left     / sps.PicWidthInMbs;
        int xBottomRight = pps.slice_groups[iGroup].bottom_right % sps.PicWidthInMbs;
        int yBottomRight = pps.slice_groups[iGroup].bottom_right / sps.PicWidthInMbs;
        for (int y = yTopLeft; y <= yBottomRight; ++y) {
            for (int x = xTopLeft; x <= xBottomRight; ++x)
                shr.MapUnitToSliceGroupMap[y * sps.PicWidthInMbs + x] = iGroup;
        }
    }
}

void slice_t::FmoGenerateType3MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        shr.MapUnitToSliceGroupMap[i] = 1;

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
        mapUnitVacant = (shr.MapUnitToSliceGroupMap[y * sps.PicWidthInMbs + x] == 1);
        if (mapUnitVacant)
            shr.MapUnitToSliceGroupMap[y * sps.PicWidthInMbs + x] = 0;

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

void slice_t::FmoGenerateType4MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    uint32_t sizeOfUpperLeftGroup =
        pps.slice_group_change_direction_flag ?
            (sps.PicSizeInMapUnits - shr.MapUnitsInSliceGroup0) :
            shr.MapUnitsInSliceGroup0;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i) {
        if (i < sizeOfUpperLeftGroup)
            shr.MapUnitToSliceGroupMap[i] = pps.slice_group_change_direction_flag;
        else
            shr.MapUnitToSliceGroupMap[i] = 1 - pps.slice_group_change_direction_flag;
    }
}

void slice_t::FmoGenerateType5MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    uint32_t sizeOfUpperLeftGroup =
        pps.slice_group_change_direction_flag ?
            (sps.PicSizeInMapUnits - shr.MapUnitsInSliceGroup0) :
            shr.MapUnitsInSliceGroup0;

    int k = 0;
    for (int j = 0; j < sps.PicWidthInMbs; ++j) {
        for (int i = 0; i < sps.PicHeightInMapUnits; ++i) {
            if (k++ < sizeOfUpperLeftGroup)
                shr.MapUnitToSliceGroupMap[i * sps.PicWidthInMbs + j] = pps.slice_group_change_direction_flag;
            else
                shr.MapUnitToSliceGroupMap[i * sps.PicWidthInMbs + j] = 1 - pps.slice_group_change_direction_flag;
        }
    }
}

void slice_t::FmoGenerateType6MapUnitMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps; 
    shr_t& shr = this->header;

    for (int i = 0; i < sps.PicSizeInMapUnits; ++i)
        shr.MapUnitToSliceGroupMap[i] = pps.slice_group_id[i];
}

void slice_t::FmoGenerateMapUnitToSliceGroupMap()
{
    sps_t& sps = *this->active_sps;
    pps_t& pps = *this->active_pps;
    shr_t& shr = this->header;

    if (shr.MapUnitToSliceGroupMap)
        delete []shr.MapUnitToSliceGroupMap;
    shr.MapUnitToSliceGroupMap = new uint8_t[sps.PicSizeInMapUnits];

    if (pps.num_slice_groups_minus1 == 0) {
        memset(shr.MapUnitToSliceGroupMap, 0, sps.PicSizeInMapUnits * sizeof(uint8_t));
        return;
    }

    switch (pps.slice_group_map_type) {
    case 0:
        return this->FmoGenerateType0MapUnitMap();
    case 1:
        return this->FmoGenerateType1MapUnitMap();
    case 2:
        return this->FmoGenerateType2MapUnitMap();
    case 3:
        return this->FmoGenerateType3MapUnitMap();
    case 4:
        return this->FmoGenerateType4MapUnitMap();
    case 5:
        return this->FmoGenerateType5MapUnitMap();
    case 6:
        return this->FmoGenerateType6MapUnitMap();
    }
}

void slice_t::FmoGenerateMbToSliceGroupMap()
{
    sps_t& sps = *this->active_sps;
    shr_t& shr = this->header;

    if (shr.MbToSliceGroupMap)
        delete []shr.MbToSliceGroupMap;
    shr.MbToSliceGroupMap = new uint8_t[shr.PicSizeInMbs];

    if (sps.frame_mbs_only_flag || shr.field_pic_flag) {
        for (int i = 0; i < shr.PicSizeInMbs; ++i)
            shr.MbToSliceGroupMap[i] = shr.MapUnitToSliceGroupMap[i];
    } else if (shr.MbaffFrameFlag) {
        for (int i = 0; i < shr.PicSizeInMbs; ++i)
            shr.MbToSliceGroupMap[i] = shr.MapUnitToSliceGroupMap[i / 2];
    } else {
        for (int i = 0; i < shr.PicSizeInMbs; ++i)
            shr.MbToSliceGroupMap[i] = shr.MapUnitToSliceGroupMap[
                (i / (2 * sps.PicWidthInMbs)) * sps.PicWidthInMbs
                + (i % sps.PicWidthInMbs)];
    }
}


void slice_t::fmo_init()
{
    this->FmoGenerateMapUnitToSliceGroupMap();
    this->FmoGenerateMbToSliceGroupMap();
}

void slice_t::fmo_close()
{
    shr_t& shr = this->header;

    if (shr.MbToSliceGroupMap) {
        delete []shr.MbToSliceGroupMap;
        shr.MbToSliceGroupMap = nullptr;
    }

    if (shr.MapUnitToSliceGroupMap) {
        delete []shr.MapUnitToSliceGroupMap;
        shr.MapUnitToSliceGroupMap = nullptr;
    }
}

int slice_t::FmoGetNextMBNr(int n)
{
    shr_t& shr = this->header;

    assert(n < shr.PicSizeInMbs);
    assert(shr.MbToSliceGroupMap);

    int i = n + 1;
    int SliceGroup = shr.MbToSliceGroupMap[n];
    while (i < shr.PicSizeInMbs && shr.MbToSliceGroupMap[i] != SliceGroup)
        ++i;

    return (i < shr.PicSizeInMbs) ? i : -1;
}
