
/*!
 *****************************************************************************
 *
 * \file fmo.c
 *
 * \brief
 *    Support for Flexible mb_t Ordering (FMO)
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *    - Stephan Wenger      stewe@cs.tu-berlin.de
 *    - Karsten Suehring
 ******************************************************************************
 */

#include "global.h"
#include "slice.h"
#include "fmo.h"


static void FmoGenerateType0MapUnitMap(VideoParameters *p_Vid)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t* pps = p_Vid->active_pps;
    int iGroup, j, i = 0;

    do {
        for (iGroup = 0;
             iGroup <= pps->num_slice_groups_minus1 && i < sps->PicSizeInMapUnits;
             i += pps->run_length_minus1[iGroup++] + 1) {
            for (j = 0; j <= pps->run_length_minus1[iGroup] && i + j < sps->PicSizeInMapUnits; j++)
                p_Vid->MapUnitToSliceGroupMap[i+j] = iGroup;
        }
    } while (i < sps->PicSizeInMapUnits);
}

static void FmoGenerateType1MapUnitMap(VideoParameters *p_Vid)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;
    int i;

    for (i = 0; i < sps->PicSizeInMapUnits; i++)
        p_Vid->MapUnitToSliceGroupMap[i] =
            ((i % sps->PicWidthInMbs) +
             (((i / sps->PicWidthInMbs) * (pps->num_slice_groups_minus1 + 1)) / 2))
            % (pps->num_slice_groups_minus1 + 1);
}

static void FmoGenerateType2MapUnitMap(VideoParameters *p_Vid)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;
    int iGroup, i, x, y;

    for (i = 0; i < sps->PicSizeInMapUnits; i++)
        p_Vid->MapUnitToSliceGroupMap[i] = pps->num_slice_groups_minus1;

    for (iGroup = pps->num_slice_groups_minus1 - 1 ; iGroup >= 0; iGroup--) {
        int xTopLeft     = pps->top_left    [iGroup] % sps->PicWidthInMbs;
        int yTopLeft     = pps->top_left    [iGroup] / sps->PicWidthInMbs;
        int xBottomRight = pps->bottom_right[iGroup] % sps->PicWidthInMbs;
        int yBottomRight = pps->bottom_right[iGroup] / sps->PicWidthInMbs;
        for (y = yTopLeft; y <= yBottomRight; y++) {
            for (x = xTopLeft; x <= xBottomRight; x++)
                p_Vid->MapUnitToSliceGroupMap[y * sps->PicWidthInMbs + x] = iGroup;
        }
    }
}

static void FmoGenerateType3MapUnitMap(VideoParameters *p_Vid, slice_t *currSlice)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;
    int i, k;
    int x, y, leftBound, topBound, rightBound, bottomBound;
    int xDir, yDir, mapUnitVacant;

    for (i = 0; i < sps->PicSizeInMapUnits; i++)
        p_Vid->MapUnitToSliceGroupMap[i] = 1;

    x = (sps->PicWidthInMbs       - pps->slice_group_change_direction_flag) / 2;
    y = (sps->PicHeightInMapUnits - pps->slice_group_change_direction_flag) / 2;

    leftBound   = x;
    topBound    = y;
    rightBound  = x;
    bottomBound = y;

    xDir = pps->slice_group_change_direction_flag - 1;
    yDir = pps->slice_group_change_direction_flag;

    for (k = 0; k < currSlice->MapUnitsInSliceGroup0; k += mapUnitVacant) {
        mapUnitVacant = (p_Vid->MapUnitToSliceGroupMap[y * sps->PicWidthInMbs + x] == 1);
        if (mapUnitVacant)
            p_Vid->MapUnitToSliceGroupMap[y * sps->PicWidthInMbs + x] = 0;

        if (xDir == -1 && x == leftBound) {
            leftBound = imax(leftBound - 1, 0);
            x = leftBound;
            xDir = 0;
            yDir = 2 * pps->slice_group_change_direction_flag - 1;
        } else if (xDir == 1 && x == rightBound) {
            rightBound = imin(rightBound + 1, sps->PicWidthInMbs - 1);
            x = rightBound;
            xDir = 0;
            yDir = 1 - 2 * pps->slice_group_change_direction_flag;
        } else if (yDir == -1 && y == topBound) {
            topBound = imax(topBound - 1, 0);
            y = topBound;
            xDir = 1 - 2 * pps->slice_group_change_direction_flag;
            yDir = 0;
        } else if (yDir == 1 && y == bottomBound) {
            bottomBound = imin(bottomBound + 1, sps->PicHeightInMapUnits - 1);
            y = bottomBound;
            xDir = 2 * pps->slice_group_change_direction_flag - 1;
            yDir = 0;
        } else {
            x += xDir;
            y += yDir;
        }
    }
}

static void FmoGenerateType4MapUnitMap(VideoParameters *p_Vid, slice_t *currSlice)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;
    int i;

    uint32_t sizeOfUpperLeftGroup =
        pps->slice_group_change_direction_flag ?
            (sps->PicSizeInMapUnits - currSlice->MapUnitsInSliceGroup0) :
            currSlice->MapUnitsInSliceGroup0;

    for (i = 0; i < sps->PicSizeInMapUnits; i++) {
        if (i < sizeOfUpperLeftGroup)
            p_Vid->MapUnitToSliceGroupMap[i] = pps->slice_group_change_direction_flag;
        else
            p_Vid->MapUnitToSliceGroupMap[i] = 1 - pps->slice_group_change_direction_flag;
    }
}

static void FmoGenerateType5MapUnitMap(VideoParameters *p_Vid, slice_t *currSlice)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;
    int i, j, k = 0;

    uint32_t sizeOfUpperLeftGroup =
        pps->slice_group_change_direction_flag ?
            (sps->PicSizeInMapUnits - currSlice->MapUnitsInSliceGroup0) :
            currSlice->MapUnitsInSliceGroup0;

    for (j = 0; j < sps->PicWidthInMbs; j++) {
        for (i = 0; i < sps->PicHeightInMapUnits; i++) {
            if (k++ < sizeOfUpperLeftGroup)
                p_Vid->MapUnitToSliceGroupMap[i * sps->PicWidthInMbs + j] = pps->slice_group_change_direction_flag;
            else
                p_Vid->MapUnitToSliceGroupMap[i * sps->PicWidthInMbs + j] = 1 - pps->slice_group_change_direction_flag;
        }
    }
}

static void FmoGenerateType6MapUnitMap(VideoParameters *p_Vid)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps; 
    int i;

    for (i = 0; i < sps->PicSizeInMapUnits; i++)
        p_Vid->MapUnitToSliceGroupMap[i] = pps->slice_group_id[i];
}

static int FmoGenerateMapUnitToSliceGroupMap(VideoParameters *p_Vid, slice_t *currSlice)
{
    sps_t *sps = p_Vid->active_sps;
    pps_t *pps = p_Vid->active_pps;

    if (pps->slice_group_map_type == 6) {
        if ((pps->pic_size_in_map_units_minus1 + 1) != sps->PicSizeInMapUnits)
            error("wrong pps->pic_size_in_map_units_minus1 for used SPS and FMO type 6", 500);
    }

    // allocate memory for p_Vid->MapUnitToSliceGroupMap
    if (p_Vid->MapUnitToSliceGroupMap)
        free(p_Vid->MapUnitToSliceGroupMap);
    if ((p_Vid->MapUnitToSliceGroupMap = (int *)malloc ((sps->PicSizeInMapUnits) * sizeof (int))) == NULL) {
        printf ("cannot allocated %d bytes for p_Vid->MapUnitToSliceGroupMap, exit\n", (int) ( (pps->pic_size_in_map_units_minus1+1) * sizeof (int)));
        exit (-1);
    }

    if (pps->num_slice_groups_minus1 == 0) {
        memset(p_Vid->MapUnitToSliceGroupMap, 0, sps->PicSizeInMapUnits * sizeof (int));
        return 0;
    }

    switch (pps->slice_group_map_type) {
    case 0:
        FmoGenerateType0MapUnitMap(p_Vid);
        break;
    case 1:
        FmoGenerateType1MapUnitMap(p_Vid);
        break;
    case 2:
        FmoGenerateType2MapUnitMap(p_Vid);
        break;
    case 3:
        FmoGenerateType3MapUnitMap(p_Vid, currSlice);
        break;
    case 4:
        FmoGenerateType4MapUnitMap(p_Vid, currSlice);
        break;
    case 5:
        FmoGenerateType5MapUnitMap(p_Vid, currSlice);
        break;
    case 6:
        FmoGenerateType6MapUnitMap(p_Vid);
        break;
    default:
        printf("Illegal slice_group_map_type %d , exit \n", (int) pps->slice_group_map_type);
        exit(-1);
    }
    return 0;
}

static int FmoGenerateMbToSliceGroupMap(VideoParameters *p_Vid, slice_t *pSlice)
{
    sps_t *sps = p_Vid->active_sps;
    int i;

    // allocate memory for p_Vid->MbToSliceGroupMap
    if (p_Vid->MbToSliceGroupMap)
        free(p_Vid->MbToSliceGroupMap);

    if ((p_Vid->MbToSliceGroupMap = (int *)malloc (pSlice->PicSizeInMbs * sizeof (int))) == NULL) {
        printf("cannot allocate %d bytes for p_Vid->MbToSliceGroupMap, exit\n", (int) ((pSlice->PicSizeInMbs) * sizeof (int)));
        exit(-1);
    }

    if (sps->frame_mbs_only_flag || pSlice->field_pic_flag) {
        for (i = 0; i < pSlice->PicSizeInMbs; i++)
            p_Vid->MbToSliceGroupMap[i] = p_Vid->MapUnitToSliceGroupMap[i];
    } else if (pSlice->MbaffFrameFlag) {
        for (i = 0; i < pSlice->PicSizeInMbs; i++)
            p_Vid->MbToSliceGroupMap[i] = p_Vid->MapUnitToSliceGroupMap[i/2];
    } else {
        for (i = 0; i < pSlice->PicSizeInMbs; i++)
            p_Vid->MbToSliceGroupMap[i] = p_Vid->MapUnitToSliceGroupMap[
                (i / (2 * sps->PicWidthInMbs)) * sps->PicWidthInMbs
                + (i % sps->PicWidthInMbs)];
    }
    return 0;
}


int fmo_init(VideoParameters *p_Vid, slice_t *pSlice)
{
    FmoGenerateMapUnitToSliceGroupMap(p_Vid, pSlice);
    FmoGenerateMbToSliceGroupMap(p_Vid, pSlice);
    return 0;
}

int FmoFinit(VideoParameters *p_Vid)
{
    if (p_Vid->MbToSliceGroupMap) {
        free(p_Vid->MbToSliceGroupMap);
        p_Vid->MbToSliceGroupMap = NULL;
    }
    if (p_Vid->MapUnitToSliceGroupMap) {
        free(p_Vid->MapUnitToSliceGroupMap);
        p_Vid->MapUnitToSliceGroupMap = NULL;
    }
    return 0;
}

int FmoGetNextMBNr(VideoParameters *p_Vid, int n)
{
    slice_t *currSlice = p_Vid->ppSliceList[0];
    int SliceGroup, i;

    assert(n < currSlice->PicSizeInMbs);
    assert(p_Vid->MbToSliceGroupMap != NULL);

    SliceGroup = p_Vid->MbToSliceGroupMap[n];

    i = n + 1;
    while (i < currSlice->PicSizeInMbs && p_Vid->MbToSliceGroupMap[i] != SliceGroup)
        i++;

    return (i < currSlice->PicSizeInMbs) ? i : -1;
}
