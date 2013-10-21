#include "global.h"
#include "memalloc.h"

using vio::h264::mb_t;
using namespace vio::h264;

#include "erc_api.h"
#include "erc_do.h"
#include "slice.h"
#include "macroblock.h"


void ercInit(VideoParameters *p_Vid, int pic_sizex, int pic_sizey, int flag)
{
    ercClose(p_Vid, p_Vid->erc_errorVar);

    p_Vid->erc_object_list = new objectBuffer_t[(pic_sizex * pic_sizey) >> 6];

    // the error concealment instance is allocated
    ercVariables_t* errorVar = new ercVariables_t;
    errorVar->nOfMBs              = 0;
    errorVar->segments            = NULL;
    errorVar->currSegment         = 0;
    errorVar->yCondition          = NULL;
    errorVar->uCondition          = NULL;
    errorVar->vCondition          = NULL;
    errorVar->prevFrameYCondition = NULL;
    errorVar->concealment         = 1;

    p_Vid->erc_errorVar = errorVar;

    // set error concealment ON
    if (p_Vid->erc_errorVar)
        p_Vid->erc_errorVar->concealment = flag;
}

void ercReset(ercVariables_t *errorVar, int nOfMBs, int numOfSegments, int picSizeX)
{
    if (errorVar && errorVar->concealment) {
        ercSegment_t *segments = NULL;
        // If frame size has been changed
        if (nOfMBs != errorVar->nOfMBs && errorVar->yCondition != NULL) {
            delete []errorVar->yCondition;
            errorVar->yCondition = NULL;
            delete []errorVar->prevFrameYCondition;
            errorVar->prevFrameYCondition = NULL;
            delete []errorVar->uCondition;
            errorVar->uCondition = NULL;
            delete []errorVar->vCondition;
            errorVar->vCondition = NULL;
            delete []errorVar->segments;
            errorVar->segments = NULL;
        }

        // If the structures are uninitialized (first frame, or frame size is changed)
        if (errorVar->yCondition == NULL) {
            errorVar->segments = new ercSegment_t[numOfSegments];
            memset(errorVar->segments, 0, numOfSegments * sizeof(ercSegment_t));
            errorVar->nOfSegments = numOfSegments;
            errorVar->yCondition = new char[4 * nOfMBs];
            errorVar->prevFrameYCondition = new char[4 * nOfMBs];
            errorVar->uCondition = new char[nOfMBs];
            errorVar->vCondition = new char[nOfMBs];
            errorVar->nOfMBs = nOfMBs;
        } else {
            // Store the yCondition struct of the previous frame
            char* tmp = errorVar->prevFrameYCondition;
            errorVar->prevFrameYCondition = errorVar->yCondition;
            errorVar->yCondition = tmp;
        }

        // Reset tables and parameters
        memset(errorVar->yCondition, 0, 4 * nOfMBs * sizeof(*errorVar->yCondition));
        memset(errorVar->uCondition, 0,     nOfMBs * sizeof(*errorVar->uCondition));
        memset(errorVar->vCondition, 0,     nOfMBs * sizeof(*errorVar->vCondition));

        if (errorVar->nOfSegments != numOfSegments) {
            delete []errorVar->segments;
            errorVar->segments = new ercSegment_t[numOfSegments];
            errorVar->nOfSegments = numOfSegments;
        }

        segments = errorVar->segments;
        for (int i = 0; i < errorVar->nOfSegments; i++) {
            segments->startMBPos = 0;
            segments->endMBPos = (short) (nOfMBs - 1);
            (segments++)->fCorrupted = 1; //! mark segments as corrupted
        }

        errorVar->currSegment = 0;
        errorVar->nOfCorruptedSegments = 0;
    }
}

void ercClose(VideoParameters *p_Vid, ercVariables_t *errorVar)
{
    if (errorVar) {
        if (errorVar->yCondition) {
            delete []errorVar->segments;
            delete []errorVar->yCondition;
            delete []errorVar->uCondition;
            delete []errorVar->vCondition;
            delete []errorVar->prevFrameYCondition;
        }
        delete errorVar;
        errorVar = NULL;
    }

    if (p_Vid->erc_object_list) {
        delete []p_Vid->erc_object_list;
        p_Vid->erc_object_list = NULL;
    }
}


static void ercStartSegment(int currMBNum, int segment, unsigned int bitPos, ercVariables_t *errorVar)
{
    if (errorVar && errorVar->concealment) {
        errorVar->currSegmentCorrupted = 0;
        errorVar->segments[segment].fCorrupted = 0;
        errorVar->segments[segment].startMBPos = (short) currMBNum;
    }
}

static void ercStopSegment(int currMBNum, int segment, unsigned int bitPos, ercVariables_t *errorVar)
{
    if (errorVar && errorVar->concealment) {
        errorVar->segments[segment].endMBPos = (short) currMBNum;
        errorVar->currSegment++;
    }
}

static void ercMarkCurrSegmentLost(int picSizeX, ercVariables_t *errorVar)
{
    int current_segment = errorVar->currSegment - 1;

    if (errorVar && errorVar->concealment) {
        if (errorVar->currSegmentCorrupted == 0) {
            errorVar->nOfCorruptedSegments++;
            errorVar->currSegmentCorrupted = 1;
        }

        for (int j = errorVar->segments[current_segment].startMBPos; j <= errorVar->segments[current_segment].endMBPos; j++) {
            errorVar->yCondition[MBNum2YBlock (j, 0, picSizeX)] = ERC_BLOCK_CORRUPTED;
            errorVar->yCondition[MBNum2YBlock (j, 1, picSizeX)] = ERC_BLOCK_CORRUPTED;
            errorVar->yCondition[MBNum2YBlock (j, 2, picSizeX)] = ERC_BLOCK_CORRUPTED;
            errorVar->yCondition[MBNum2YBlock (j, 3, picSizeX)] = ERC_BLOCK_CORRUPTED;
            errorVar->uCondition[j] = ERC_BLOCK_CORRUPTED;
            errorVar->vCondition[j] = ERC_BLOCK_CORRUPTED;
        }
        errorVar->segments[current_segment].fCorrupted = 1;
    }
}

static void ercMarkCurrSegmentOK(int picSizeX, ercVariables_t *errorVar)
{
    int current_segment = errorVar->currSegment - 1;

    if (errorVar && errorVar->concealment) {
        // mark all the Blocks belonging to the segment as OK */
        for (int j = errorVar->segments[current_segment].startMBPos; j <= errorVar->segments[current_segment].endMBPos; j++) {
            errorVar->yCondition[MBNum2YBlock (j, 0, picSizeX)] = ERC_BLOCK_OK;
            errorVar->yCondition[MBNum2YBlock (j, 1, picSizeX)] = ERC_BLOCK_OK;
            errorVar->yCondition[MBNum2YBlock (j, 2, picSizeX)] = ERC_BLOCK_OK;
            errorVar->yCondition[MBNum2YBlock (j, 3, picSizeX)] = ERC_BLOCK_OK;
            errorVar->uCondition[j] = ERC_BLOCK_OK;
            errorVar->vCondition[j] = ERC_BLOCK_OK;
        }
        errorVar->segments[current_segment].fCorrupted = 0;
    }
}


void erc_picture(VideoParameters *p_Vid)
{
    storable_picture* dec_picture = p_Vid->dec_picture;
    sps_t* sps = p_Vid->active_sps;
    frame recfr;
    recfr.p_Vid = p_Vid;
    recfr.yptr = &dec_picture->imgY[0][0];
    if (sps->chroma_format_idc != YUV400) {
        recfr.uptr = &dec_picture->imgUV[0][0][0];
        recfr.vptr = &dec_picture->imgUV[1][0][0];
    }

    //! this is always true at the beginning of a picture
    int ercStartMB = 0;
    int ercSegment = 0;

    //! mark the start of the first segment
    if (!dec_picture->slice_headers[0]->header.MbaffFrameFlag) {
        int i;
        ercStartSegment(0, ercSegment, 0 , p_Vid->erc_errorVar);
        //! generate the segments according to the macroblock map
        for (i = 1; i < (int) dec_picture->slice_headers[0]->header.PicSizeInMbs; ++i) {
            if (p_Vid->mb_data[i].ei_flag != p_Vid->mb_data[i-1].ei_flag) {
                ercStopSegment(i-1, ercSegment, 0, p_Vid->erc_errorVar); //! stop current segment

                //! mark current segment as lost or OK
                if (p_Vid->mb_data[i-1].ei_flag)
                    ercMarkCurrSegmentLost(dec_picture->size_x, p_Vid->erc_errorVar);
                else
                    ercMarkCurrSegmentOK(dec_picture->size_x, p_Vid->erc_errorVar);

                ++ercSegment;  //! next segment
                ercStartSegment(i, ercSegment, 0 , p_Vid->erc_errorVar); //! start new segment
                ercStartMB = i;//! save start MB for this segment
            }
        }
        //! mark end of the last segment
        ercStopSegment(dec_picture->slice_headers[0]->header.PicSizeInMbs-1, ercSegment, 0, p_Vid->erc_errorVar);
        if (p_Vid->mb_data[i-1].ei_flag)
            ercMarkCurrSegmentLost(dec_picture->size_x, p_Vid->erc_errorVar);
        else
            ercMarkCurrSegmentOK(dec_picture->size_x, p_Vid->erc_errorVar);

        //! call the right error concealment function depending on the frame type.
        p_Vid->erc_mvperMB /= dec_picture->slice_headers[0]->header.PicSizeInMbs;

        p_Vid->erc_img = p_Vid;

        if (dec_picture->slice.slice_type == I_slice || dec_picture->slice.slice_type == SI_slice) // I-frame
            ercConcealIntraFrame(p_Vid, &recfr, dec_picture->size_x, dec_picture->size_y, p_Vid->erc_errorVar);
        else
            ercConcealInterFrame(&recfr, p_Vid->erc_object_list, dec_picture->size_x, dec_picture->size_y, p_Vid->erc_errorVar, sps->chroma_format_idc);
    }
}

void ercWriteMBMODEandMV(mb_t* currMB)
{
    VideoParameters *p_Vid = currMB->p_Slice->p_Vid;
    int currMBNum = currMB->mbAddrX;
    storable_picture *dec_picture = p_Vid->dec_picture;
    int mbx = xPosMB(currMBNum, dec_picture->size_x), mby = yPosMB(currMBNum, dec_picture->size_x);
    objectBuffer_t *currRegion, *pRegion;

    currRegion = p_Vid->erc_object_list + (currMBNum<<2);

    if (p_Vid->type != B_slice) { //non-B frame
        for (int i = 0; i < 4; ++i) {
            pRegion             = currRegion + i;
            pRegion->regionMode = currMB->mb_type      == I_16x16 ? REGMODE_INTRA :
                                  currMB->SubMbType[i] == I_4x4   ? REGMODE_INTRA_8x8  :
                                  currMB->SubMbType[i] == 0       ? REGMODE_INTER_COPY :
                                  currMB->SubMbType[i] == 1       ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8;
            if (currMB->SubMbType[i] == 0 || currMB->SubMbType[i] == I_4x4) { // INTRA OR COPY
                pRegion->mv[0] = 0;
                pRegion->mv[1] = 0;
                pRegion->mv[2] = 0;
            } else {
                int ii = 4 * mbx + (i & 0x01) * 2;
                int jj = 4 * mby + (i >> 1  ) * 2;
                if (currMB->SubMbType[i] >= 5 && currMB->SubMbType[i] <= 7) { // SMALL BLOCKS
                    pRegion->mv[0] = (dec_picture->mv_info[jj    ][ii    ].mv[LIST_0].mv_x +
                                      dec_picture->mv_info[jj    ][ii + 1].mv[LIST_0].mv_x +
                                      dec_picture->mv_info[jj + 1][ii    ].mv[LIST_0].mv_x +
                                      dec_picture->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_x + 2) / 4;
                    pRegion->mv[1] = (dec_picture->mv_info[jj][ii].mv[LIST_0].mv_y + dec_picture->mv_info[jj][ii + 1].mv[LIST_0].mv_y + dec_picture->mv_info[jj + 1][ii].mv[LIST_0].mv_y + dec_picture->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_y + 2)/4;
                } else { // 16x16, 16x8, 8x16, 8x8
                    pRegion->mv[0] = dec_picture->mv_info[jj][ii].mv[LIST_0].mv_x;
                    pRegion->mv[1] = dec_picture->mv_info[jj][ii].mv[LIST_0].mv_y;
                }
                currMB->p_Slice->erc_mvperMB += abs(pRegion->mv[0]) + abs(pRegion->mv[1]);
                pRegion->mv[2] = dec_picture->mv_info[jj][ii].ref_idx[LIST_0];
            }
        }
    } else { //B-frame
        for (int i = 0; i < 4; ++i) {
            int ii = 4 * mbx + (i % 2) * 2;
            int jj = 4 * mby + (i / 2) * 2;
            pRegion = currRegion + i;
            pRegion->regionMode = currMB->mb_type      == I_16x16 ? REGMODE_INTRA :
                                  currMB->SubMbType[i] == I_4x4   ? REGMODE_INTRA_8x8 : REGMODE_INTER_PRED_8x8;
            if (currMB->mb_type == I_16x16 || currMB->SubMbType[i] == I_4x4) { // INTRA
                pRegion->mv[0] = 0;
                pRegion->mv[1] = 0;
                pRegion->mv[2] = 0;
            } else {
                int idx = (dec_picture->mv_info[jj][ii].ref_idx[0] < 0) ? 1 : 0;
                pRegion->mv[0] = (dec_picture->mv_info[jj    ][ii    ].mv[idx].mv_x + 
                                  dec_picture->mv_info[jj    ][ii + 1].mv[idx].mv_x + 
                                  dec_picture->mv_info[jj + 1][ii    ].mv[idx].mv_x + 
                                  dec_picture->mv_info[jj + 1][ii + 1].mv[idx].mv_x + 2) / 4;
                pRegion->mv[1] = (dec_picture->mv_info[jj    ][ii    ].mv[idx].mv_y + 
                                  dec_picture->mv_info[jj    ][ii + 1].mv[idx].mv_y + 
                                  dec_picture->mv_info[jj + 1][ii    ].mv[idx].mv_y + 
                                  dec_picture->mv_info[jj + 1][ii + 1].mv[idx].mv_y + 2) / 4;
                currMB->p_Slice->erc_mvperMB += abs(pRegion->mv[0]) + abs(pRegion->mv[1]);
                pRegion->mv[2] = (dec_picture->mv_info[jj][ii].ref_idx[idx]);
            }
        }
    }
}
