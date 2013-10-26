#include "global.h"
#include "slice.h"
#include "macroblock.h"
#include "erc_api.h"

using namespace vio::h264;


//! region structure stores information about a region that is needed for concealment
struct objectBuffer_t {
    uint8_t     regionMode;  //!< region mode as above
    int         xMin;        //!< X coordinate of the pixel position of the top-left corner of the region
    int         yMin;        //!< Y coordinate of the pixel position of the top-left corner of the region
    int         mv[3];       //!< motion vectors in 1/4 pixel units: mvx = mv[0], mvy = mv[1],
                             //!< and ref_frame = mv[2]
};

/* segment data structure */
struct ercSegment_t {
    short       startMBPos;
    short       endMBPos;
    char        fCorrupted;
};


/* "block" means an 8x8 pixel area */

/* Region modes */
#define REGMODE_INTER_COPY       0  //!< Copy region
#define REGMODE_INTER_PRED       1  //!< Inter region with motion vectors
#define REGMODE_INTRA            2  //!< Intra region
#define REGMODE_SPLITTED         3  //!< Any region mode higher than this indicates that the region
                                    //!< is splitted which means 8x8 block
#define REGMODE_INTER_COPY_8x8   4
#define REGMODE_INTER_PRED_8x8   5
#define REGMODE_INTRA_8x8        6


#define xPosYBlock(currYBlockNum, picSizeX) \
    ((currYBlockNum) % ((picSizeX) >> 3))

#define yPosYBlock(currYBlockNum, picSizeX) \
    ((currYBlockNum) / ((picSizeX) >> 3))

#define xPosMB(currMBNum, picSizeX) \
    ((currMBNum) % ((picSizeX) >> 4))

#define yPosMB(currMBNum, picSizeX) \
    ((currMBNum) / ((picSizeX) >> 4))

#define MBxy2YBlock(currXPos, currYPos, comp, picSizeX) \
    ((((currYPos) << 1) + ((comp) >> 1)) * ((picSizeX) >> 3) + ((currXPos) << 1) + ((comp) & 1))

#define MBNum2YBlock(currMBNum, comp, picSizeX) \
    MBxy2YBlock(xPosMB((currMBNum), (picSizeX)), yPosMB((currMBNum), (picSizeX)), (comp), (picSizeX))


/* If the average motion vector of the correctly received macroblocks is less than the
threshold, concealByCopy is used, otherwise concealByTrial is used. */
#define MVPERMB_THR 8

#define ERC_BLOCK_OK                3
#define ERC_BLOCK_CONCEALED         2
#define ERC_BLOCK_CORRUPTED         1
#define ERC_BLOCK_EMPTY             0


#define isSplitted(object_list, currMBNum) \
    ((object_list + ((currMBNum) << 2))->regionMode >= REGMODE_SPLITTED)

/* this can be used as isBlock(...,INTRA) or isBlock(...,INTER_COPY) */
#define isBlock(object_list, currMBNum, comp, regMode) \
    (isSplitted(object_list, currMBNum) ? \
     ((object_list + ((currMBNum) << 2) + (comp))->regionMode == REGMODE_##regMode##_8x8) : \
     ((object_list + ((currMBNum) << 2))->regionMode == REGMODE_##regMode))

/* this can be used as getParam(...,mv) or getParam(...,xMin) or getParam(...,yMin) */
#define getParam(object_list, currMBNum, comp, param) \
    (isSplitted(object_list, currMBNum) ? \
     ((object_list + ((currMBNum) << 2) + (comp))->param) : \
     ((object_list + ((currMBNum) << 2))->param))


static const int uv_div[2][4] = { //[x/y][yuv_format]
    {0, 1, 1, 0},
    {0, 1, 0, 0}
};

//! look up tables for FRExt_chroma support
static const uint8_t subblk_offset_x[3][8][4] = {
    {{  0,  4,  0,  4 },
     {  0,  4,  0,  4 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 }},
    {{  0,  4,  0,  4 },
     {  0,  4,  0,  4 },
     {  0,  4,  0,  4 },
     {  0,  4,  0,  4 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 }},
    {{  0,  4,  0,  4 },
     {  8, 12,  8, 12 },
     {  0,  4,  0,  4 },
     {  8, 12,  8, 12 },
     {  0,  4,  0,  4 },
     {  8, 12,  8, 12 },
     {  0,  4,  0,  4 },
     {  8, 12,  8, 12 }}
};

static const uint8_t subblk_offset_y[3][8][4] = {
    {{  0,  0,  4,  4 },
     {  0,  0,  4,  4 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 }},
    {{  0,  0,  4,  4 },
     {  8,  8, 12, 12 },
     {  0,  0,  4,  4 },
     {  8,  8, 12, 12 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 },
     {  0,  0,  0,  0 }},
    {{  0,  0,  4,  4 },
     {  0,  0,  4,  4 },
     {  8,  8, 12, 12 },
     {  8,  8, 12, 12 },
     {  0,  0,  4,  4 },
     {  0,  0,  4,  4 },
     {  8,  8, 12, 12 },
     {  8,  8, 12, 12 }}
};


ercVariables_t::ercVariables_t(int pic_sizex, int pic_sizey, bool flag)
{
    // the error concealment instance is allocated
    this->nOfMBs      = 0;
    this->segments    = nullptr;
    this->currSegment = 0;
    this->yCondition  = nullptr;
    this->uCondition  = nullptr;
    this->vCondition  = nullptr;
    this->concealment = flag;

    this->erc_object_list = new objectBuffer_t[(pic_sizex * pic_sizey) >> 6];

    this->erc_mvperMB = 0;
}

ercVariables_t::~ercVariables_t()
{
    if (this->erc_object_list)
        delete []this->erc_object_list;
    if (this->yCondition) {
        delete []this->segments;
        delete []this->yCondition;
        delete []this->uCondition;
        delete []this->vCondition;
    }
}

void ercVariables_t::reset(int nOfMBs, int numOfSegments)
{
    if (this->concealment) {
        // If frame size has been changed
        if (nOfMBs != this->nOfMBs && this->yCondition) {
            delete []this->segments;
            delete []this->yCondition;
            delete []this->uCondition;
            delete []this->vCondition;
            this->segments   = nullptr;
            this->yCondition = nullptr;
            this->uCondition = nullptr;
            this->vCondition = nullptr;
        }

        // If the structures are uninitialized (first frame, or frame size is changed)
        if (!this->yCondition) {
            this->nOfMBs      = nOfMBs;
            this->segments    = new ercSegment_t[numOfSegments];
            this->nOfSegments = numOfSegments;
            this->yCondition  = new char[4 * nOfMBs];
            this->uCondition  = new char[nOfMBs];
            this->vCondition  = new char[nOfMBs];
            memset(this->segments, 0, numOfSegments * sizeof(ercSegment_t));
        }

        // Reset tables and parameters
        memset(this->yCondition, 0, 4 * nOfMBs * sizeof(*this->yCondition));
        memset(this->uCondition, 0,     nOfMBs * sizeof(*this->uCondition));
        memset(this->vCondition, 0,     nOfMBs * sizeof(*this->vCondition));

        if (this->nOfSegments != numOfSegments) {
            delete []this->segments;
            this->segments = new ercSegment_t[numOfSegments];
            this->nOfSegments = numOfSegments;
        }

        ercSegment_t* segments = this->segments;
        for (int i = 0; i < this->nOfSegments; ++i) {
            segments->startMBPos = 0;
            segments->endMBPos = (short)(nOfMBs - 1);
            (segments++)->fCorrupted = 1; //! mark segments as corrupted
        }

        this->currSegment = 0;
        this->nOfCorruptedSegments = 0;
    }

    this->erc_mvperMB = 0;
}


void ercVariables_t::ercWriteMBMODEandMV(mb_t& mb, uint8_t slice_type, storable_picture* pic)
{
    objectBuffer_t* currRegion = this->erc_object_list + (mb.mbAddrX << 2);
    int mbx = xPosMB(mb.mbAddrX, pic->size_x);
    int mby = yPosMB(mb.mbAddrX, pic->size_x);

    if (slice_type != B_slice) { //non-B frame
        for (int i = 0; i < 4; ++i) {
            uint8_t mb_type     = mb.mb_type;
            uint8_t sub_mb_type = mb.SubMbType[i];

            objectBuffer_t* pRegion = currRegion + i;
            pRegion->regionMode = mb_type     == I_16x16 ? REGMODE_INTRA :
                                  sub_mb_type == I_4x4   ? REGMODE_INTRA_8x8  :
                                  sub_mb_type == 0       ? REGMODE_INTER_COPY :
                                  sub_mb_type == 1       ? REGMODE_INTER_PRED :
                                                           REGMODE_INTER_PRED_8x8;
            if (sub_mb_type == 0 || sub_mb_type == I_4x4) { // INTRA OR COPY
                pRegion->mv[0] = 0;
                pRegion->mv[1] = 0;
                pRegion->mv[2] = 0;
            } else {
                int ii = 4 * mbx + (i % 2) * 2;
                int jj = 4 * mby + (i / 2) * 2;
                if (sub_mb_type >= 5 && sub_mb_type <= 7) { // SMALL BLOCKS
                    pRegion->mv[0] = (pic->mv_info[jj    ][ii    ].mv[LIST_0].mv_x +
                                      pic->mv_info[jj    ][ii + 1].mv[LIST_0].mv_x +
                                      pic->mv_info[jj + 1][ii    ].mv[LIST_0].mv_x +
                                      pic->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_x + 2) / 4;
                    pRegion->mv[1] = (pic->mv_info[jj    ][ii    ].mv[LIST_0].mv_y +
                                      pic->mv_info[jj    ][ii + 1].mv[LIST_0].mv_y +
                                      pic->mv_info[jj + 1][ii    ].mv[LIST_0].mv_y +
                                      pic->mv_info[jj + 1][ii + 1].mv[LIST_0].mv_y + 2)/4;
                } else { // 16x16, 16x8, 8x16, 8x8
                    pRegion->mv[0] = pic->mv_info[jj][ii].mv[LIST_0].mv_x;
                    pRegion->mv[1] = pic->mv_info[jj][ii].mv[LIST_0].mv_y;
                }
                pRegion->mv[2] = pic->mv_info[jj][ii].ref_idx[LIST_0];

                this->erc_mvperMB += abs(pRegion->mv[0]) + abs(pRegion->mv[1]);
            }
        }
    } else { //B-frame
        for (int i = 0; i < 4; ++i) {
            uint8_t mb_type     = mb.mb_type;
            uint8_t sub_mb_type = mb.SubMbType[i];

            objectBuffer_t* pRegion = currRegion + i;
            pRegion->regionMode = mb_type     == I_16x16 ? REGMODE_INTRA :
                                  sub_mb_type == I_4x4   ? REGMODE_INTRA_8x8 :
                                                           REGMODE_INTER_PRED_8x8;
            if (mb_type == I_16x16 || sub_mb_type == I_4x4) { // INTRA
                pRegion->mv[0] = 0;
                pRegion->mv[1] = 0;
                pRegion->mv[2] = 0;
            } else {
                int ii = 4 * mbx + (i % 2) * 2;
                int jj = 4 * mby + (i / 2) * 2;
                int idx = (pic->mv_info[jj][ii].ref_idx[0] < 0) ? 1 : 0;
                pRegion->mv[0] = (pic->mv_info[jj    ][ii    ].mv[idx].mv_x + 
                                  pic->mv_info[jj    ][ii + 1].mv[idx].mv_x + 
                                  pic->mv_info[jj + 1][ii    ].mv[idx].mv_x + 
                                  pic->mv_info[jj + 1][ii + 1].mv[idx].mv_x + 2) / 4;
                pRegion->mv[1] = (pic->mv_info[jj    ][ii    ].mv[idx].mv_y + 
                                  pic->mv_info[jj    ][ii + 1].mv[idx].mv_y + 
                                  pic->mv_info[jj + 1][ii    ].mv[idx].mv_y + 
                                  pic->mv_info[jj + 1][ii + 1].mv[idx].mv_y + 2) / 4;
                pRegion->mv[2] = (pic->mv_info[jj][ii].ref_idx[idx]);

                this->erc_mvperMB += abs(pRegion->mv[0]) + abs(pRegion->mv[1]);
            }
        }
    }
}

void ercVariables_t::erc_picture(storable_picture* pic)
{
    slice_t& slice = *pic->slice_headers[0];
    shr_t& shr = slice.header;

    if (shr.MbaffFrameFlag)
        return;

    VideoParameters* p_Vid = slice.p_Vid;

    //! this is always true at the beginning of a picture
    int ercStartMB = 0;
    int ercSegment = 0;

    //! mark the start of the first segment
    int i = 0;
    this->ercStartSegment(i, ercSegment);
    //! generate the segments according to the macroblock map
    for (i = 1; i < shr.PicSizeInMbs; ++i) {
        if (p_Vid->mb_data[i].ei_flag != p_Vid->mb_data[i - 1].ei_flag) {
            this->ercStopSegment(i - 1, ercSegment); //! stop current segment

            //! mark current segment as lost or OK
            if (p_Vid->mb_data[i - 1].ei_flag)
                this->ercMarkCurrSegmentLost(pic->size_x);
            else
                this->ercMarkCurrSegmentOK(pic->size_x);

            ++ercSegment;  //! next segment
            this->ercStartSegment(i, ercSegment); //! start new segment
            ercStartMB = i;//! save start MB for this segment
        }
    }
    //! mark end of the last segment
    this->ercStopSegment(i - 1, ercSegment);
    if (p_Vid->mb_data[i - 1].ei_flag)
        this->ercMarkCurrSegmentLost(pic->size_x);
    else
        this->ercMarkCurrSegmentOK(pic->size_x);

    //! call the right error concealment function depending on the frame type.
    this->erc_mvperMB /= shr.PicSizeInMbs;

    if (shr.slice_type == I_slice || shr.slice_type == SI_slice) // I-frame
        this->ercConcealIntraFrame(pic);
    else
        this->ercConcealInterFrame(pic);
}


void ercVariables_t::ercStartSegment(int currMBNum, int segment)
{
    if (this->concealment) {
        this->currSegmentCorrupted = 0;
        this->segments[segment].fCorrupted = 0;
        this->segments[segment].startMBPos = (short)currMBNum;
    }
}

void ercVariables_t::ercStopSegment(int currMBNum, int segment)
{
    if (this->concealment) {
        this->segments[segment].endMBPos = (short)currMBNum;
        this->currSegment++;
    }
}

void ercVariables_t::ercMarkCurrSegmentLost(int picSizeX)
{
    int current_segment = this->currSegment - 1;

    if (this->concealment) {
        if (this->currSegmentCorrupted == 0) {
            this->nOfCorruptedSegments++;
            this->currSegmentCorrupted = 1;
        }

        for (int j = this->segments[current_segment].startMBPos; j <= this->segments[current_segment].endMBPos; ++j) {
            this->yCondition[MBNum2YBlock(j, 0, picSizeX)] = ERC_BLOCK_CORRUPTED;
            this->yCondition[MBNum2YBlock(j, 1, picSizeX)] = ERC_BLOCK_CORRUPTED;
            this->yCondition[MBNum2YBlock(j, 2, picSizeX)] = ERC_BLOCK_CORRUPTED;
            this->yCondition[MBNum2YBlock(j, 3, picSizeX)] = ERC_BLOCK_CORRUPTED;
            this->uCondition[j] = ERC_BLOCK_CORRUPTED;
            this->vCondition[j] = ERC_BLOCK_CORRUPTED;
        }
        this->segments[current_segment].fCorrupted = 1;
    }
}

void ercVariables_t::ercMarkCurrSegmentOK(int picSizeX)
{
    int current_segment = this->currSegment - 1;

    if (this->concealment) {
        // mark all the Blocks belonging to the segment as OK */
        for (int j = this->segments[current_segment].startMBPos; j <= this->segments[current_segment].endMBPos; ++j) {
            this->yCondition[MBNum2YBlock(j, 0, picSizeX)] = ERC_BLOCK_OK;
            this->yCondition[MBNum2YBlock(j, 1, picSizeX)] = ERC_BLOCK_OK;
            this->yCondition[MBNum2YBlock(j, 2, picSizeX)] = ERC_BLOCK_OK;
            this->yCondition[MBNum2YBlock(j, 3, picSizeX)] = ERC_BLOCK_OK;
            this->uCondition[j] = ERC_BLOCK_OK;
            this->vCondition[j] = ERC_BLOCK_OK;
        }
        this->segments[current_segment].fCorrupted = 0;
    }
}


int ercVariables_t::ercConcealIntraFrame(storable_picture* pic)
{
    // if concealment is on
    if (this->concealment) {
        // if there are segments to be concealed
        if (this->nOfCorruptedSegments) {
            this->concealBlocks(pic, 0);
            this->concealBlocks(pic, 1);
            this->concealBlocks(pic, 2);
        }
        return 1;
    }

    return 0;
}

int ercVariables_t::ercConcealInterFrame(storable_picture* pic)
{
    sps_t* sps = pic->sps;
    int picSizeX = pic->size_x;
    int picSizeY = pic->size_y;

    storable_picture* ref_pic = pic->slice_headers[0]->RefPicList[0][0];

    int predBlocks[8];

    /* if concealment is on */
    if (this->concealment) {
        /* if there are segments to be concealed */
        if (this->nOfCorruptedSegments) {
            px_t* predMB;
            if (sps->chroma_format_idc != CHROMA_FORMAT_400)
                predMB = new px_t[256 + sps->MbWidthC * sps->MbHeightC * 2];
            else
                predMB = new px_t[256];

            int lastRow    = (picSizeY >> 4);
            int lastColumn = (picSizeX >> 4);

            for (int columnInd = 0; columnInd < lastColumn; ++columnInd) {
                int column = (columnInd % 2) ? (lastColumn - columnInd / 2 - 1) : (columnInd / 2);

                for (int row = 0; row < lastRow; ++row) {
                    if (this->yCondition[MBxy2YBlock(column, row, 0, picSizeX)] <= ERC_BLOCK_CORRUPTED) {
                        // ERC_BLOCK_CORRUPTED (1) or ERC_BLOCK_EMPTY (0)
                        int firstCorruptedRow = row;
                        int lastCorruptedRow = -1;
                        /* find the last row which has corrupted blocks (in same continuous area) */
                        for (lastCorruptedRow = row + 1; lastCorruptedRow < lastRow; ++lastCorruptedRow) {
                            /* check blocks in the current column */
                            if (this->yCondition[MBxy2YBlock(column, lastCorruptedRow, 0, picSizeX)] > ERC_BLOCK_CORRUPTED) {
                                /* current one is already OK, so the last was the previous one */
                                --lastCorruptedRow;
                                break;
                            }
                        }
                        if (lastCorruptedRow >= lastRow ) {
                            row = lastRow;
                            lastCorruptedRow = lastRow - 1;

                            /* correct only from above */
                            for (int currRow = firstCorruptedRow; currRow < lastRow; ++currRow) {
                                this->ercCollect8PredBlocks(predBlocks, (currRow << 1), (column << 1),
                                                            this->yCondition, (lastRow << 1), (lastColumn << 1), 2, 0);

                                if (this->erc_mvperMB >= MVPERMB_THR)
                                    this->concealByTrial(pic, predMB, currRow * lastColumn + column, predBlocks);
                                else
                                    this->concealByCopy(pic, ref_pic, currRow * lastColumn + column);

                                this->ercMarkCurrMBConcealed(currRow * lastColumn + column, -1, picSizeX);
                            }
                        } else if (firstCorruptedRow == 0) {
                            row = lastCorruptedRow + 1;

                            /* correct only from below */
                            for (int currRow = lastCorruptedRow; currRow >= 0; --currRow) {
                                this->ercCollect8PredBlocks(predBlocks, (currRow << 1), (column << 1),
                                                            this->yCondition, (lastRow << 1), (lastColumn << 1), 2, 0);

                                if (this->erc_mvperMB >= MVPERMB_THR)
                                    this->concealByTrial(pic, predMB, currRow * lastColumn + column, predBlocks);
                                else
                                    this->concealByCopy(pic, ref_pic, currRow * lastColumn + column);

                                this->ercMarkCurrMBConcealed(currRow * lastColumn + column, -1, picSizeX);
                            }
                        } else {
                            /* correct bi-directionally */
                            row = lastCorruptedRow + 1;
                            int areaHeight = lastCorruptedRow - firstCorruptedRow + 1;

                            /*
                            *  Conceal the corrupted area switching between the up and the bottom rows
                            */
                            for (int i = 0; i < areaHeight; ++i) {
                                int currRow = i % 2 ? lastCorruptedRow : firstCorruptedRow;
                                if (i % 2)
                                    --lastCorruptedRow;
                                else
                                    ++firstCorruptedRow;

                                this->ercCollect8PredBlocks(predBlocks, (currRow << 1), (column << 1),
                                                            this->yCondition, (lastRow << 1), (lastColumn << 1), 2, 0);

                                if (this->erc_mvperMB >= MVPERMB_THR)
                                    this->concealByTrial(pic, predMB, currRow * lastColumn + column, predBlocks);
                                else
                                    this->concealByCopy(pic, ref_pic, currRow * lastColumn + column);

                                this->ercMarkCurrMBConcealed(currRow * lastColumn + column, -1, picSizeX);
                            }
                        }
                    }
                }
            }

            delete []predMB;
        }
        return 1;
    }

    return 0;
}


/*!
 ************************************************************************
 * \brief
 *      This function checks the neighbors of a mb_t for usability in
 *      concealment. First the OK macroblocks are marked, and if there is not
 *      enough of them, then the CONCEALED ones as well.
 *      A "1" in the the output array means reliable, a "0" non reliable MB.
 *      The block order in "predBlocks":
 *              1 4 0
 *              5 x 7
 *              2 6 3
 *      i.e., corners first.
 ************************************************************************
 */
int ercVariables_t::ercCollect8PredBlocks(int predBlocks[], int currRow, int currColumn, char* condition,
                                          int maxRow, int maxColumn, int step, uint8_t fNoCornerNeigh)
{
    int srcCounter  = 0;
    int srcCountMin = (fNoCornerNeigh ? 2 : 4);
    int threshold   = ERC_BLOCK_OK;

    memset(predBlocks, 0, 8 * sizeof(int));

    // collect the reliable neighboring blocks
    do {
        srcCounter = 0;
        // top
        if (currRow > 0 && condition[(currRow - 1) * maxColumn + currColumn] >= threshold) {
            //ERC_BLOCK_OK (3) or ERC_BLOCK_CONCEALED (2)
            predBlocks[4] = condition[(currRow - 1) * maxColumn + currColumn];
            srcCounter++;
        }
        // bottom
        if (currRow < (maxRow - step) && condition[(currRow + step) * maxColumn + currColumn] >= threshold) {
            predBlocks[6] = condition[(currRow + step) * maxColumn + currColumn];
            srcCounter++;
        }

        if (currColumn > 0) {
            // left
            if (condition[currRow*maxColumn + currColumn - 1] >= threshold) {
                predBlocks[5] = condition[currRow * maxColumn + currColumn - 1];
                srcCounter++;
            }

            if (!fNoCornerNeigh) {
                // top-left
                if (currRow > 0 && condition[(currRow - 1) * maxColumn + currColumn - 1] >= threshold) {
                    predBlocks[1] = condition[(currRow - 1) * maxColumn + currColumn - 1];
                    srcCounter++;
                }
                // bottom-left
                if (currRow < (maxRow - step) && condition[(currRow + step) * maxColumn + currColumn - 1] >= threshold) {
                    predBlocks[2] = condition[(currRow + step) * maxColumn + currColumn - 1];
                    srcCounter++;
                }
            }
        }

        if (currColumn < (maxColumn - step)) {
            // right
            if (condition[currRow * maxColumn + currColumn + step] >= threshold) {
                predBlocks[7] = condition[currRow * maxColumn + currColumn + step];
                srcCounter++;
            }

            if (!fNoCornerNeigh) {
                // top-right
                if (currRow > 0 && condition[(currRow - 1) * maxColumn + currColumn + step] >= threshold) {
                    predBlocks[0] = condition[(currRow - 1) * maxColumn + currColumn + step];
                    srcCounter++;
                }
                // bottom-right
                if (currRow < (maxRow - step) && condition[(currRow + step) * maxColumn + currColumn + step] >= threshold) {
                    predBlocks[3] = condition[(currRow + step) * maxColumn + currColumn + step];
                    srcCounter++;
                }
            }
        }
        // prepare for the next round
        threshold--;
        if (threshold < ERC_BLOCK_CONCEALED)
            break;
    } while (srcCounter < srcCountMin);

    return srcCounter;
}


int ercVariables_t::ercCollectColumnBlocks(int predBlocks[], int currRow, int currColumn, char* condition, int maxRow, int maxColumn, int step)
{
    int srcCounter = 0, threshold = ERC_BLOCK_CORRUPTED;

    memset(predBlocks, 0, 8 * sizeof(int));

    // in this case, row > 0 and row < 17
    if (condition[(currRow - 1) * maxColumn + currColumn] > threshold) {
        predBlocks[4] = 1;
        srcCounter++;
    }
    if (condition[(currRow + step) * maxColumn + currColumn] > threshold) {
        predBlocks[6] = 1;
        srcCounter++;
    }

    return srcCounter;
}

void ercVariables_t::pixMeanInterpolateBlock(px_t* src[], px_t* block, int blockSize, int frameWidth, int BitDepth)
{
    int bmax = blockSize - 1;
    int k = 0;
    for (int row = 0; row < blockSize; ++row) {
        for (int column = 0; column < blockSize; ++column) {
            int tmp = 0;
            int srcCounter = 0;
            int weight = 0;
            // above
            if (src[4]) {
                weight = blockSize - row;
                tmp += weight * src[4][bmax * frameWidth + column];
                srcCounter += weight;
            }
            // left
            if (src[5]) {
                weight = blockSize - column;
                tmp += weight * src[5][row * frameWidth + bmax];
                srcCounter += weight;
            }
            // below
            if (src[6]) {
                weight = row + 1;
                tmp += weight * src[6][column];
                srcCounter += weight;
            }
            // right
            if (src[7]) {
                weight = column + 1;
                tmp += weight * src[7][row * frameWidth];
                srcCounter += weight;
            }

            if (srcCounter > 0)
                block[k + column] = (px_t)(tmp / srcCounter);
            else
                block[k + column] = (px_t)(1 << (BitDepth - 1));
        }
        k += frameWidth;
    }
}

void ercVariables_t::ercPixConcealIMB(storable_picture* pic, int comp, int row, int column, int predBlocks[])
{
    sps_t* sps = pic->sps;
    px_t* currFrame     = comp == 0 ? &pic->imgY[0][0] : comp == 1 ? &pic->imgUV[0][0][0] : &pic->imgUV[1][0][0];
    int frameWidth      = comp == 0 ? pic->size_x : pic->size_x >> 1;
    int mbWidthInBlocks = comp == 0 ? 2 : 1;
    int BitDepth        = comp == 0 ? sps->BitDepthY : sps->BitDepthC;

    px_t* src[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
    px_t* currBlock = NULL;

    // collect the reliable neighboring blocks
    if (predBlocks[0])
        src[0] = currFrame + (row - mbWidthInBlocks) * frameWidth * 8 + (column + mbWidthInBlocks) * 8;
    if (predBlocks[1])
        src[1] = currFrame + (row - mbWidthInBlocks) * frameWidth * 8 + (column - mbWidthInBlocks) * 8;
    if (predBlocks[2])
        src[2] = currFrame + (row + mbWidthInBlocks) * frameWidth * 8 + (column - mbWidthInBlocks) * 8;
    if (predBlocks[3])
        src[3] = currFrame + (row + mbWidthInBlocks) * frameWidth * 8 + (column + mbWidthInBlocks) * 8;
    if (predBlocks[4])
        src[4] = currFrame + (row - mbWidthInBlocks) * frameWidth * 8 + column * 8;
    if (predBlocks[5])
        src[5] = currFrame + row * frameWidth * 8 + (column - mbWidthInBlocks) * 8;
    if (predBlocks[6])
        src[6] = currFrame + (row + mbWidthInBlocks) * frameWidth * 8 + column * 8;
    if (predBlocks[7])
        src[7] = currFrame + row * frameWidth * 8 + (column + mbWidthInBlocks) * 8;

    currBlock = currFrame + row * frameWidth * 8 + column * 8;
    this->pixMeanInterpolateBlock(src, currBlock, mbWidthInBlocks * 8, frameWidth, BitDepth);
}

void ercVariables_t::concealBlocks(storable_picture* pic, int comp)
{
    int lastColumn  = comp == 0 ? pic->size_x >> 3 : pic->size_x >> 4;
    int lastRow     = comp == 0 ? pic->size_y >> 3 : pic->size_y >> 4;
    char* condition = comp == 0 ? this->yCondition : comp == 1 ? this->uCondition : this->vCondition;

    int predBlocks[8];

    // in the Y component do the concealment MB-wise (not block-wise):
    // this is useful if only whole MBs can be damaged or lost
    int step = comp == 0 ? 2 : 1;

    for (int column = 0; column < lastColumn; column += step) {
        for (int row = 0; row < lastRow; row += step) {
            if (condition[row * lastColumn + column] <= ERC_BLOCK_CORRUPTED) {
                int firstCorruptedRow = row;
                int lastCorruptedRow = -1;
                // find the last row which has corrupted blocks (in same continuous area)
                for (lastCorruptedRow = row + step; lastCorruptedRow < lastRow; lastCorruptedRow += step) {
                    // check blocks in the current column
                    if (condition[lastCorruptedRow * lastColumn + column] > ERC_BLOCK_CORRUPTED) {
                        // current one is already OK, so the last was the previous one
                        lastCorruptedRow -= step;
                        break;
                    }
                }

                if (lastCorruptedRow >= lastRow) {
                    // correct only from above
                    row = lastRow;
                    lastCorruptedRow = lastRow - step;

                    for (int currRow = firstCorruptedRow; currRow < lastRow; currRow += step) {
                        this->ercCollect8PredBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step, 1);
                        this->ercPixConcealIMB(pic, comp, currRow, column, predBlocks);

                        if (comp == 0) {
                            condition[currRow * lastColumn + column                 ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column              + 1] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn    ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn + 1] = ERC_BLOCK_CONCEALED;
                        } else
                            condition[currRow * lastColumn + column] = ERC_BLOCK_CONCEALED;
                    }
                } else if (firstCorruptedRow == 0) {
                    // correct only from below
                    row = lastCorruptedRow + step;

                    for (int currRow = lastCorruptedRow; currRow >= 0; currRow -= step) {
                        this->ercCollect8PredBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step, 1);
                        this->ercPixConcealIMB(pic, comp, currRow, column, predBlocks);

                        condition[currRow * lastColumn + column] = ERC_BLOCK_CONCEALED;
                        if (comp == 0) {
                            condition[currRow * lastColumn + column              + 1] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn    ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn + 1] = ERC_BLOCK_CONCEALED;
                        }
                    }
                } else {
                    // correct bi-directionally

                    row = lastCorruptedRow + step;
                    int areaHeight = lastCorruptedRow - firstCorruptedRow+step;

                    // Conceal the corrupted area switching between the up and the bottom rows
                    for (int i = 0; i < areaHeight; i += step) {
                        int currRow = i % 2 ? lastCorruptedRow : firstCorruptedRow;
                        if (i % 2)
                            lastCorruptedRow -= step;
                        else
                            firstCorruptedRow += step;

                        int smoothColumn = 0;
                        if (smoothColumn > 0)
                            this->ercCollectColumnBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step);
                        else
                            this->ercCollect8PredBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step, 1);
                        this->ercPixConcealIMB(pic, comp, currRow, column, predBlocks);

                        condition[currRow * lastColumn + column] = ERC_BLOCK_CONCEALED;
                        if (comp == 0) {
                            condition[currRow * lastColumn + column              + 1] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn    ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn + 1] = ERC_BLOCK_CONCEALED;
                        }
                    }
                }
            }
        }
    }
}


void ercVariables_t::buildPredRegionYUV(storable_picture* pic, int* mv, int x, int y, px_t* predMB)
{
    sps_t* sps = pic->sps;
    int yuv = sps->chroma_format_idc - 1;
    int ref_frame = max(mv[2], 0); // !!KS: quick fix, we sometimes seem to get negative ref_pic here, so restrict to zero and above
    int mb_nr = (y / 16) * (sps->PicWidthInMbs) + (x / 16);

    // This should be allocated only once. 
    px_t tmp_block[16][16];

    /* Update coordinates of the current concealed macroblock */
    VideoParameters* p_Vid = pic->slice_headers[0]->p_Vid;
    mb_t& mb = p_Vid->mb_data[mb_nr];   // intialization code deleted, see below, StW  
    mb.mb.x = (short)(x / 16);
    mb.mb.y = (short)(y / 16);
    slice_t& slice = *mb.p_Slice;
    auto ref_pic = slice.RefPicList[0][ref_frame];
  
    int mv_mul = 4;

    int num_uv_blocks;
    if (sps->chroma_format_idc != CHROMA_FORMAT_400)
        num_uv_blocks = (((1 << sps->chroma_format_idc) & (~0x1)) >> 1);
    else
        num_uv_blocks = 0;

    // luma *******************************************************
    for (int j = 0; j < 16/4; ++j) {
        int joff = j * 4;
        int j4 = mb.mb.y * 4 + j;
        for (int i = 0; i < 16/4; i++) {
            int ioff = i * 4;
            int i4 = mb.mb.x * 4 + i;

            int vec1_x = i4 * 4 * mv_mul + mv[0];
            int vec1_y = j4 * 4 * mv_mul + mv[1];

            slice.decoder.get_block_luma(ref_pic,
                vec1_x, vec1_y, 4, 4, tmp_block, pic->iLumaStride, pic->size_x - 1,
                mb.mb_field_decoding_flag ? (pic->size_y >> 1) - 1 : pic->size_y - 1,
                PLANE_Y, &mb);

            for (int ii = 0; ii < 4; ++ii) {
                for (int jj = 0; jj < 16/4; ++jj)
                    slice.mb_pred[PLANE_Y][jj + joff][ii + ioff] = tmp_block[jj][ii];
            }
        }
    }

    px_t* pMB = predMB;
    for (int j = 0; j < 16; ++j) {
        for (int i = 0; i < 16; ++i)
            pMB[j * 16 + i] = slice.mb_pred[PLANE_Y][j][i];
    }
    pMB += 256;

    if (sps->chroma_format_idc != CHROMA_FORMAT_400) {
        // chroma *******************************************************
        int f1_x = 64 / sps->MbWidthC;
        int f2_x = f1_x - 1;

        int f1_y = 64 / sps->MbHeightC;
        int f2_y = f1_y - 1;

        int f3 = f1_x * f1_y;
        int f4 = f3 >> 1;

        for (int uv = 0; uv < 2; ++uv) {
            for (int b8 = 0; b8 < num_uv_blocks; ++b8) {
                for (int b4 = 0; b4 < 4; ++b4) {
                    int joff = subblk_offset_y[yuv][b8][b4];
                    int j4   = mb.mb.y * sps->MbHeightC + joff;
                    int ioff = subblk_offset_x[yuv][b8][b4];
                    int i4   = mb.mb.x * sps->MbWidthC + ioff;

                    for (int jj = 0; jj < 4; ++jj) {
                        for (int ii = 0; ii < 4; ++ii) {
                            int i1 = (i4 + ii) * f1_x + mv[0];
                            int j1 = (j4 + jj) * f1_y + mv[1];

                            int ii0 = clip3(0, pic->size_x_cr - 1, i1 / f1_x);
                            int jj0 = clip3(0, pic->size_y_cr - 1, j1 / f1_y);
                            int ii1 = clip3(0, pic->size_x_cr - 1, (i1 + f2_x) / f1_x);
                            int jj1 = clip3(0, pic->size_y_cr - 1, (j1 + f2_y) / f1_y);

                            int if1 = (i1 & f2_x);
                            int jf1 = (j1 & f2_y);
                            int if0 = (f1_x - if1);
                            int jf0 = (f1_y - jf1);

                            slice.mb_pred[uv + 1][jj + joff][ii + ioff] = (px_t) 
                                ((if0 * jf0 * ref_pic->imgUV[uv][jj0][ii0] +
                                  if1 * jf0 * ref_pic->imgUV[uv][jj0][ii1] +
                                  if0 * jf1 * ref_pic->imgUV[uv][jj1][ii0] +
                                  if1 * jf1 * ref_pic->imgUV[uv][jj1][ii1] + f4) / f3);
                        }
                    }
                }
            }

            for (int j = 0; j < 8; ++j) {
                for (int i = 0; i < 8; ++i)
                    pMB[j * 8 + i] = slice.mb_pred[uv + 1][j][i];
            }
            pMB += 64;
        }
    }
}

int ercVariables_t::edgeDistortion(storable_picture* pic, int predBlocks[], int currYBlockNum, px_t* predMB, int regionSize)
{
    px_t* recY = &pic->imgY[0][0];
    int picSizeX = pic->size_x;

    int threshold = ERC_BLOCK_OK;
    px_t* currBlock = recY + (yPosYBlock(currYBlockNum,picSizeX) << 3) * picSizeX +
                             (xPosYBlock(currYBlockNum,picSizeX) << 3);
    int distortion, numOfPredBlocks;

    do {
        distortion = 0;
        numOfPredBlocks = 0;

        // loop the 4 neighbors
        for (int j = 4; j < 8; ++j) {
            /* if reliable, count boundary pixel difference */
            if (predBlocks[j] >= threshold) {
                px_t* neighbor = NULL;
                int currBlockOffset = 0;
                switch (j) {
                case 4:
                    neighbor = currBlock - picSizeX;
                    for (int i = 0; i < regionSize; ++i)
                        distortion += abs((int)(predMB[i] - neighbor[i]));
                    break;
                case 5:
                    neighbor = currBlock - 1;
                    for (int i = 0; i < regionSize; ++i)
                        distortion += abs((int)(predMB[i * 16] - neighbor[i * picSizeX]));
                    break;
                case 6:
                    neighbor = currBlock + regionSize*picSizeX;
                    currBlockOffset = (regionSize - 1) * 16;
                    for (int i = 0; i < regionSize; ++i)
                        distortion += abs((int)(predMB[i + currBlockOffset] - neighbor[i]));
                    break;
                case 7:
                    neighbor = currBlock + regionSize;
                    currBlockOffset = regionSize-1;
                    for (int i = 0; i < regionSize; ++i)
                        distortion += abs((int)(predMB[i * 16 + currBlockOffset] - neighbor[i * picSizeX]));
                    break;
                }
                numOfPredBlocks++;
            }
        }

        threshold--;
        if (threshold < ERC_BLOCK_CONCEALED)
            break;
    } while (numOfPredBlocks == 0);

    if (numOfPredBlocks == 0)
        return 0;
        // assert (numOfPredBlocks != 0); !!!KS hmm, trying to continue...
    return (distortion / numOfPredBlocks);
}

void ercVariables_t::copyPredMB(storable_picture* pic, int currYBlockNum, px_t* predMB, int regionSize)
{
    int picSizeX = pic->size_x;
    sps_t* sps = pic->sps;

    int uv_x = uv_div[0][sps->chroma_format_idc];
    int uv_y = uv_div[1][sps->chroma_format_idc];

    int xmin = (xPosYBlock(currYBlockNum, picSizeX) << 3);
    int ymin = (yPosYBlock(currYBlockNum, picSizeX) << 3);
    int xmax = xmin + regionSize - 1;
    int ymax = ymin + regionSize - 1;

    for (int j = ymin; j <= ymax; ++j) {
        for (int k = xmin; k <= xmax; ++k) {
            int locationTmp = (j - ymin) * 16 + (k - xmin);
            pic->imgY[j][k] = predMB[locationTmp];
        }
    }

    if (sps->chroma_format_idc != CHROMA_FORMAT_400) {
        for (int j = (ymin >> uv_y); j <= (ymax >> uv_y); ++j) {
            for (int k = (xmin >> uv_x); k <= (xmax >> uv_x); ++k) {
                int locationTmp = (j - (ymin >> uv_y)) * sps->MbWidthC + (k - (xmin >> 1)) + 256;
                pic->imgUV[0][j][k] = predMB[locationTmp];

                locationTmp += 64;

                pic->imgUV[1][j][k] = predMB[locationTmp];
            }
        }
    }
}

int ercVariables_t::concealByTrial(storable_picture* pic, px_t* predMB, int currMBNum, int predBlocks[])
{
    int picSizeX = pic->size_x;

    int compLeft = 1;
    int threshold = ERC_BLOCK_OK;
    int minDist = 0;
    int fInterNeighborExists = 0;
    int fZeroMotionChecked = 0;
    int mvBest[3] = {0, 0, 0}, mvPred[3] = {0, 0, 0}, *mvptr;

    int numMBPerLine = (picSizeX >> 4);

    int comp = 0;
    int regionSize = 16;

    do { /* 4 blocks loop */
        objectBuffer_t* currRegion = this->erc_object_list + (currMBNum << 2) + comp;

        /* set the position of the region to be concealed */
        currRegion->xMin = (xPosYBlock(MBNum2YBlock(currMBNum, comp, picSizeX), picSizeX) << 3);
        currRegion->yMin = (yPosYBlock(MBNum2YBlock(currMBNum, comp, picSizeX), picSizeX) << 3);

        do { /* reliability loop */
            int numIntraNeighbours = 0;
            minDist = 0;
            fInterNeighborExists = 0;
            fZeroMotionChecked = 0;

            /* loop the 4 neighbours */
            for (int i = 4; i < 8; ++i) {
                /* if reliable, try it */
                if (predBlocks[i] >= threshold) {
                    int predMBNum = 0;
                    int compSplit1 = 0;
                    int compSplit2 = 0;
                    switch (i) {
                    case 4:
                        predMBNum = currMBNum - numMBPerLine;
                        compSplit1 = 2;
                        compSplit2 = 3;
                        break;
                    case 5:
                        predMBNum = currMBNum - 1;
                        compSplit1 = 1;
                        compSplit2 = 3;
                        break;
                    case 6:
                        predMBNum = currMBNum + numMBPerLine;
                        compSplit1 = 0;
                        compSplit2 = 1;
                        break;
                    case 7:
                        predMBNum = currMBNum + 1;
                        compSplit1 = 0;
                        compSplit2 = 2;
                        break;
                    }

                    /* try the concealment with the Motion Info of the current neighbour
                    only try if the neighbour is not Intra */
                    if (isBlock(this->erc_object_list, predMBNum, compSplit1, INTRA) ||
                        isBlock(this->erc_object_list, predMBNum, compSplit2, INTRA))
                        numIntraNeighbours++;
                    else {
                        /* if neighbour MB is splitted, try both neighbour blocks */
                        for (int predSplitted = isSplitted(this->erc_object_list, predMBNum), compPred = compSplit1;
                             predSplitted >= 0;
                             compPred = compSplit2, predSplitted -= (compSplit1 == compSplit2 ? 2 : 1)) {
                            /* if Zero Motion Block, do the copying. This option is tried only once */
                            if (isBlock(this->erc_object_list, predMBNum, compPred, INTER_COPY)) {
                                if (fZeroMotionChecked)
                                    continue;
                                else {
                                    fZeroMotionChecked = 1;
                                    mvPred[0] = mvPred[1] = mvPred[2] = 0;

                                    this->buildPredRegionYUV(pic, mvPred, currRegion->xMin, currRegion->yMin, predMB);
                                }
                            } else if (isBlock(this->erc_object_list, predMBNum, compPred, INTRA))
                                /* build motion using the neighbour's Motion Parameters */
                                continue;
                            else {
                                mvptr = getParam(this->erc_object_list, predMBNum, compPred, mv);

                                mvPred[0] = mvptr[0];
                                mvPred[1] = mvptr[1];
                                mvPred[2] = mvptr[2];

                                this->buildPredRegionYUV(pic, mvPred, currRegion->xMin, currRegion->yMin, predMB);
                            }

                            /* measure absolute boundary pixel difference */
                            int currDist = this->edgeDistortion(pic, predBlocks,
                                MBNum2YBlock(currMBNum, comp, picSizeX), predMB, regionSize);

                            /* if so far best -> store the pixels as the best concealment */
                            if (currDist < minDist || !fInterNeighborExists) {
                                minDist = currDist;
                                for (int k = 0; k < 3; ++k)
                                    mvBest[k] = mvPred[k];

                                currRegion->regionMode = isBlock(this->erc_object_list, predMBNum, compPred, INTER_COPY) ?
                                    (regionSize == 16 ? REGMODE_INTER_COPY : REGMODE_INTER_COPY_8x8) :
                                    (regionSize == 16 ? REGMODE_INTER_PRED : REGMODE_INTER_PRED_8x8);

                                this->copyPredMB(pic, MBNum2YBlock(currMBNum, comp, picSizeX), predMB, regionSize);
                            }

                            fInterNeighborExists = 1;
                        }
                    }
                }
            }

            threshold--;

        } while ((threshold >= ERC_BLOCK_CONCEALED) && (fInterNeighborExists == 0));

        /* always try zero motion */
        if (!fZeroMotionChecked) {
            mvPred[0] = mvPred[1] = mvPred[2] = 0;

            this->buildPredRegionYUV(pic, mvPred, currRegion->xMin, currRegion->yMin, predMB);

            int currDist = this->edgeDistortion(pic, predBlocks,
                MBNum2YBlock(currMBNum, comp, picSizeX), predMB, regionSize);

            if (currDist < minDist || !fInterNeighborExists) {
                minDist = currDist;
                for (int k = 0; k < 3; ++k)
                    mvBest[k] = mvPred[k];

                currRegion->regionMode = (regionSize == 16 ? REGMODE_INTER_COPY : REGMODE_INTER_COPY_8x8);

                this->copyPredMB(pic, MBNum2YBlock(currMBNum, comp, picSizeX), predMB, regionSize);
            }
        }

        for (int i = 0; i < 3; ++i)
            currRegion->mv[i] = mvBest[i];

        int order = 1;
        this->yCondition[MBNum2YBlock(currMBNum, comp, picSizeX)] = ERC_BLOCK_CONCEALED;
        comp = (comp + order + 4) % 4;
        --compLeft;

    } while (compLeft);

    return 0;
}

void ercVariables_t::copyBetweenFrames(storable_picture* dec_pic, storable_picture* ref_pic, int currYBlockNum, int regionSize)
{
    int picSizeX = dec_pic->size_x;
    sps_t* sps = dec_pic->sps;
    px_t* yptr = &dec_pic->imgY[0][0];
    px_t* uptr = &dec_pic->imgUV[0][0][0];
    px_t* vptr = &dec_pic->imgUV[1][0][0];

    /* set the position of the region to be copied */
    int xmin = (xPosYBlock(currYBlockNum, picSizeX) << 3);
    int ymin = (yPosYBlock(currYBlockNum, picSizeX) << 3);

    for (int j = ymin; j < ymin + regionSize; ++j) {
        for (int k = xmin; k < xmin + regionSize; ++k) {
            int location   = j * picSizeX + k;
            yptr[location] = ref_pic->imgY[j][k];
        }
    }

    for (int j = ymin >> uv_div[1][sps->chroma_format_idc];
         j < (ymin + regionSize) >> uv_div[1][sps->chroma_format_idc]; ++j) {
        for (int k = xmin >> uv_div[0][sps->chroma_format_idc];
             k < (xmin + regionSize) >> uv_div[0][sps->chroma_format_idc]; ++k) {
            int location   = ((j * picSizeX) >> uv_div[0][sps->chroma_format_idc]) + k;
            uptr[location] = ref_pic->imgUV[0][j][k];
            vptr[location] = ref_pic->imgUV[1][j][k];
        }
    }
}

int ercVariables_t::concealByCopy(storable_picture* dec_pic, storable_picture* ref_pic, int currMBNum)
{
    int picSizeX = dec_pic->size_x;
    objectBuffer_t* currRegion = this->erc_object_list + (currMBNum << 2);

    currRegion->regionMode = REGMODE_INTER_COPY;
    currRegion->xMin       = (xPosMB(currMBNum, picSizeX) << 4);
    currRegion->yMin       = (yPosMB(currMBNum, picSizeX) << 4);

    this->copyBetweenFrames(dec_pic, ref_pic, MBNum2YBlock(currMBNum, 0, picSizeX), 16);

    return 0;
}

void ercVariables_t::ercMarkCurrMBConcealed(int currMBNum, int comp, int picSizeX)
{
    if (this->concealment) {
        int setAll = 0;
        if (comp < 0) {
            setAll = 1;
            comp = 0;
        }

        switch (comp) {
        case 0:
            this->yCondition[MBNum2YBlock(currMBNum, 0, picSizeX)] = ERC_BLOCK_CONCEALED;
            this->yCondition[MBNum2YBlock(currMBNum, 1, picSizeX)] = ERC_BLOCK_CONCEALED;
            this->yCondition[MBNum2YBlock(currMBNum, 2, picSizeX)] = ERC_BLOCK_CONCEALED;
            this->yCondition[MBNum2YBlock(currMBNum, 3, picSizeX)] = ERC_BLOCK_CONCEALED;
            if (!setAll)
                break;
        case 1:
            this->uCondition[currMBNum] = ERC_BLOCK_CONCEALED;
            if (!setAll)
                break;
        case 2:
            this->vCondition[currMBNum] = ERC_BLOCK_CONCEALED;
        }
    }
}
