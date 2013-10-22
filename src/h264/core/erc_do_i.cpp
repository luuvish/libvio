#include "global.h"
#include "erc_do.h"


static void pixMeanInterpolateBlock(VideoParameters* p_Vid, px_t* src[], px_t* block, int blockSize, int frameWidth)
{
    sps_t* sps = p_Vid->active_sps;

    int bmax = blockSize - 1;
    int k = 0;
    for (int row = 0; row < blockSize; ++row) {
        for (int column = 0; column < blockSize; ++column) {
            int tmp = 0;
            int srcCounter = 0;
            int weight = 0;
            // above
            if (src[4]) {
                weight = blockSize-row;
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
                block[k + column] = (px_t)(blockSize == 8 ? (1 << (sps->BitDepthC - 1)) : (1 << (sps->BitDepthY - 1)));
        }
        k += frameWidth;
    }
}

static void ercPixConcealIMB(VideoParameters* p_Vid, px_t* currFrame, int row, int column, int predBlocks[], int frameWidth, int mbWidthInBlocks)
{
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
    pixMeanInterpolateBlock(p_Vid, src, currBlock, mbWidthInBlocks * 8, frameWidth);
}

static int ercCollectColumnBlocks(int predBlocks[], int currRow, int currColumn, char* condition, int maxRow, int maxColumn, int step)
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
int ercCollect8PredBlocks(int predBlocks[], int currRow, int currColumn, char* condition,
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


static void concealBlocks(VideoParameters* p_Vid, int lastColumn, int lastRow, int comp, frame* recfr, int picSizeX, char* condition)
{
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
                        ercCollect8PredBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step, 1);

                        switch (comp) {
                        case 0:
                            ercPixConcealIMB(p_Vid, recfr->yptr, currRow, column, predBlocks, picSizeX, 2);
                            break;
                        case 1:
                            ercPixConcealIMB(p_Vid, recfr->uptr, currRow, column, predBlocks, (picSizeX >> 1), 1);
                            break;
                        case 2 :
                            ercPixConcealIMB(p_Vid, recfr->vptr, currRow, column, predBlocks, (picSizeX >> 1), 1);
                            break;
                        }

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
                        ercCollect8PredBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step, 1);

                        switch (comp) {
                        case 0:
                            ercPixConcealIMB(p_Vid, recfr->yptr, currRow, column, predBlocks, picSizeX, 2);
                            break;
                        case 1:
                            ercPixConcealIMB(p_Vid, recfr->uptr, currRow, column, predBlocks, (picSizeX >> 1), 1);
                            break;
                        case 2:
                            ercPixConcealIMB(p_Vid, recfr->vptr, currRow, column, predBlocks, (picSizeX >> 1), 1);
                            break;
                        }

                        if (comp == 0) {
                            condition[currRow * lastColumn + column                 ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column              + 1] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn    ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn + 1] = ERC_BLOCK_CONCEALED;
                        } else
                            condition[currRow * lastColumn + column] = ERC_BLOCK_CONCEALED;
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
                            ercCollectColumnBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step);
                        else
                            ercCollect8PredBlocks(predBlocks, currRow, column, condition, lastRow, lastColumn, step, 1);

                        switch (comp) {
                        case 0:
                            ercPixConcealIMB(p_Vid, recfr->yptr, currRow, column, predBlocks, picSizeX, 2);
                            break;
                        case 1:
                            ercPixConcealIMB(p_Vid, recfr->uptr, currRow, column, predBlocks, (picSizeX >> 1), 1);
                            break;
                        case 2:
                            ercPixConcealIMB(p_Vid, recfr->vptr, currRow, column, predBlocks, (picSizeX >> 1), 1);
                            break;
                        }

                        if (comp == 0) {
                            condition[currRow * lastColumn + column                 ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column              + 1] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn    ] = ERC_BLOCK_CONCEALED;
                            condition[currRow * lastColumn + column + lastColumn + 1] = ERC_BLOCK_CONCEALED;
                        } else
                            condition[currRow * lastColumn + column] = ERC_BLOCK_CONCEALED;
                    }
                }

                lastCorruptedRow = -1;
                firstCorruptedRow = -1;
            }
        }
    }
}

int ercConcealIntraFrame(VideoParameters* p_Vid, frame* recfr, int picSizeX, int picSizeY, ercVariables_t* errorVar)
{
    // if concealment is on
    if (errorVar && errorVar->concealment) {
        // if there are segments to be concealed
        if (errorVar->nOfCorruptedSegments) {
            concealBlocks(p_Vid, (picSizeX >> 3), (picSizeY >> 3), 0, recfr, picSizeX, errorVar->yCondition);
            concealBlocks(p_Vid, (picSizeX >> 4), (picSizeY >> 4), 1, recfr, picSizeX, errorVar->uCondition);
            concealBlocks(p_Vid, (picSizeX >> 4), (picSizeY >> 4), 2, recfr, picSizeX, errorVar->vCondition);
        }
        return 1;
    }

    return 0;
}
