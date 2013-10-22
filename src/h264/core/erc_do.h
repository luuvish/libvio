#ifndef _ERC_DO_H_
#define _ERC_DO_H_


#include "erc_api.h"

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

/* used to determine the size of the allocated memory for a temporal Region (MB) */
#define DEF_REGION_SIZE 384  /* 8*8*6 */

#define ERC_BLOCK_OK                3
#define ERC_BLOCK_CONCEALED         2
#define ERC_BLOCK_CORRUPTED         1
#define ERC_BLOCK_EMPTY             0


//! YUV pixel domain image arrays for a video frame
struct frame {
    VideoParameters* p_Vid;
    px_t*       yptr;
    px_t*       uptr;
    px_t*       vptr;
};

int ercConcealIntraFrame(VideoParameters *p_Vid, frame *recfr, int picSizeX, int picSizeY, ercVariables_t *errorVar);
int ercConcealInterFrame(frame *recfr, objectBuffer_t *object_list, int picSizeX, int picSizeY, ercVariables_t *errorVar, int chroma_format_idc);


int ercCollect8PredBlocks(int predBlocks[], int currRow, int currColumn, char* condition,
                          int maxRow, int maxColumn, int step, uint8_t fNoCornerNeigh);

#define isSplitted(object_list,currMBNum) \
    ((object_list+((currMBNum)<<2))->regionMode >= REGMODE_SPLITTED)

/* this can be used as isBlock(...,INTRA) or isBlock(...,INTER_COPY) */
#define isBlock(object_list,currMBNum,comp,regMode) \
    (isSplitted(object_list,currMBNum) ? \
     ((object_list+((currMBNum)<<2)+(comp))->regionMode == REGMODE_##regMode##_8x8) : \
     ((object_list+((currMBNum)<<2))->regionMode == REGMODE_##regMode))

/* this can be used as getParam(...,mv) or getParam(...,xMin) or getParam(...,yMin) */
#define getParam(object_list,currMBNum,comp,param) \
    (isSplitted(object_list,currMBNum) ? \
     ((object_list+((currMBNum)<<2)+(comp))->param) : \
     ((object_list+((currMBNum)<<2))->param))

#endif
