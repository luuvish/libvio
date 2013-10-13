#ifndef _ERC_API_H_
#define _ERC_API_H_


#include "defines.h"

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

//! YUV pixel domain image arrays for a video frame
typedef struct frame_s
{
  VideoParameters *p_Vid;
  imgpel *yptr;
  imgpel *uptr;
  imgpel *vptr;
} frame;

//! region structure stores information about a region that is needed for concealment
typedef struct object_buffer
{
  byte regionMode;  //!< region mode as above
  int xMin;         //!< X coordinate of the pixel position of the top-left corner of the region
  int yMin;         //!< Y coordinate of the pixel position of the top-left corner of the region
  int mv[3];        //!< motion vectors in 1/4 pixel units: mvx = mv[0], mvy = mv[1],
                    //!< and ref_frame = mv[2]
} objectBuffer_t;

/*
* Defines
*/

/* If the average motion vector of the correctly received macroblocks is less than the
threshold, concealByCopy is used, otherwise concealByTrial is used. */
#define MVPERMB_THR 8

/* used to determine the size of the allocated memory for a temporal Region (MB) */
#define DEF_REGION_SIZE 384  /* 8*8*6 */

#define ERC_BLOCK_OK                3
#define ERC_BLOCK_CONCEALED         2
#define ERC_BLOCK_CORRUPTED         1
#define ERC_BLOCK_EMPTY             0


/*
* Functions to convert MBNum representation to blockNum
*/

#define xPosYBlock(currYBlockNum,picSizeX) \
((currYBlockNum)%((picSizeX)>>3))

#define yPosYBlock(currYBlockNum,picSizeX) \
((currYBlockNum)/((picSizeX)>>3))

#define xPosMB(currMBNum,picSizeX) \
((currMBNum)%((picSizeX)>>4))

#define yPosMB(currMBNum,picSizeX) \
((currMBNum)/((picSizeX)>>4))

#define MBxy2YBlock(currXPos,currYPos,comp,picSizeX) \
((((currYPos)<<1)+((comp)>>1))*((picSizeX)>>3)+((currXPos)<<1)+((comp)&1))

#define MBNum2YBlock(currMBNum,comp,picSizeX) \
MBxy2YBlock(xPosMB((currMBNum),(picSizeX)),yPosMB((currMBNum),(picSizeX)),(comp),(picSizeX))


/*
* typedefs
*/

/* segment data structure */
typedef struct ercSegment_s
{
  short     startMBPos;
  short     endMBPos;
  char      fCorrupted;
} ercSegment_t;

/* Error detector & concealment instance data structure */
typedef struct ercVariables_s
{
  /*  Number of macroblocks (size or size/4 of the arrays) */
  int   nOfMBs;
  /* Number of segments (slices) in frame */
  int     nOfSegments;

  /*  Array for conditions of Y blocks */
  char     *yCondition;
  /*  Array for conditions of U blocks */
  char     *uCondition;
  /*  Array for conditions of V blocks */
  char     *vCondition;

  /* Array for slice_t level information */
  ercSegment_t *segments;
  int     currSegment;

  /* Conditions of the MBs of the previous frame */
  char   *prevFrameYCondition;

  /* Flag telling if the current segment was found to be corrupted */
  int   currSegmentCorrupted;
  /* Counter for corrupted segments per picture */
  int   nOfCorruptedSegments;

  /* State variables for error detector and concealer */
  int   concealment;

} ercVariables_t;

/*
* External function interface
*/

void ercInit (VideoParameters *p_Vid, int pic_sizex, int pic_sizey, int flag);
ercVariables_t *ercOpen( void );
void ercReset( ercVariables_t *errorVar, int nOfMBs, int numOfSegments, int picSizeX );
void ercClose( VideoParameters *p_Vid, ercVariables_t *errorVar );
void ercSetErrorConcealment( ercVariables_t *errorVar, int value );

void ercStartSegment( int currMBNum, int segment, unsigned int bitPos, ercVariables_t *errorVar );
void ercStopSegment( int currMBNum, int segment, unsigned int bitPos, ercVariables_t *errorVar );
void ercMarkCurrSegmentLost(int picSizeX, ercVariables_t *errorVar );
void ercMarkCurrSegmentOK(int picSizeX, ercVariables_t *errorVar );
void ercMarkCurrMBConcealed( int currMBNum, int comp, int picSizeX, ercVariables_t *errorVar );

int ercConcealIntraFrame( VideoParameters *p_Vid, frame *recfr, int picSizeX, int picSizeY, ercVariables_t *errorVar );
int ercConcealInterFrame( frame *recfr, objectBuffer_t *object_list,
                          int picSizeX, int picSizeY, ercVariables_t *errorVar, int chroma_format_idc );


/* Thomson APIs for concealing entire frame loss */

#include "dpb.h"
#include "output.h"

struct concealment_node {
    storable_picture* picture;
    int missingpocs;
    concealment_node* next;
};

extern struct concealment_node * init_node(storable_picture* , int );
extern void init_lists_for_non_reference_loss(dpb_t *p_Dpb, int , bool );

extern void conceal_lost_frames (dpb_t *p_Dpb, struct slice_t *pSlice);

extern int comp(const void *, const void *);

void erc_picture(VideoParameters* p_Vid);
void ercWriteMBMODEandMV(mb_t *currMB);


#endif
