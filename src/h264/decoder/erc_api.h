#ifndef _ERC_API_H_
#define _ERC_API_H_


struct VideoParameters;
namespace vio { namespace h264 {
struct macroblock_t;
}}
using mb_t = vio::h264::macroblock_t;


struct ercSegment_t;
struct objectBuffer_t;
struct frame;

/* Error detector & concealment instance data structure */
struct ercVariables_t {
    int         nOfMBs;               /* Number of macroblocks (size or size/4 of the arrays) */
    int         nOfSegments;          /* Number of segments (slices) in frame */
    char*       yCondition;           /* Array for conditions of Y blocks */
    char*       uCondition;           /* Array for conditions of U blocks */
    char*       vCondition;           /* Array for conditions of V blocks */
    ercSegment_t* segments;           /* Array for slice_t level information */
    int         currSegment;
    char*       prevFrameYCondition;  /* Conditions of the MBs of the previous frame */
    int         currSegmentCorrupted; /* Flag telling if the current segment was found to be corrupted */
    int         nOfCorruptedSegments; /* Counter for corrupted segments per picture */
    int         concealment;          /* State variables for error detector and concealer */

    objectBuffer_t* erc_object_list;

    int         erc_mvperMB;

    ercVariables_t(int pic_sizex, int pic_sizey, int flag);
    ~ercVariables_t();

    void ercReset(int nOfMBs, int numOfSegments, int picSizeX);
    void ercWriteMBMODEandMV(mb_t* currMB, storable_picture* pic);
    void erc_picture(storable_picture* pic);

private:
    void ercStartSegment(int currMBNum, int segment);
    void ercStopSegment (int currMBNum, int segment);
    void ercMarkCurrSegmentLost(int picSizeX);
    void ercMarkCurrSegmentOK  (int picSizeX);

    int  ercConcealIntraFrame(frame* recfr, int picSizeX, int picSizeY);
    int  ercConcealInterFrame(frame* recfr, int picSizeX, int picSizeY, int chroma_format_idc);

    int  ercCollect8PredBlocks(int predBlocks[], int currRow, int currColumn, char* condition,
                               int maxRow, int maxColumn, int step, uint8_t fNoCornerNeigh);

    void pixMeanInterpolateBlock(px_t* src[], px_t* block, int blockSize, int frameWidth, int BitDepth);
    void ercPixConcealIMB       (px_t* currFrame, int row, int column, int predBlocks[], int frameWidth, int mbWidthInBlocks, int BitDepth);
    int  ercCollectColumnBlocks (int predBlocks[], int currRow, int currColumn, char* condition, int maxRow, int maxColumn, int step);
    void concealBlocks          (int lastColumn, int lastRow, int comp, frame* recfr, int picSizeX, char* condition);

    void buildPredRegionYUV    (VideoParameters* p_Vid, int* mv, int x, int y, px_t* predMB);
    int  edgeDistortion        (int predBlocks[], int currYBlockNum, px_t* predMB, px_t* recY, int picSizeX, int regionSize);
    void copyPredMB            (int currYBlockNum, px_t* predMB, frame* recfr, int picSizeX, int regionSize);
    int  concealByTrial        (frame* recfr, px_t* predMB, int currMBNum, int predBlocks[], int picSizeX, int picSizeY);
    void copyBetweenFrames     (frame* recfr, int currYBlockNum, int picSizeX, int regionSize);
    int  concealByCopy         (frame* recfr, int currMBNum, int picSizeX);
    void ercMarkCurrMBConcealed(int currMBNum, int comp, int picSizeX);
};


#endif
