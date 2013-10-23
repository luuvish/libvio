#ifndef _ERC_API_H_
#define _ERC_API_H_


struct VideoParameters;
namespace vio { namespace h264 {
struct macroblock_t;
}}
using mb_t = vio::h264::macroblock_t;

struct storable_picture;

struct ercSegment_t;
struct objectBuffer_t;

/* Error detector & concealment instance data structure */
struct ercVariables_t {
    int         nOfMBs;               /* Number of macroblocks (size or size/4 of the arrays) */
    int         nOfSegments;          /* Number of segments (slices) in frame */
    char*       yCondition;           /* Array for conditions of Y blocks */
    char*       uCondition;           /* Array for conditions of U blocks */
    char*       vCondition;           /* Array for conditions of V blocks */
    ercSegment_t* segments;           /* Array for slice_t level information */
    int         currSegment;
    int         currSegmentCorrupted; /* Flag telling if the current segment was found to be corrupted */
    int         nOfCorruptedSegments; /* Counter for corrupted segments per picture */
    bool        concealment;          /* State variables for error detector and concealer */

    objectBuffer_t* erc_object_list;

    int         erc_mvperMB;

    ercVariables_t(int pic_sizex, int pic_sizey, bool flag);
    ~ercVariables_t();

    void reset(int nOfMBs, int numOfSegments);
    void ercWriteMBMODEandMV(mb_t& mb, uint8_t slice_type, storable_picture* pic);
    void erc_picture(storable_picture* pic);

private:
    void ercStartSegment(int currMBNum, int segment);
    void ercStopSegment (int currMBNum, int segment);
    void ercMarkCurrSegmentLost(int picSizeX);
    void ercMarkCurrSegmentOK  (int picSizeX);

    int  ercConcealIntraFrame(storable_picture* pic);
    int  ercConcealInterFrame(storable_picture* pic);

    int  ercCollect8PredBlocks(int predBlocks[], int currRow, int currColumn, char* condition,
                               int maxRow, int maxColumn, int step, uint8_t fNoCornerNeigh);

    int  ercCollectColumnBlocks (int predBlocks[], int currRow, int currColumn, char* condition, int maxRow, int maxColumn, int step);
    void pixMeanInterpolateBlock(px_t* src[], px_t* block, int blockSize, int frameWidth, int BitDepth);
    void ercPixConcealIMB       (storable_picture* pic, int comp, int row, int column, int predBlocks[]);
    void concealBlocks          (storable_picture* pic, int comp);

    void buildPredRegionYUV    (storable_picture* pic, int* mv, int x, int y, px_t* predMB);
    int  edgeDistortion        (storable_picture* pic, int predBlocks[], int currYBlockNum, px_t* predMB, int regionSize);
    void copyPredMB            (storable_picture* pic, int currYBlockNum, px_t* predMB, int regionSize);
    int  concealByTrial        (storable_picture* pic, px_t* predMB, int currMBNum, int predBlocks[]);
    void copyBetweenFrames     (storable_picture* dec_pic, storable_picture* ref_pic, int currYBlockNum, int regionSize);
    int  concealByCopy         (storable_picture* dec_pic, storable_picture* ref_pic, int currMBNum);
    void ercMarkCurrMBConcealed(int currMBNum, int comp, int picSizeX);
};


#endif
