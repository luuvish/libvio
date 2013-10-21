#ifndef _ERC_API_H_
#define _ERC_API_H_


#include "defines.h"

struct VideoParameters;
struct decoded_picture_buffer_t;
using dpb_t = decoded_picture_buffer_t;
struct slice_t;
namespace vio { namespace h264 {
struct macroblock_t;
}}
using mb_t = vio::h264::macroblock_t;

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
};

/*
* External function interface
*/

extern void ercInit(VideoParameters *p_Vid, int pic_sizex, int pic_sizey, int flag);
extern void ercReset(ercVariables_t *errorVar, int nOfMBs, int numOfSegments, int picSizeX);
extern void ercClose(VideoParameters *p_Vid, ercVariables_t *errorVar);

extern void erc_picture(VideoParameters* p_Vid);
extern void ercWriteMBMODEandMV(mb_t *currMB);


/* Thomson APIs for concealing entire frame loss */

struct concealment_node {
    storable_picture* picture;
    int missingpocs;
    concealment_node* next;
};


#endif
