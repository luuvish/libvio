/*
 * ===========================================================================
 *
 *   This confidential and proprietary software may be used only
 *  as authorized by a licensing agreement from Thumb o'Cat Inc.
 *  In the event of publication, the following notice is applicable:
 * 
 *       Copyright (C) 2013 - 2013 Thumb o'Cat
 *                     All right reserved.
 * 
 *   The entire notice above must be reproduced on all authorized copies.
 *
 * ===========================================================================
 *
 *  File      : deblock.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _DEBLOCK_H_
#define _DEBLOCK_H_


#include "global.h"
#include "dpb.h"

struct slice_t;

struct deblock_t {
	void init();
	void deblock(VideoParameters *p_Vid, storable_picture *p);

private:
	void init_neighbors(VideoParameters *p_Vid);
	void make_frame_picture_JV(VideoParameters *p_Vid);

	void deblock_pic(VideoParameters *p_Vid, storable_picture *p);
	void deblock_mb(macroblock_t* mb);

	void make_bS();

  	void edge_loop_luma_ver  (ColorPlane pl, uint8_t* Strength, macroblock_t* MbQ, int edge);
  	void edge_loop_luma_hor  (ColorPlane pl, uint8_t* Strength, macroblock_t* MbQ, int edge);
  	void edge_loop_chroma_ver(ColorPlane pl, uint8_t* Strength, macroblock_t* MbQ, int edge);
  	void edge_loop_chroma_hor(ColorPlane pl, uint8_t* Strength, macroblock_t* MbQ, int edge);

	void edge_loop(bool verticalEdgeFlag, bool chromaEdgeFlag, bool chromaStyleFilteringFlag,
                   macroblock_t* MbQ, uint8_t* Strength, ColorPlane pl, int edge);

	void deblock_strong(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag);
	void deblock_normal(imgpel *pixP, imgpel *pixQ, int widthP, int widthQ, int alpha, int beta, int bS, bool chromaStyleFilteringFlag, int chromaEdgeFlag, int BitDepthY, int BitDepthC, int indexA);
};


extern deblock_t deblock;


#endif /* _DEBLOCK_H_ */
