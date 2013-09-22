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
 *  File      : neighbour.h
 *  Author(s) : Luuvish
 *  Version   : 1.0
 *  Revision  :
 *      1.0 June 16, 2013    first release
 *
 * ===========================================================================
 */

#ifndef _NEIGHBOUR_H_
#define _NEIGHBOUR_H_


struct position_t {
    int x;
    int y;
};

using pos_t = position_t;

struct location_t {
    int x;
    int y;
};

using loc_t = location_t;

loc_t operator + (const loc_t& l, const loc_t& r);
loc_t operator - (const loc_t& l, const loc_t& r);


void CheckAvailabilityOfNeighbors(mb_t* mb);
void CheckAvailabilityOfNeighborsCABAC(mb_t* mb);

void getNeighbour      (mb_t* mb, int xN, int yN, int mb_size[2], PixelPos *pix);
void get4x4Neighbour   (mb_t* mb, int xN, int yN, int mb_size[2], PixelPos *pix);

void get_neighbors     (mb_t* mb, PixelPos *block, int mb_x, int mb_y, int blockshape_x);
void check_dp_neighbors(mb_t* mb);


struct neighbour_t {
    loc_t get_location(slice_t* slice, int mbAddr, const pos_t& offset={0,0});
    mb_t* get_mb      (slice_t* slice, const loc_t& loc);
    int   get_mbaddr  (slice_t* slice, const loc_t& loc);
    pos_t get_blkpos  (slice_t* slice, const loc_t& loc);

    pos_t get_position(slice_t* slice, int mbAddr, int blkIdx=0);
    int   get_blkIdx  (slice_t* slice, const pos_t& pos);

    void  get_mb2pos(slice_t* slice, int mbAddr, int& xI, int& yI);
    void  get_mb2blk(slice_t* slice, int mbAddr, int& xpos, int& ypos);
    mb_t* get_pos2mb(slice_t* slice, int xP, int yP, int& mbAddr);
    int   get_pos2mb(slice_t* slice, int xP, int yP);
    int   predict_nnz(mb_t* mb, int pl, int i, int j);
};


#endif /* _NEIGHBOUR_H_ */
