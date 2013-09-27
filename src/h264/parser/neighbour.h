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


struct neighbour_t {
    loc_t get_location(slice_t* slice, int mbAddr, const pos_t& offset={0,0});
    loc_t get_location_c(slice_t* slice, int mbAddr, const pos_t& offset={0,0});
    mb_t* get_mb      (slice_t* slice, const loc_t& loc);
    mb_t* get_mb_c    (slice_t* slice, const loc_t& loc);
    int   get_mbaddr  (slice_t* slice, const loc_t& loc);
    pos_t get_blkpos  (slice_t* slice, const loc_t& loc);
    pos_t get_blkpos_c(slice_t* slice, const loc_t& loc);

    pos_t get_position(slice_t* slice, int mbAddr, int blkIdx);
    mb_t* get_mb      (slice_t* slice, int mbAddr, const pos_t& offset={0,0});
    pos_t get_blkpos  (slice_t* slice, int mbAddr, const pos_t& offset={0,0});

    int   predict_nnz(mb_t* mb, int pl, int i, int j);
};


#endif /* _NEIGHBOUR_H_ */
